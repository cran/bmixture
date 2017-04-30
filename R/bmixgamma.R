
## Main function: BDMCMC algorithm for finite mixture of Gamma distribution
################################################################################
# INPUT for bdmcmc funciton 
# 1) data:         the data with posetive and no missing values
# 2) k             number of components of mixture distribution. Defult is unknown
# 2) iter:         nuber of iteration of the algorithm 
# 3) burnin:       number of burn-in iteration
# 4) lambda:       rate for birth and parameter of prior distribution of k
# 5) mu and nu:    alpha parameters of alpha in mixture distribution
# 6) kesi and tau: alpha parameters of alpha in mixture distribution 
# 7) k, alpha, beta, and pa: initial values for parameters respectively k, alpha, beta and pi_r
################################################################################
bmixgamma = function( data, k = "unknown", iter = 1000, burnin = iter / 2, lambda = 1, 
                    mu = NULL, nu = NULL, kesi = NULL, tau = NULL, 
                    k.start = NULL, alpha.start = NULL, beta.start = NULL, pi.start = NULL, 
                    k_max = 30, trace = TRUE )
{
	if( any( is.na( data ) ) ) stop( "Data should contain no missing data" ) 
	if( any( data < 0  ) )     stop( "Data should have positive value" ) 
	if( iter <= burnin )       stop( "Number of iteration must be more than number of burn-in" )	

	burnin   = floor( burnin )
	n        = length( data )
	lambda_r = lambda
	
	# Values for paprameters of prior distributon of alpha
	if( is.null( mu ) ) mu <- ( mean(data) ^ 2 / var( data ) ) ^ ( 2 / 3 ) 
	if( is.null( nu ) ) nu <- 1 / sqrt( mu )
	
	# Values for paprameters of prior distributon of beta
	if( is.null( kesi ) ) kesi <- ( mean( data ) / var( data ) ) ^ ( 2 / 3 )
	if( is.null( tau ) )  tau  <- 1 / sqrt( kesi ) 
		
	# initial value
	if( k == "unknown" )
	{
		component_size = "unknown"
		if( is.null( k.start ) ) { k = 3 } else { k = k.start }
	}else{
		component_size = "fixed"
	}

	if( is.null( pi.start ) )    pi.start    <- c( rep( 1 / k, k  ) )
	if( is.null( alpha.start ) ) alpha.start <- rgamma( k, mean( data ) ^ 2, var( data ) ) # moment estimation of alpha
	if( is.null( beta.start  ) ) beta.start  <- rgamma( k, mean( data )    , var( data ) ) # moment estimation of data	

	pi_r  = pi.start
	alpha = alpha.start
    beta  = beta.start

	# Sort parameters based on pi_r
	order_pi = order( pi_r )
	pi_r     = pi_r[  order_pi ]
	alpha    = alpha[ order_pi ]
	beta     = beta[  order_pi ]
    
############### MCMC 
	if( component_size == "unknown" )
	{
		pi_sample    = matrix( 0, nrow = iter - burnin, ncol = k_max ) 
		alpha_sample = pi_sample
		beta_sample  = pi_sample
		all_k        = c( rep( 0, iter ) )
		all_weights  = all_k

		data_r  = data
		k_max_r = k_max
		
 		result = .C( "bmix_gamma_unknown_k", as.double(data_r), as.integer(n), as.integer(k), as.integer(k_max_r), as.integer(iter), as.integer(burnin), as.double(lambda_r),
  						pi_sample = as.double(pi_sample), alpha_sample = as.double(alpha_sample), beta_sample = as.double(beta_sample),
  						all_k = as.integer(all_k), all_weights = as.double(all_weights),
 						as.double(mu), as.double(nu), as.double(kesi), as.double(tau),
 						as.double(alpha), as.double(beta), as.double(pi_r) #)
 						, PACKAGE = "bmixture" )
 		
 		all_k       = result $ all_k
 		all_weights = result $ all_weights
 
 		pi_sample    = matrix( result $ pi_sample ,   nrow = iter - burnin, ncol = k_max_r )
 		alpha_sample = matrix( result $ alpha_sample, nrow = iter - burnin, ncol = k_max_r )
 		beta_sample  = matrix( result $ beta_sample,  nrow = iter - burnin, ncol = k_max_r )
 		   		
 		mcmc_sample = list( all_k = all_k, all_weights = all_weights, pi_sample = pi_sample, alpha_sample = alpha_sample, beta_sample = beta_sample, data = data_r, component_size = "unknown" )    
	} else {
		pi_sample    = matrix( 0, nrow = iter - burnin, ncol = k ) 
		alpha_sample = pi_sample
		beta_sample  = pi_sample
		
		data_r = data
  
 		result = .C( "bmix_gamma_fixed_k", as.double(data_r), as.integer(n), as.integer(k), as.integer(iter), as.integer(burnin), 
 						pi_sample = as.double(pi_sample), alpha_sample = as.double(alpha_sample), beta_sample = as.double(beta_sample),
 						as.double(mu), as.double(nu), as.double(kesi), as.double(tau),
 						as.double(alpha), as.double(beta), as.double(pi_r) #)
 						, PACKAGE = "bmixture" )
 
 		pi_sample    = matrix( result $ pi_sample   , nrow = iter - burnin, ncol = k )
 		alpha_sample = matrix( result $ alpha_sample, nrow = iter - burnin, ncol = k )
 		beta_sample  = matrix( result $ beta_sample , nrow = iter - burnin, ncol = k )
 
 		mcmc_sample = list( pi_sample = pi_sample, alpha_sample = alpha_sample, beta_sample = beta_sample, data = data, component_size = "fixed" )    
	}

	if( trace == TRUE )
	{
		mes <- paste( c(" ", iter," iteration done.                               " ), collapse = "" )
		cat( mes, "\r" )
		cat( "\n" )
		flush.console()
	}    

	class( mcmc_sample ) = "bmixgamma"
	return( mcmc_sample )
}
   
# summary of bmixgamma output
summary.bmixgamma = function( object, ... )
{
	component_size = object $ component_size
	pi_sample      = object $ pi_sample
	alpha_sample   = object $ alpha_sample
	beta_sample    = object $ beta_sample
	data           = object $ data
	
	sample_size    = nrow( beta_sample )
   
	if( component_size == "unknown" )
	{
		all_k       = object $ all_k
		all_weights = object $ all_weights
		iter        = length( all_k )
		burnin      = iter - sample_size
		k           = all_k[( burnin + 1 ):iter]
		weights     = all_weights[( burnin + 1 ):iter]
		
		op = par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) ) 
	}

	plot( object )
	
	if( component_size == "unknown" )
	{
		# plot estimated distribution of k
		max_k = max( k )
		y     = vector( mode = "numeric", length = max_k )
		for ( i in 1:max_k ) y[i] <- sum( weights[k == i] ) 
		
		plot( x = 1:max_k, y, type = "h", main = "", ylab = "Pr(k|data)", xlab = "k(Number of components)" )

		# plot k based on iterations
		sum_weights = 0 * weights
		sum_weights[1] = weights[1]
		for( i in 2:length( k ) ) sum_weights[i] = sum_weights[i - 1] + weights[i]

		plot( sum_weights, k, type = "l", xlab = "iteration", ylab = "Number of componants" )
	}
	else
	{
		print( object )
	}
}  
     
# plot for class bmixgamma
plot.bmixgamma = function( x, ... )
{
	component_size = x $ component_size
	pi_sample      = x $ pi_sample
	alpha_sample   = x $ alpha_sample
	beta_sample    = x $ beta_sample
	data           = x $ data
	
	sample_size    = nrow( beta_sample )

	# plot for estimated distribution
	hist( data, prob = T, nclass = 25, col = "gray", border = "white"  )
	
	x_seq       <- seq( min( data ) * 0.9, max( data ) * 1.2, length = 500 )
	f_hat_x_seq <- 0 * x_seq
	size_x_seq_r = length( x_seq )	
	
	if( component_size == "unknown" )
	{
		all_k  = x $ all_k
		iter   = length( all_k )
		burnin = iter - sample_size
		k      = all_k[( burnin + 1 ):iter]

		result = .C( "dmixgamma_hat_x_seq_unknow_k", as.double(x_seq), f_hat_x_seq = as.double(f_hat_x_seq), 
					 as.double(pi_sample), as.double(alpha_sample), as.double(beta_sample),
					 as.integer(k), as.integer(sample_size), as.integer(size_x_seq_r)
					 , PACKAGE = "bmixture" )
					 
		f_hat_x_seq = result $ f_hat_x_seq
	
		lines( x_seq, f_hat_x_seq / sample_size, col = "black", lty = 2, lw = 1 )
	}
	else
	{
		size_mix = ncol( beta_sample )

		result = .C( "dmixgamma_hat_x_seq_fixed_k", as.double(x_seq), f_hat_x_seq = as.double(f_hat_x_seq), 
					 as.double(pi_sample), as.double(alpha_sample), as.double(beta_sample),
					 as.integer(size_mix), as.integer(sample_size), as.integer(size_x_seq_r)
					 , PACKAGE = "bmixture" )
					 
		f_hat_x_seq = result $ f_hat_x_seq

		lines( x_seq, f_hat_x_seq / sample_size, col = "black", lty = 2, lw = 1 )
	}
	
    legend( "topright", c( "predictive density" ), lty = 2, col = "black", lwd = 1 )
}
     
# print of the bmixgamma output
print.bmixgamma = function( x, ... )
{
	component_size = x $ component_size
	pi_sample      = x $ pi_sample
	alpha_sample   = x $ alpha_sample
	beta_sample    = x $ beta_sample

	if( component_size == "unknown" )
	{
		all_k       = x $ all_k
		all_weights = x $ all_weights
		iter        = length( all_k )
		sample_size = nrow( pi_sample )
		burnin      = iter - sample_size
		k           = all_k[( burnin + 1 ):iter]
		weights     = all_weights[( burnin + 1 ):iter]

		max_k = max( k )
		y     = vector( mode = "numeric", length = max_k )
		for ( i in 1:max_k ) y[i] <- sum( weights[k == i] ) 

		cat( paste( "" ), fill = TRUE )
		cat( paste( "Estimated of number of components = ", which( y == max( y ) ) ), fill = TRUE ) 
		cat( paste( "" ), fill = TRUE )
	} else {
		cat( paste( "" ), fill = TRUE )
		cat( paste( "Number of mixture components = ", ncol( alpha_sample ) ), fill = TRUE ) 
		cat( paste( "Estimated pi    = "), paste( round( apply( pi_sample , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated alpha = "), paste( round( apply( alpha_sample , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated beta  = "), paste( round( apply( beta_sample, 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "" ), fill = TRUE )				
	}
} 
   


















