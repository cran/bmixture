
## Main function: BDMCMC algorithm for finite mixture of t-distribution
################################################################################
# INPUT for bdmcmc funciton 
# 1) data:         the data with posetive and no missing values
# 2) k             number of components of mixture distribution. Defult is unknown
# 2) iter:         nuber of iteration of the algorithm 
# 3) burnin:       number of burn-in iteration
# 4) lambda_r:       rate for birth and parameter of prior distribution of k
# 7) k, mu, sig, and pa: initial values for parameters respectively k, mu, sig and pi
################################################################################
bmixt = function( data, k = "unknown", iter = 1000, burnin = iter / 2, lambda = 1, 
				  df = 1,
                  k.start = NULL, mu.start = NULL, sig.start = NULL, pi.start = NULL, 
                  k_max = 30, trace = TRUE )
{
	if( any( is.na( data ) ) ) stop( "Data should contain no missing data" ) 
	if( iter <= burnin )       stop( "Number of iteration must be more than number of burn-in" )	

	burnin   = floor( burnin )
	n        = length( data )
	lambda_r = lambda
	df_t     = df
	
	max_data = max( data )	
	min_data = min( data )
	R        = max_data - min_data
	
	# Values for paprameters of prior distributon of mu
	epsilon = R / 2   # midpoint of the observed range of the data
	kappa_r = 1 / R ^ 2
	# Values for paprameters of prior distributon of sigma
	alpha = 2
	# Values for paprameters of prior distributon of Beta ( Beta is a hyperparameter )
	g = 0.2
	h = ( 100 * g ) / ( alpha * R ^ 2 )
		
	# initial value
	if( k == "unknown" )
	{
		component_size = "unknown"
		if( is.null( k.start ) ) { k = 3 } else { k = k.start }
	}else{
		component_size = "fixed"
	}

	beta_r = rgamma( 1, g, h )
	if( is.null( pi.start  ) ) pi.start  = c( rep( 1 / k, k  ) )
	if( is.null( mu.start  ) ) mu.start  = rnorm( k, epsilon, sqrt( 1 / kappa_r ) ) 
	if( is.null( sig.start ) ) sig.start = 1 / rgamma( k, alpha, beta_r ) 	

	pi_r = pi.start
	mu   = mu.start
    sig  = sig.start
    q_t  = rgamma( n, df_t / 2, df_t / 2 )

	# Sort parameters based on mu
	order_pi = order( pi_r )
	pi_r     = pi_r[order_pi]
	mu       = mu[order_pi]
	sig      = sig[order_pi]
    
############### MCMC 
	if( component_size == "unknown" )
	{
		pi_sample  = matrix( 0, nrow = iter - burnin, ncol = k_max ) 
		mu_sample  = pi_sample
		sig_sample = pi_sample
		all_k       = vector( mode = "numeric", length = iter )
		all_weights = all_k

		data_r  = data
		k_max_r = k_max
		q_t_r   = q_t
		df_t_r  = df_t
		
		result = .C( "bmix_t_unknown_k", as.double(data_r), as.integer(n), as.integer(k), as.integer(k_max_r), as.integer(iter), as.integer(burnin), as.double(lambda_r),
 						pi_sample = as.double(pi_sample), mu_sample = as.double(mu_sample), sig_sample = as.double(sig_sample),
 						all_k = as.integer(all_k), all_weights = as.double(all_weights),
						as.double(epsilon), as.double(kappa_r), as.double(alpha), as.double(beta_r), as.double(g), as.double(h),
						as.double(mu), as.double(sig), as.double(pi_r), 
						as.double(q_t_r), as.integer (df_t_r), PACKAGE = "bmixture" )
		
		all_k       = result $ all_k
		all_weights = result $ all_weights

		pi_sample  = matrix( result $ pi_sample , nrow = iter - burnin, ncol = k_max_r )
		mu_sample  = matrix( result $ mu_sample , nrow = iter - burnin, ncol = k_max_r )
		sig_sample = matrix( result $ sig_sample, nrow = iter - burnin, ncol = k_max_r )
		   		
		mcmc_sample = list( all_k = all_k, all_weights = all_weights, pi_sample = pi_sample, mu_sample = mu_sample, sig_sample = sig_sample, data = data_r, df_t = df_t, component_size = "unknown" )    
	} else {
		pi_sample  = matrix( 0, nrow = iter - burnin, ncol = k ) 
		mu_sample  = pi_sample
		sig_sample = pi_sample
		
		data_r = data
		q_t_r  = q_t
		df_t_r = df_t
      
		result = .C( "bmix_t_fixed_k", as.double(data_r), as.integer(n), as.integer(k), as.integer(iter), as.integer(burnin), 
						pi_sample = as.double(pi_sample), mu_sample = as.double(mu_sample), sig_sample = as.double(sig_sample),
						as.double(epsilon), as.double(kappa_r), as.double(alpha), as.double(g), as.double(h),
						as.double(mu), as.double(sig), as.double(pi_r),  
						as.double(q_t_r), as.integer (df_t_r) #)
						, PACKAGE = "bmixture" )

		pi_sample  = matrix( result $ pi_sample , nrow = iter - burnin, ncol = k )
		mu_sample  = matrix( result $ mu_sample , nrow = iter - burnin, ncol = k )
		sig_sample = matrix( result $ sig_sample, nrow = iter - burnin, ncol = k )

		mcmc_sample = list( pi_sample = pi_sample, mu_sample = mu_sample, sig_sample = sig_sample, data = data, df_t = df_t, component_size = "fixed" )    
	}

	if( trace == TRUE )
	{
		mes <- paste( c(" ", iter," iteration done.                               " ), collapse = "" )
		cat( mes, "\r" )
		cat( "\n" )
		flush.console()
	}    

	class( mcmc_sample ) = "bmixt"
	return( mcmc_sample )
}
   
# summary of bmixt output
summary.bmixt = function( object, ... )
{
	component_size = object $ component_size
	pi_sample      = object $ pi_sample
	mu_sample      = object $ mu_sample
	sig_sample     = object $ sig_sample
	data           = object $ data
	df_t           = object $ df_t
	
	sample_size    = nrow( sig_sample )
   
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
     
# plot for class bmixt
plot.bmixt = function( x, ... )
{
	component_size = x $ component_size
	pi_sample      = x $ pi_sample
	mu_sample      = x $ mu_sample
	sig_sample     = x $ sig_sample
	data           = x $ data
	df_t           = x $ df_t
	
	sample_size    = nrow( sig_sample )

	# plot for estimated distribution
	hist( data, prob = T, nclass = 25, col = "gray", border = "white"  )
	
	x_seq       <- seq( min( data ) * 0.9, max( data ) * 1.2, length = 500 )
	f_hat_x_seq <- 0 * x_seq
	size_x_seq_r = length( x_seq )	
	
	if( component_size == "unknown" )
	{
		all_k       = x $ all_k
		iter        = length( all_k )
		burnin      = iter - sample_size
		k           = all_k[( burnin + 1 ):iter]

		result = .C( "dmixt_hat_x_seq_unknow_k", as.double(x_seq), f_hat_x_seq = as.double(f_hat_x_seq), as.integer(df_t), 
					 as.double(pi_sample), as.double(mu_sample), as.double(sig_sample),
					 as.integer(k), as.integer(sample_size), as.integer(size_x_seq_r)
					 , PACKAGE = "bmixture" )
					 
		f_hat_x_seq = result $ f_hat_x_seq
	
		lines( x_seq, f_hat_x_seq / sample_size, col = "black", lty = 2, lw = 1 )
	}
	else
	{
		size_mix = ncol( sig_sample )

		result = .C( "dmixt_hat_x_seq_fixed_k", as.double(x_seq), f_hat_x_seq = as.double(f_hat_x_seq), as.integer(df_t), 
					 as.double(pi_sample), as.double(mu_sample), as.double(sig_sample),
					 as.integer(size_mix), as.integer(sample_size), as.integer(size_x_seq_r)
					 , PACKAGE = "bmixture" )
					 
		f_hat_x_seq = result $ f_hat_x_seq

		lines( x_seq, f_hat_x_seq / sample_size, col = "black", lty = 2, lw = 1 )
	}
	
    legend( "topright", c( "predictive density" ), lty = 2, col = "black", lwd = 1 )
}
     
# print of the bmixt output
print.bmixt = function( x, ... )
{
	component_size = x $ component_size
	pi_sample      = x $ pi_sample
	mu_sample      = x $ mu_sample
	sig_sample     = x $ sig_sample
	df_t           = x $ df_t

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
		cat( paste( "Number of mixture components = ", ncol( mu_sample ) ), fill = TRUE ) 
		cat( paste( "Estimated pi  = "), paste( round( apply( pi_sample , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated mu  = "), paste( round( apply( mu_sample , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated sig = "), paste( round( apply( sig_sample, 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated df  = " ), df_t                                            , fill = TRUE ) 
		cat( paste( "" ), fill = TRUE )				
	}
} 
   






