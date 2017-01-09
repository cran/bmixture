
# MCMC sampling algorithm based on Birth-Death MCMC scheme
# for mixture of Gamma distribution with an unknown number of components  
bmixgamma_unknown_k = function( data, n, k, iter, burnin, lambda, 
							mu, nu, kesi, tau, 
							alpha, beta, pi, 
							k_max, trace )
{
	all_k        = vector( mode = "numeric", length = iter )
	all_weights  = all_k
	pi_sample    = list()  
    alpha_sample = pi_sample
    beta_sample  = pi_sample

	counter = 1
    for( i_mcmc in 1:iter )
    {
		if ( trace == TRUE && i_mcmc %% 10 == 0 )
		{
			mes <- paste( c( " iteration: ", i_mcmc, " from ", iter, ". Number of componets = ", k ), collapse = "" )
			cat( mes, "\r" )
			flush.console()	
		}    
    
		# first we obtain death_rates (vector of death rates for each component)
		if ( k != 1 ) 
		{ 
			death_rates = vector( mode = "numeric", length = k )
			
			for ( j in 1:k )
			{
				log_death_rates <- 0
				for ( t in 1:n )
					log_death_rates = log_death_rates + log( ( 1 - ( pi[j] * dgamma( data[t], alpha[j], beta[j] ) / sum( pi * dgamma( data[t], alpha, beta ) ) ) ) ) - log( 1 - pi[j] )
				
				death_rates[j] <- exp( log_death_rates )
			}
			
			# checking for infinite value
			death_rates[is.infinite( death_rates )] <- gamma( 171 )
		}			
		else 
		{
			death_rates = 0
		}
		
		# simulate the time s to the next jump, from an exponential distribution
		weight = 1 / ( lambda + sum( death_rates ) )
		# first Birth and death run and modify the parameters
		birthp <- lambda * weight
		
		if( ( runif( 1 ) < birthp ) && ( k < k_max ) )
		{
			k      <- k + 1
			pi_new <- rbeta( 1, 1, k - 1 )
			pi     <- c( pi * ( 1 - pi_new ), pi_new )
			alpha  <- c( alpha, rgamma( 1, mu  , nu  ) )
			beta   <- c( beta , rgamma( 1, kesi, tau ) )
		} 
		else 
		{
			j     <- sample( 1:k, 1, prob = death_rates )
			k     <- k - 1
			pi    <- pi[-j] / ( 1 - pi[j] )
			alpha <- alpha[-j]
			beta  <- beta[-j]
		}
		
		# step 2: updating z
		z <- matrix( NA, k, n )
		for ( i in 1:n )
			z[,i] = rmultinom( 1, size = 1, prob = pi * dgamma( data[i], alpha, beta ) )
   
		# step 3: updating pi
		sum_z = apply( z, 1, sum )
		pi    = rdirichlet( 1, 1 + sum_z )
				
		# step 4 : updating beta and alpha using a Metropolis-Hastings
		for ( i in 1:k )
		{
			# step 4.1 : updating beta
			beta[i] <- rgamma( 1, kesi + sum_z[i] * alpha[i], tau + sum( data * z[i,] ) ) 
			
			# step 4.2 : updating alpha using a Metropolis-Hastings
			pr             <- ( ( beta[i] * data ) * z[i,] )
			proposed_alpha <- rgamma( 1, mu, nu )
			accept_alpha   <- exp( sum_z[i] * ( lgamma( alpha[i] ) - lgamma( proposed_alpha ) ) + ( proposed_alpha - alpha[i] ) * sum( log( pr[pr > 0] ) ) )
			
			if ( accept_alpha > runif( 1 ) ) alpha[i] <- proposed_alpha 
		}
		
		all_k[i_mcmc]       = k # saving all values of k for chicking convergency
		all_weights[i_mcmc] = weight
		
		if ( i_mcmc > burnin )
		{
#			k_sample                             <- c( k_sample, k )
			pi_sample[[    counter ]] <- pi
			alpha_sample[[ counter ]] <- alpha
			beta_sample[[  counter ]] <- beta
			counter = counter + 1
		}
    }
    
	mcmc_sample = list( all_k = all_k, all_weights = all_weights, pi_sample = pi_sample, alpha_sample = alpha_sample, beta_sample = beta_sample, data = data, component_size = "unknown" )	
	return( mcmc_sample )
}
   
# MCMC sampling algorithm 
# for mixture of Gamma distribution with an fixed number of components  
bmixgamma_fixed_k = function( data, n, k, iter, burnin, 
							mu, nu, kesi, tau, 
							alpha, beta, pi, 
							trace )
{
	pi_sample    = matrix( NA, nrow = iter - burnin, ncol = k ) 
    alpha_sample = pi_sample
    beta_sample  = pi_sample

	counter = 1
    for( i_mcmc in 1:iter )
    {
		if( trace == TRUE && i_mcmc %% 10 == 0 )
		{
			mes <- paste( c( " iteration: ", i_mcmc, " from ", iter, ". Number of componets = ", k ), collapse = "" )
			cat( mes, "\r" )
			flush.console()	
		}    
    		
		# step 1: updating z
		z <- matrix( NA, k, n )
		for ( i in 1:n )
			z[,i] = rmultinom( 1, size = 1, prob = pi * dgamma( data[i], alpha, beta ) )
   
		# step 2: updating pi
		sum_z = apply( z, 1, sum )
		pi    = rdirichlet( 1, 1 + sum_z )
				
		# step 3 : updating beta and alpha using a Metropolis-Hastings
		for ( i in 1:k )
		{
			# step 3.1 : updating beta
			beta[i] <- rgamma( 1, kesi + sum_z[i] * alpha[i], tau + sum( data * z[i,] ) ) 
			
			# step 3.2 : updating alpha using a Metropolis-Hastings
			pr             <- ( ( beta[i] * data ) * z[i,] )
			proposed_alpha <- rgamma( 1, mu, nu )
			accept_alpha   <- exp( sum_z[i] * ( lgamma( alpha[i] ) - lgamma( proposed_alpha ) ) + ( proposed_alpha - alpha[i] ) * sum( log( pr[pr > 0] ) ) )
			
			if( accept_alpha > runif( 1 ) ) alpha[i] <- proposed_alpha 
		}
		
		if( i_mcmc > burnin )
		{
			pi_sample[counter, ]    <- pi
			alpha_sample[counter, ] <- alpha
			beta_sample[counter, ]  <- beta
			counter = counter + 1
		}
    }
        
	mcmc_sample = list( pi_sample = pi_sample, alpha_sample = alpha_sample, beta_sample = beta_sample, data = data, component_size = "fixed" )
	return( mcmc_sample )
}
   
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
# 7) k, alpha, beta, and pa: initial values for parameters respectively k, alpha, beta and pi
################################################################################
bmixgamma = function( data, k = "unknown", iter = 1000, burnin = iter / 2, lambda = 1, 
                    mu = NULL, nu = NULL, kesi = NULL, tau = NULL, 
                    k.start = NULL, alpha.start = NULL, beta.start = NULL, pi.start = NULL, 
                    k_max = 30, trace = TRUE )
{
	if( any( is.na( data ) ) ) stop( "Data should contain no missing data" ) 
	if( iter <= burnin )       stop( "Number of iteration must be more than number of burn-in" )	

	burnin = floor( burnin )
	n      = length( data )
		
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
	}
	else
	{
		component_size = "fixed"
	}

	if( is.null( pi.start ) )    pi.start    <- c( rep( 1 / k, k  ) )
	if( is.null( alpha.start ) ) alpha.start <- rgamma( k, mean( data ) ^ 2, var( data ) ) # moment estimation of alpha
	if( is.null( beta.start  ) ) beta.start  <- rgamma( k, mean( data )    , var( data ) ) # moment estimation of data	

	pi    = pi.start
	alpha = alpha.start
    beta  = beta.start
    
	if( component_size == "unknown" )
	{
		mcmc_sample = bmixgamma_unknown_k( data = data, n = n, k = k, iter = iter, burnin = burnin, lambda = lambda, 
									  mu = mu, nu = nu, kesi = kesi, tau = tau, 
									  alpha = alpha, beta = beta, pi = pi, 
									  k_max = k_max, trace = trace )		
	}
	else
	{
		mcmc_sample = bmixgamma_fixed_k( data = data, n = n, k = k, iter = iter, burnin = burnin, 
									  mu = mu, nu = nu, kesi = kesi, tau = tau, 
									  alpha = alpha, beta = beta, pi = pi, 
									  trace = trace )			
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
    
# summary of bestmix output
summary.bmixgamma = function( object, ... )
{
	data           = object $ data
	pi             = object $ pi_sample
	alpha          = object $ alpha_sample
	beta           = object $ beta_sample
	component_size = object $ component_size

	if( component_size == "unknown" )
	{
		all_k       = object $ all_k
		all_weights = object $ all_weights
		iter        = length( all_k )
		burnin      = iter - length( beta )
		k           = all_k[( burnin + 1 ):iter]
		weights     = all_weights[( burnin + 1 ):iter]
		
		op = par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) ) 
	}

	# plot for estimated distribution
	hist( data, prob = T, nclass = 25, col = "gray", border = "white"  )
	
	tt <- seq( min( data ) * 0.9, max( data ) * 1.2, length = 500 )
	f  <- 0 * tt

	if( component_size == "unknown" )
	{
		for ( i in 1:length(tt) )
			for ( j in 1:length( k ) ) 
				f[i] = f[i] + sum( pi[[j]] * dgamma( tt[i], alpha[[j]], beta[[j]] ) )
	
		lines( tt, f / length(k), col = "black", lty = 2, lw = 1 )
	}
	else
	{
		for ( i in 1:length( tt ) )
			for ( j in 1:nrow( pi ) ) 
				f[i] = f[i] + sum( pi[j,] * dgamma( tt[i], alpha[j,], beta[j,] ) )

		lines( tt, f / nrow(pi), col = "black", lty = 2, lw = 1 )
	}

    legend( "topright", c( "predictive density" ), lty = 2, col = "black", lwd = 1 )

	if( component_size == "unknown" )
	{
		# plot estimated distribution of k
		max_k = max( k )
		y     = vector( mode = "numeric", length = max_k )
		for ( i in 1:max_k ) y[i] <- sum( weights[k == i] ) 
		
		plot( x = 1:max_k, y, type = "h", main = "", ylab = "Pr(k|data)", xlab = "k(Number of components)" )

		# plot k based on iterations
		sum_weights    = 0 * weights
		sum_weights[1] = weights[1]
		for( i in 2:length( k ) ) sum_weights[i] = sum_weights[i - 1] + weights[i]

		plot( sum_weights, k, type = "l", xlab = "iteration", ylab = "Number of componants" )
	}
	else
	{
		cat( paste( "" ), fill = TRUE )
		cat( paste( "Number of mixture components = ", ncol( alpha ) ), fill = TRUE ) 
		cat( paste( "Estimated pi    = "), paste( round( apply( pi   , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated alpha = "), paste( round( apply( alpha, 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated beta  = "), paste( round( apply( beta , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "" ), fill = TRUE )				
	}
}  
     
# plot for class bestmix
plot.bmixgamma = function( x, ... )
{
	data           = x $ data
	pi             = x $ pi_sample
	alpha          = x $ alpha_sample
	beta           = x $ beta_sample
	component_size = x $ component_size
	
	# plot for estimated distribution
	hist( data, prob = T, nclass = 25, col = "gray", border = "white"  )
	
	tt <- seq( min( data ) * 0.9, max( data ) * 1.2, length = 500 )
	f  <- 0 * tt
	
	if( component_size == "unknown" )
	{
		all_k       = x $ all_k
		all_weights = x $ all_weights
		iter        = length( all_k )
		burnin      = iter - length( beta )
		k           = all_k[( burnin + 1 ):iter]
		weights     = all_weights[( burnin + 1 ):iter]

		for ( i in 1:length( tt ) )
			for ( j in 1:length( k ) ) 
				f[i] = f[i] + sum( pi[[j]] * dgamma( tt[i], alpha[[j]], beta[[j]] ) )
	
		lines( tt, f / length( k ), col = "black", lty = 2, lw = 1 )
	}
	else
	{
		for ( i in 1:length( tt ) )
			for ( j in 1:nrow( pi ) ) 
				f[i] = f[i] + sum( pi[j,] * dgamma( tt[i], alpha[j,], beta[j,] ) )

		lines( tt, f / nrow( pi ), col = "black", lty = 2, lw = 1 )
	}
	
    legend( "topright", c( "predictive density" ), lty = 2, col = "black", lwd = 1 )
}
     
# print of the bestmix output
print.bmixgamma = function( x, ... )
{
	if( ( x $ component_size ) == "unknown" )
	{
		all_k       = x $ all_k
		all_weights = x $ all_weights
		iter        = length( all_k )
		burnin      = iter - length( beta )
		k           = all_k[( burnin + 1 ):iter]
		weights     = all_weights[( burnin + 1 ):iter]

		max_k = max( k )
		y     = vector( mode = "numeric", length = max_k )
		for ( i in 1:max_k ) y[i] <- sum( weights[k == i] ) 

		cat( paste( "" ), fill = TRUE )
		cat( paste( "Estimated of number of components = ", which( y == max( y ) ) ), fill = TRUE ) 
		cat( paste( "" ), fill = TRUE )
	}
	else
	{
		cat( paste( "" ), fill = TRUE )
		cat( paste( "Number of mixture components = ",    ncol( x $ alpha_sample ) )              , fill = TRUE ) 
		cat( paste( "Estimated pi    = " ), paste( round( apply( x $ pi_sample   , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated alpha = " ), paste( round( apply( x $ alpha_sample, 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated beta  = " ), paste( round( apply( x $ beta_sample , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "" ), fill = TRUE )		
	}
} 
   
