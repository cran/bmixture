
# MCMC sampling algorithm based on Birth-Death MCMC scheme
# for mixture of Normal distribution with an unknown number of components  
bmixnorm_unknown_k = function( data, n, k, iter, burnin, lambda, 
									  epsilon, kappa, alpha, g, h,
									  mu, sig, pi, 
									  k_max, trace )
{
	all_k       = vector( mode = "numeric", length = iter )
	all_weights = all_k
	pi_sample   = list()  
    mu_sample   = pi_sample
    sig_sample  = pi_sample
 
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
			
			for( j in 1:k )
			{
				log_death_rates <- 0
				for ( t in 1:n )
				{
					likelhood = sum( pi * dnorm( data[t], mu, sqrt( sig ) ) )
					if( likelhood == 0 ) likelhood = 10 ^ ( -320 )
					log_death_rates = log_death_rates + log( likelhood - pi[j] * dnorm( data[t], mu[j], sqrt( sig[j] ) ) ) - log( likelhood ) - log( 1 - pi[j] )
				}
				death_rates[j] <- exp( log_death_rates )
			}
			
			# checking for infinite value
			death_rates[is.infinite( death_rates )] <- gamma( 171 )
		}else{
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
			mu     <- c( mu,  rnorm( 1, epsilon, sqrt( 1 / kappa ) ) )
			sig    <- c( sig, 1 / rgamma( 1, alpha, beta ) )
		}else{
			j     <- sample( 1:k, 1, prob = death_rates )
			k     <- k - 1
			pi    <- pi[-j] / ( 1 - pi[j] )
			mu    <- mu[-j]
			sig   <- sig[-j]
		}

		# Sample group memberships from multinominal distribution
		z  = matrix( 0, nrow = k, ncol = n ) 
		## from multinomial distribution
		for( i in 1:n ) 
			z[, i] = rmultinom( n = 1, size = 1, prob = pi * dnorm( data[i], mean = mu, sd = sqrt( sig ) ) )
		
		## --- Sample parameters
		
		# step 2: updating beta
		beta = rgamma( 1, g + k * alpha, h + sum( 1 / sig ) )
		
		# step 3: updating pi
		n_i = apply( z, 1, sum )
		pi  = rdirichlet( 1, 1 + n_i )
				
		# spep 4: updating mu
		for( i in 1:k ) 
		{
			## Set mixture sample mean
#~ 			sum_data = sum( data[z[i,] == 1] )
			sum_data = sum( z[i,] * data )
			## Set posterior mean and standard deviation
			mu_sig  = 1 / ( n_i[i] / sig[i] + kappa ) 
			mu_mean = ( sum_data / sig[i] + kappa * epsilon ) * mu_sig

			## Sample new mean from normal distribution
			mu[i] = rnorm( n = 1, mean = mu_mean, sd = sqrt( mu_sig ) )
		}
		
		# spep 5: updating sigma
		for( i in 1:k ) 
		{
			## Set inverse gamma posterior parameters
			alpha_sig = alpha + n_i[i] / 2
			beta_sig  = beta  + ( sum( z[i,] * ( data - mu[i] ) ^ 2 ) ) / 2 
			## Sample variance from gamma random variable and invert it
			sig[i] = 1 / rgamma( n = 1, alpha_sig, beta_sig )
		}

		# Sort parameters based on mu
		order_mu = order( mu )
		pi       = pi[order_mu]
		mu       = mu[order_mu]
		sig      = sig[order_mu]
		z        = z[order_mu, ]
	
		all_k[i_mcmc]       = k # saving all values of k for chicking convergency
		all_weights[i_mcmc] = weight
		
		if ( i_mcmc > burnin )
		{
			pi_sample[[    counter ]] <- pi
			mu_sample[[ counter ]]    <- mu
			sig_sample[[  counter ]]  <- sig
			counter = counter + 1
		}
    }

	mcmc_sample = list( all_k = all_k, all_weights = all_weights, pi_sample = pi_sample, mu_sample = mu_sample, sig_sample = sig_sample, data = data, component_size = "unknown" )    
	return( mcmc_sample )
}
        
# MCMC sampling algorithm 
# for mixture of Normal distribution with an fixed number of components  
bmixnorm_fixed_k = function( data, n, k, iter, burnin, 
									  epsilon, kappa, alpha, g, h,
									  mu, sig, pi, 
									  trace )
{
	pi_sample  = matrix( NA, nrow = iter - burnin, ncol = k ) 
    mu_sample  = pi_sample
    sig_sample = pi_sample
  
	## Intialise matrix of mixing proportion times normal likelihood
	z  = matrix( 0, nrow = k, ncol = n ) 
  
	counter = 1
    for( i_mcmc in 1:iter )
    {
		if ( trace == TRUE && i_mcmc %% 10 == 0 )
		{
			mes <- paste( c( " iteration: ", i_mcmc, " from ", iter, ". Number of componets = ", k ), collapse = "" )
			cat( mes, "\r" )
			flush.console()	
		}    

		# Sample group memberships from multinominal distribution
		for( i in 1:n ) 
			z[, i] = rmultinom( n = 1, size = 1, prob = pi * dnorm( data[i], mean = mu, sd = sqrt( sig ) ) )
		
		## --- Sample parameters
		
		# step 2: updating beta
		beta = rgamma( 1, g + k * alpha, h + sum( 1 / sig ) )
		
		# step 3: updating pi
		n_i = apply( z, 1, sum )
		pi  = rdirichlet( 1, 1 + n_i )
				
		# spep 4: updating mu
		for( i in 1:k ) 
		{
			## Set mixture sample mean
#~ 			sum_data = sum( data[z[i,] == 1] )
			sum_data = sum( z[i,] * data )
			
			## Set posterior mean and standard deviation
			mu_sig  = 1 / ( n_i[i] / sig[i] + kappa ) 
			mu_mean = ( sum_data / sig[i] + kappa * epsilon ) * mu_sig

			## Sample new mean from normal distribution
			mu[i] = rnorm( n = 1, mean = mu_mean, sd = sqrt( mu_sig ) )
		}
		
		# spep 5: updating sigma
		for( i in 1:k ) 
		{
			## Set inverse gamma posterior parameters
			alpha_sig = alpha + n_i[i] / 2
			beta_sig  = beta  + ( sum( z[i,] * ( data - mu[i] ) ^ 2 ) ) / 2 
			## Sample variance from gamma random variable and invert it
			sig[i] = 1 / rgamma( n = 1, alpha_sig, beta_sig )
		}

		# Sort parameters based on mu
		order_mu = order( mu )
		pi       = pi[order_mu]
		mu       = mu[order_mu]
		sig      = sig[order_mu]
		z        = z[order_mu, ]
	
		if( i_mcmc > burnin )
		{
			pi_sample[counter, ]  = pi
			mu_sample[counter, ]  = mu
			sig_sample[counter, ] = sig
			counter = counter + 1
		}
    }

	mcmc_sample = list( pi_sample = pi_sample, mu_sample = mu_sample, sig_sample = sig_sample, data = data, component_size = "fixed" )    
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
# 7) k, mu, sig, and pa: initial values for parameters respectively k, mu, sig and pi
################################################################################
bmixnorm = function( data, k = "unknown", iter = 1000, burnin = iter / 2, lambda = 1, 
                    k.start = NULL, mu.start = NULL, sig.start = NULL, pi.start = NULL, 
                    k_max = 30, trace = TRUE )
{
	if( any( is.na( data ) ) ) stop( "Data should contain no missing data" ) 
	if( iter <= burnin )       stop( "Number of iteration must be more than number of burn-in" )	

	burnin = floor( burnin )
	n      = length( data )
	
	max_data = max( data )	
	min_data = min( data )
	R        = max_data - min(data)	
	
	# Values for paprameters of prior distributon of mu
	epsilon = R / 2   # midpoint of the observed range of the data
	kappa   = 1 / R ^ 2
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

	beta = rgamma( 1, g, h )
	if( is.null( pi.start  ) ) pi.start  = c( rep( 1 / k, k  ) )
	if( is.null( mu.start  ) ) mu.start  = rnorm( k, epsilon, sqrt( 1 / kappa ) ) 
	if( is.null( sig.start ) ) sig.start = 1 / rgamma( k, alpha, beta ) 	

	pi  = pi.start
	mu  = mu.start
    sig = sig.start

	# Sort parameters based on mu
	order_mu = order( mu )
	pi       = pi[order_mu]
	mu       = mu[order_mu]
	sig      = sig[order_mu]
    
############### MCMC 
	if( component_size == "unknown" )
	{
		mcmc_sample = bmixnorm_unknown_k( data = data, n = n, k = k, iter = iter, burnin = burnin, lambda = lambda, 
									  epsilon = epsilon, kappa = kappa, alpha = alpha, g = g, h = h,
									  mu = mu, sig = sig, pi = pi, 
									  k_max = k_max, trace = trace )		
	}
	else
	{
		mcmc_sample = bmixnorm_fixed_k( data = data, n = n, k = k, iter = iter, burnin = burnin, 
									  epsilon = epsilon, kappa = kappa, alpha = alpha, g = g, h = h,
									  mu = mu, sig = sig, pi = pi, 
									  trace = trace )			
	}

	if( trace == TRUE )
	{
		mes <- paste( c(" ", iter," iteration done.                               " ), collapse = "" )
		cat( mes, "\r" )
		cat( "\n" )
		flush.console()
	}    

	class( mcmc_sample ) = "bmixnorm"
	return( mcmc_sample )
}
   
# summary of bestmix output
summary.bmixnorm = function( object, ... )
{
	data           = object $ data
	pi             = object $ pi_sample
	mu             = object $ mu_sample
	sig            = object $ sig_sample
	component_size = object $ component_size

	if( component_size == "unknown" )
	{
		all_k       = object $ all_k
		all_weights = object $ all_weights
		iter        = length( all_k )
		burnin      = iter - length( object $ pi_sample )
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
			for ( j in 1:length(k) ) 
				f[i] = f[i] + sum( pi[[j]] * dnorm( tt[i], mu[[j]], sqrt( sig[[j]] ) ) )
	
		lines( tt, f / length(k), col = "black", lty = 2, lw = 1 )
	}
	else
	{
		for ( i in 1:length(tt) )
			for ( j in 1:nrow( pi ) ) 
				f[i] = f[i] + sum( pi[j,] * dnorm( tt[i], mu[j,], sqrt( sig[j,] ) ) )

		lines( tt, f / nrow( pi ), col = "black", lty = 2, lw = 1 )
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
		sum_weights = 0 * weights
		sum_weights[1] = weights[1]
		for( i in 2:length( k ) ) sum_weights[i] = sum_weights[i - 1] + weights[i]

		plot( sum_weights, k, type = "l", xlab = "iteration", ylab = "Number of componants" )
	}
	else
	{
		cat( paste( "" ), fill = TRUE )
		cat( paste( "Number of mixture components = ", ncol( mu ) ), fill = TRUE ) 
		cat( paste( "Estimated pi  = "), paste( round( apply( pi   , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated mu  = "), paste( round( apply( mu, 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated sig = "), paste( round( apply( sig , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "" ), fill = TRUE )				
	}
}  
     
# plot for class bestmix
plot.bmixnorm = function( x, ... )
{
	data           = x $ data
	pi             = x $ pi_sample
	mu          = x $ mu_sample
	sig           = x $ sig_sample
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
		sample_size = length( sig )
		burnin      = iter - sample_size
		k           = all_k[( burnin + 1 ):iter]
		weights     = all_weights[( burnin + 1 ):iter]

		for ( i in 1:length(tt) )
			for ( j in 1:sample_size ) 
				f[i] = f[i] + sum( pi[[j]] * dnorm( tt[i], mu[[j]], sqrt( sig[[j]] ) ) )
	
		lines( tt, f / sample_size, col = "black", lty = 2, lw = 1 )
	}
	else
	{
		for ( i in 1:length(tt) )
			for ( j in 1:nrow(pi) ) 
				f[i] = f[i] + sum( pi[j,] * dnorm( tt[i], mu[j,], sqrt( sig[j,] ) ) )

		lines( tt, f / nrow( pi ), col = "black", lty = 2, lw = 1 )
	}
	
    legend( "topright", c( "predictive density" ), lty = 2, col = "black", lwd = 1 )
}
     
# print of the bestmix output
print.bmixnorm = function( x, ... )
{
	if( x $ component_size == "unknown" )
	{
		all_k       = x $ all_k
		all_weights = x $ all_weights
		iter        = length( all_k )
		burnin      = iter - length( x $ pi_sample )
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
		cat( paste( "Number of mixture components = ", ncol( x $ mu_sample ) ), fill = TRUE ) 
		cat( paste( "Estimated pi  = " ), paste( round( apply( x $ pi_sample  , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated mu  = " ), paste( round( apply( x $ mu_sample  , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "Estimated sig = " ), paste( round( apply( x $ sig_sample , 2, mean ), 2 ) ), fill = TRUE ) 
		cat( paste( "" ), fill = TRUE )		
	}
} 
   
