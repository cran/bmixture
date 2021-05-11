## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2017 - 2021  Reza Mohammadi                                |
#                                                                              |
#     This file is part of ssgraph package.                                    |
#                                                                              |
#     "bmixture" is free software: you can redistribute it and/or modify it    |
#     under the terms of the GNU General Public License as published by the    |
#     Free Software Foundation: <https://cran.r-project.org/web/licenses/GPL-3>|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     MCMC sampling algorithm based on Birth-Death MCMC scheme
#     for mixture of Normal distribution with an unknown number of components  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bmixnorm_unknown_k = function( data, n, k, iter, burnin, lambda, 
                               epsilon, kappa, alpha, g, h,
                               mu, sig, pi, 
                               k.max, trace )
{
    all_k       = vector( mode = "numeric", length = iter )
    all_weights = all_k
    pi_sample   = list()  
    mu_sample   = pi_sample
    sig_sample  = pi_sample
    
    # - - - BDMCMC agorithm  - - - - - - - - - - - - - - - - - - - - - - - - - |
    counter = 1
    for( i_mcmc in 1:iter )
    {
        if( trace == TRUE && i_mcmc %% 100 == 0 )
        {
            mes <- paste( c( " iteration: ", i_mcmc, " from ", iter, ". Number of componets = ", k ), collapse = "" )
            cat( mes, "\r" )
            utils::flush.console()	
        }    
        
        # STEP 1: computing birth-death rates - - - - - - - - - - - - - - - - -|
        # first we obtain death_rates (vector of death rates for each component)
        if( k > 1 ) 
        { 
            death_rates = vector( mode = "numeric", length = k )
            
            for( j in 1:k )
            {
                log_death_rates <- 0
                for( t in 1:n )
                {
                    likelhood = sum( pi * stats::dnorm( data[ t ], mu, sqrt( sig ) ) )
                    if( likelhood == 0 ) likelhood <- .Machine$double.xmin # 10 ^ ( -320 )
                    log_death_rates = log_death_rates + log( likelhood - pi[ j ] * stats::dnorm( data[ t ], mu[ j ], sqrt( sig[ j ] ) ) ) - log( likelhood ) - log( 1 - pi[ j ] )
                }
                death_rates[ j ] <- exp( log_death_rates )
            }
            
            # checking for infinite value
            death_rates[ is.infinite( death_rates ) ] <- .Machine$double.xmax / 1000  # gamma( 171 )
        }else{
            death_rates = 0
        }
        
        # simulate the time s to the next jump, from an exponential distribution
        weight = 1 / ( lambda + sum( death_rates ) )
        # first Birth and death run and modify the parameters
        birthp <- lambda * weight
        
        if( ( stats::runif( 1 ) < birthp ) && ( k < k.max ) )
        {
            beta   <- stats::rgamma( 1, g, h )
            pi_new <- stats::rbeta( 1, 1, k )
            k      <- k + 1
            pi     <- c( pi * ( 1 - pi_new ), pi_new )
            mu     <- c( mu,  stats::rnorm( 1, epsilon, sqrt( 1 / kappa ) ) )
            sig    <- c( sig, 1 / stats::rgamma( 1, alpha, beta ) )
        }else{
            j     <- sample( 1:k, 1, prob = death_rates )
            k     <- k - 1
            pi    <- pi[  -j ] / ( 1 - pi[ j ] )
            mu    <- mu[  -j ]
            sig   <- sig[ -j ]
        }
        
        # STEP 2: Sample group memberships from multinominal distribution - - -|
        z  = matrix( 0, nrow = k, ncol = n ) 
        ## from multinomial distribution
        for( i in 1:n ) 
            z[ , i ] = stats::rmultinom( n = 1, size = 1, prob = pi * stats::dnorm( data[ i ], mean = mu, sd = sqrt( sig ) ) )
        
        # STEP 3: Sample parameters  - - - - - - - - - - - - - - - - - - - - - |
        # updating beta
        beta = stats::rgamma( 1, g + k * alpha, h + sum( 1 / sig ) )
        
        # updating pi
        n_i = apply( z, 1, sum )
        pi  = rdirichlet( 1, 1 + n_i )
        
        # updating mu
        for( i in 1:k ) 
        {
            ## Set mixture sample mean
            #~ 			sum_data = sum( data[z[i,] == 1] )
            sum_data = sum( z[ i, ] * data )
            ## Set posterior mean and standard deviation
            mu_sig  = 1 / ( n_i[ i ] / sig[ i ] + kappa ) 
            mu_mean = ( sum_data / sig[ i ] + kappa * epsilon ) * mu_sig
            
            ## Sample new mean from normal distribution
            mu[ i ] = stats::rnorm( n = 1, mean = mu_mean, sd = sqrt( mu_sig ) )
        }
        
        # updating sigma
        for( i in 1:k ) 
        {
            ## Set inverse gamma posterior parameters
            alpha_sig = alpha + n_i[ i ] / 2
            beta_sig  = beta  + ( sum( z[ i, ] * ( data - mu[ i ] ) ^ 2 ) ) / 2 
            ## Sample variance from gamma random variable and invert it
            sig[ i ] = 1 / stats::rgamma( n = 1, alpha_sig, beta_sig )
        }
        
        # STEP 4: Sort parameters based on mu - - - - - - - - - - - - - - - - -|
        order_mu = order( mu )
        pi       = pi[  order_mu ]
        mu       = mu[  order_mu ]
        sig      = sig[ order_mu ]
        z        = z[   order_mu, ]
        
        # Saving MCMC samples - - - - - - - - - - - - - - - - - - - - - - - - -|
        all_k[ i_mcmc ]       = k # saving all values of k for chicking convergency
        all_weights[ i_mcmc ] = weight
        
        if( i_mcmc > burnin )
        {
            pi_sample[[  counter ]] <- pi
            mu_sample[[  counter ]] <- mu
            sig_sample[[ counter ]] <- sig
            counter = counter + 1
        }
    } # End of BDMCMC sampling algorithm - - - - - - - - - - - - - - - - - - - |
    
    mcmc_sample = list( all_k = all_k, all_weights = all_weights, pi_sample = pi_sample, mu_sample = mu_sample, sig_sample = sig_sample, data = data, component_size = "unknown" )    
    
    return( mcmc_sample )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  MCMC sampling algorithm 
#  for mixture of Normal distribution with an fixed number of components  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

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
        if( trace == TRUE && i_mcmc %% 10 == 0 )
        {
            mes <- paste( c( " iteration: ", i_mcmc, " from ", iter, ". Number of componets = ", k ), collapse = "" )
            cat( mes, "\r" )
            utils::flush.console()	
        }    
        
        # Sample group memberships from multinominal distribution
        for( i in 1:n ) 
            z[ , i ] = stats::rmultinom( n = 1, size = 1, prob = pi * stats::dnorm( data[ i ], mean = mu, sd = sqrt( sig ) ) )
        
        ## --- Sample parameters
        
        # step 2: updating beta
        beta = stats::rgamma( 1, g + k * alpha, h + sum( 1 / sig ) )
        
        # step 3: updating pi
        n_i = apply( z, 1, sum )
        pi  = rdirichlet( 1, 1 + n_i )
        
        # spep 4: updating mu
        for( i in 1:k ) 
        {
            ## Set mixture sample mean
            #~ 			sum_data = sum( data[z[i,] == 1] )
            sum_data = sum( z[ i, ] * data )
            
            ## Set posterior mean and standard deviation
            mu_sig  = 1 / ( n_i[ i ] / sig[ i ] + kappa ) 
            mu_mean = ( sum_data / sig[ i ] + kappa * epsilon ) * mu_sig
            
            ## Sample new mean from normal distribution
            mu[ i ] = stats::rnorm( n = 1, mean = mu_mean, sd = sqrt( mu_sig ) )
        }
        
        # spep 5: updating sigma
        for( i in 1:k ) 
        {
            ## Set inverse gamma posterior parameters
            alpha_sig = alpha + n_i[ i ] / 2
            beta_sig  = beta  + ( sum( z[ i, ] * ( data - mu[ i ] ) ^ 2 ) ) / 2 
            ## Sample variance from gamma random variable and invert it
            sig[ i ]  = 1 / stats::rgamma( n = 1, alpha_sig, beta_sig )
        }
        
        # Sort parameters based on mu
        order_mu = order( mu )
        pi       = pi[  order_mu ]
        mu       = mu[  order_mu ]
        sig      = sig[ order_mu ]
        z        = z[   order_mu, ]
        
        if( i_mcmc > burnin )
        {
            pi_sample[  counter, ] = pi
            mu_sample[  counter, ] = mu
            sig_sample[ counter, ] = sig
            counter = counter + 1
        }
    }
    
    mcmc_sample = list( pi_sample = pi_sample, mu_sample = mu_sample, sig_sample = sig_sample, data = data, component_size = "fixed" )    
    return( mcmc_sample )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
