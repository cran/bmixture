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
#     Main function: BDMCMC algorithm for finite mixture of Gamma distribution
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# INPUT for bdmcmc funciton 
# 1) data:         the data with posetive and no missing values
# 2) k             number of components of mixture distribution. Defult is unknown
# 2) iter:         nuber of iteration of the algorithm 
# 3) burnin:       number of burn-in iteration
# 4) lambda_r:     rate for birth and parameter of prior distribution of k
# 7) k, mu, sig, and pa: initial values for parameters respectively k, mu, sig and pi
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bmixnorm = function( data, k = "unknown", iter = 1000, burnin = iter / 2, lambda = 1, 
                     k.start = NULL, mu.start = NULL, sig.start = NULL, pi.start = NULL, 
                     k.max = 30, trace = TRUE )
{
    if( any( is.na( data ) ) ) stop( "'data' must contain no missing values" ) 
    if( iter <= burnin )       stop( "'iter' must be higher than 'burnin'" )	
    
    burnin = floor( burnin )
    n      = length( data )
    
    max_data = max( data )	
    min_data = min( data )
    R        = max_data - min_data	
    
    # Values for paprameters of prior distributon of mu
    epsilon = R / 2        # midpoint of the observed range of the data
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
        if( is.null( k.start ) ){ k = 3 }else{ k = k.start }
    }else{
        component_size = "fixed"
    }
    
    beta = stats::rgamma( 1, g, h )
    if( is.null( pi.start  ) ) pi.start  = c( rep( 1 / k, k  ) )
    if( is.null( mu.start  ) ) mu.start  = stats::rnorm( k, epsilon, sqrt( 1 / kappa ) ) 
    if( is.null( sig.start ) ) sig.start = 1 / stats::rgamma( k, alpha, beta ) 	
    
    pi  = pi.start
    mu  = mu.start
    sig = sig.start
    
    # Sort parameters based on mu
    order_mu = order( mu )
    pi       = pi[  order_mu ]
    mu       = mu[  order_mu ]
    sig      = sig[ order_mu ]
    
    ############### MCMC 
    if( component_size == "unknown" )
    {
        mcmc_sample = bmixnorm_unknown_k( data = data, n = n, k = k, iter = iter, burnin = burnin, lambda = lambda, 
                                          epsilon = epsilon, kappa = kappa, alpha = alpha, g = g, h = h,
                                          mu = mu, sig = sig, pi = pi, 
                                          k.max = k.max, trace = trace )		
    }else{
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
        utils::flush.console()
    }    
    
    class( mcmc_sample ) = "bmixnorm"
    return( mcmc_sample )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#   summary of bestmix output
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
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
        k           = all_k[ ( burnin + 1 ):iter ]
        weights     = all_weights[ ( burnin + 1 ):iter ]
        
        op = graphics::par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) ) 
    }
    
    # plot for estimated distribution
    graphics::hist( data, prob = T, nclass = 25, col = "gray", border = "white"  )
    
    tt <- seq( min( data ) * 0.9, max( data ) * 1.2, length = 500 )
    f  <- 0 * tt
    
    if( component_size == "unknown" )
    {
        for( i in 1:length( tt ) )
            for( j in 1:length( k ) ) 
                f[ i ] = f[ i ] + sum( pi[[ j ]] * stats::dnorm( tt[ i ], mu[[ j ]], sqrt( sig[[ j ]] ) ) )
            
            graphics::lines( tt, f / length(k), col = "black", lty = 2, lw = 1 )
    }else{
        for( i in 1:length( tt ) )
            for( j in 1:nrow( pi ) ) 
                f[ i ] = f[ i ] + sum( pi[ j, ] * stats::dnorm( tt[ i ], mu[ j, ], sqrt( sig[ j, ] ) ) )
            
            graphics::lines( tt, f / nrow( pi ), col = "black", lty = 2, lw = 1 )
    }
    
    graphics::legend( "topright", c( "predictive density" ), lty = 2, col = "black", lwd = 1 )
    
    if( component_size == "unknown" )
    {
        # plot estimated distribution of k
        max_k = max( k )
        y     = vector( mode = "numeric", length = max_k )
        for( i in 1:max_k ) y[ i ] <- sum( weights[ k == i ] ) 
        
        graphics::plot( x = 1:max_k, y, type = "h", main = "", ylab = "Pr(k|data)", xlab = "k(Number of components)" )
        
        # plot k based on iterations
        sum_weights      = 0 * weights
        sum_weights[ 1 ] = weights[ 1 ]
        for( i in 2:length( k ) ) sum_weights[ i ] = sum_weights[ i - 1 ] + weights[ i ]
        
        graphics::plot( sum_weights, k, type = "l", xlab = "iteration", ylab = "Number of componants" )
    }else{
        cat( paste( "" ), fill = TRUE )
        cat( paste( "Number of mixture components = ", ncol( mu ) ), fill = TRUE ) 
        cat( paste( "Estimated 'pi'  = "), paste( round( apply( pi , 2, mean ), 2 ) ), fill = TRUE ) 
        cat( paste( "Estimated 'mu'  = "), paste( round( apply( mu , 2, mean ), 2 ) ), fill = TRUE ) 
        cat( paste( "Estimated 'sig' = "), paste( round( apply( sig, 2, mean ), 2 ) ), fill = TRUE ) 
        cat( paste( "" ), fill = TRUE )				
    }
}  

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#   plot for class bestmix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

plot.bmixnorm = function( x, ... )
{
    data           = x $ data
    pi             = x $ pi_sample
    mu             = x $ mu_sample
    sig            = x $ sig_sample
    component_size = x $ component_size
    
    # plot for estimated distribution
    graphics::hist( data, prob = T, nclass = 25, col = "gray", border = "white"  )
    
    tt <- seq( min( data ) * 0.9, max( data ) * 1.2, length = 500 )
    f  <- 0 * tt
    
    if( component_size == "unknown" )
    {
        all_k       = x $ all_k
        all_weights = x $ all_weights
        iter        = length( all_k )
        sample_size = length( sig )
        burnin      = iter - sample_size
        k           = all_k[ ( burnin + 1 ):iter ]
        weights     = all_weights[ ( burnin + 1 ):iter ]
        
        for( i in 1:length( tt ) )
            for( j in 1:sample_size ) 
                f[ i ] = f[ i ] + sum( pi[[ j ]] * stats::dnorm( tt[ i ], mu[[ j ]], sqrt( sig[[ j ]] ) ) )
        
        graphics::lines( tt, f / sample_size, col = "black", lty = 2, lw = 1 )
    }else{
        for( i in 1:length( tt ) )
            for( j in 1:nrow( pi ) ) 
                f[ i ] = f[ i ] + sum( pi[ j, ] * stats::dnorm( tt[ i ], mu[ j, ], sqrt( sig[ j, ] ) ) )
            
            graphics::lines( tt, f / nrow( pi ), col = "black", lty = 2, lw = 1 )
    }
    
    graphics::legend( "topright", c( "predictive density" ), lty = 2, col = "black", lwd = 1 )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#   print of the bestmix output
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

print.bmixnorm = function( x, ... )
{
    if( x $ component_size == "unknown" )
    {
        all_k       = x $ all_k
        all_weights = x $ all_weights
        iter        = length( all_k )
        burnin      = iter - length( x $ pi_sample )
        k           = all_k[ ( burnin + 1 ):iter ]
        weights     = all_weights[ ( burnin + 1 ):iter ]
        
        max_k = max( k )
        y     = vector( mode = "numeric", length = max_k )
        for( i in 1:max_k ) y[ i ] <- sum( weights[ k == i ] ) 
        
        cat( paste( "" ), fill = TRUE )
        cat( paste( "Estimation for the number of components = ", which( y == max( y ) ) ), fill = TRUE ) 
        cat( paste( "" ), fill = TRUE )
    }else{
        cat( paste( "" ), fill = TRUE )
        cat( paste( "Number of mixture components = ", ncol( x $ mu_sample ) ), fill = TRUE ) 
        cat( paste( "Estimated 'pi'  = " ), paste( round( apply( x $ pi_sample  , 2, mean ), 2 ) ), fill = TRUE ) 
        cat( paste( "Estimated 'mu'  = " ), paste( round( apply( x $ mu_sample  , 2, mean ), 2 ) ), fill = TRUE ) 
        cat( paste( "Estimated 'sig' = " ), paste( round( apply( x $ sig_sample , 2, mean ), 2 ) ), fill = TRUE ) 
        cat( paste( "" ), fill = TRUE )		
    }
} 

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
