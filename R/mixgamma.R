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
#    Random generation for the mixture of normal distribution 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

rmixgamma = function( n = 10, weight = 1, alpha = 1, beta = 1 ) 
{
	if( length( alpha ) != length( beta ) ) stop( " 'alpha' and 'beta' lengths differ" )

	component = sample( length( weight ), n, replace = TRUE, prob = weight ) 
	
	return( stats::rgamma( n, alpha[ component ], beta[ component ] ) ) 
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#   Density, distribution function, for the mixture of Gamma distribution 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

dmixgamma = function( x, weight = 1, alpha = 1, beta = 1 ) 
{
	if( length( alpha ) != length( beta ) ) stop( " 'alpha' and 'beta' lengths differ" )

	densmixgamma = 0
	for( i in 1:length( alpha ) )
		densmixgamma = densmixgamma + weight[ i ] * stats::dgamma( x, alpha[ i ], beta[ i ] )
	
	return( densmixgamma ) 
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
