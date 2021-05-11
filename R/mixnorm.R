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
#   Random generation for the mixture of normal distribution 
#   with mean equal to mean and standard deviation equal to sd.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

rmixnorm = function( n = 10, weight = 1, mean = 0, sd = 1 ) 
{
	if( length( mean ) != length( sd ) ) stop( " 'mean' and 'sd' lengths differ" )

	component = sample( length( weight ), n, replace = TRUE, prob = weight ) 
	
	return( stats::rnorm( n, mean = mean[ component ], sd = sd[ component ] ) ) 
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#   Density, distribution function, for the mixture of normal distribution 
#   with mean equal to mean and standard deviation equal to sd.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

dmixnorm = function( x, weight = 1, mean = 0, sd = 1 ) 
{
	if( length( mean ) != length( sd ) ) stop( " 'mean' and 'sd' lengths differ" )
	
	densmixnorm = 0
	for( i in 1:length( mean ) )
		densmixnorm = densmixnorm + weight[ i ] * stats::dnorm( x, mean[ i ], sd[ i ] )
	
	return( densmixnorm ) 
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
