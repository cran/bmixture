## ------------------------------------------------------------------------------------------------|
#     Copyright (C) 2017 - 2018  Reza Mohammadi                                                    |
#                                                                                                  |
#     This file is part of ssgraph package.                                                        |
#                                                                                                  |
#     "bmixture" is free software: you can redistribute it and/or modify it under                  |
#     the terms of the GNU General Public License as published by the Free                         |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                  |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                              |
## ------------------------------------------------------------------------------------------------|
## Random generation for the mixture of normal distribution 
# with mean equal to mean and standard deviation equal to sd.
## ------------------------------------------------------------------------------------------------|
rmixnorm = function( n = 10, weight = 1, mean = 0, sd = 1 ) 
{
	if ( length( mean ) != length( sd ) ) stop( "mean and sd must be in the same size" )

	component = sample( length( weight ), n, replace = TRUE, prob = weight ) 
	
	return( rnorm( n, mean = mean[component], sd = sd[component] ) ) 
}
    
# Density, distribution function, for the mixture of normal distribution 
# with mean equal to mean and standard deviation equal to sd.
dmixnorm = function( x, weight = 1, mean = 0, sd = 1 ) 
{
	if ( length( mean ) != length( sd ) ) stop( "mean and sd must be in the same size" )
	
	densmixnorm = 0
	for( i in 1:length( mean ) )
		densmixnorm = densmixnorm + weight[i] * dnorm( x, mean[i], sd[i] )
	
	return( densmixnorm ) 
}
    
