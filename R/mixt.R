## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2017 - 2019  Reza Mohammadi                                                    |
#                                                                                                  |
#     This file is part of ssgraph package.                                                        |
#                                                                                                  |
#     "bmixture" is free software: you can redistribute it and/or modify it under                  |
#     the terms of the GNU General Public License as published by the Free                         |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                  |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#   Random generation for the mixture of t-distribution 
#   with mean equal to mean and standard deviation equal to sd.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

rmixt = function( n = 10, weight = 1, df = 1, mean = 0, sd = 1 ) 
{
	length_mean = length( mean )
	if( length_mean != length( sd ) ) stop( "mean and sd must be in the same size" )
	if( length( df ) == 1 ) df = rep( df, length_mean )
   
	component = sample( length( weight ), n, replace = TRUE, prob = weight ) 
	
	return( mean[ component ] + sd[ component ] * stats::rt( n, df = df[ component ] ) ) 
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#   Density, distribution function, for the mixture of t-distribution 
#   with mean equal to mean and standard deviation equal to sd.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
dmixt = function( x, weight = 1, df = 1, mean = 0, sd = 1 ) 
{
	length_mean = length( mean )
	if( length_mean != length( sd ) ) stop( "mean and sd must be in the same size" )
	if( length( df ) == 1 ) df = rep( df, length_mean )
	
	densmixt = 0
	for( i in 1:length_mean )
		densmixt = densmixt + weight[ i ] * stats::dt( ( x - mean[ i ] ) / sd[ i ], df[ i ] )
	
	return( densmixt ) 
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
