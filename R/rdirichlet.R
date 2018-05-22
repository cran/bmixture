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
## Random generator from Dirichlet distribution
## ------------------------------------------------------------------------------------------------|
rdirichlet = function( n = 10, alpha = c( 1, 1 ) ) 
{
    length_alpha = length( alpha )
    sample       = matrix( rgamma( length_alpha * n, alpha ), ncol = length_alpha, byrow = TRUE )
    
    return( sample / apply( sample, 1, sum ) )
}
   
