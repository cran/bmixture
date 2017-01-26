# Random generator from Dirichlet distribution
rdirichlet = function( n = 10, alpha = c( 1, 1 ) ) 
{
    length_alpha = length( alpha )
    sample       = matrix( rgamma( length_alpha * n, alpha ), ncol = length_alpha, byrow = TRUE )
    
    return( sample / apply( sample, 1, sum ) )
}
   
