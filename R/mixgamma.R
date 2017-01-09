# Random generation for the mixture of normal distribution 
rmixgamma = function( n = 10, weight = 1, alpha = 1, beta = 1 ) 
{
	if ( length( alpha ) != length( beta ) ) stop( "alpha and beta must be in the same size" )

	component = sample( length( weight ), n, replace = TRUE, prob = weight ) 
	
	return( rgamma( n, alpha[component], beta[component] ) ) 
}
    
# Density, distribution function, for the mixture of Gamma distribution 
dmixgamma = function( x, weight = 1, alpha = 1, beta = 1 ) 
{
	if ( length( alpha ) != length( beta ) ) stop( "alpha and beta must be in the same size" )

	densmixgamma = 0
	for( i in 1:length( alpha ) )
		densmixgamma = densmixgamma + weight[i] * dgamma( x, alpha[i], beta[i] )
	
	return( densmixgamma ) 
}
    
