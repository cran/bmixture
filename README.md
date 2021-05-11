# **bmixture** 
  
![](https://www.r-pkg.org/badges/version/bmixture) 
![](https://www.r-pkg.org/badges/last-release/bmixture) 
![](https://cranlogs.r-pkg.org/badges/bmixture) 
![](https://cranlogs.r-pkg.org/badges/grand-total/bmixture) 


The `R` package **bmixture** provides statistical tools for Bayesian estimation for the mixture of distributions. The package implemented the improvements in the Bayesian literature, including [Mohammadi et al. (2013)](https://link.springer.com/article/10.1007/s00180-012-0323-3) and [Mohammadi and Salehi-Rad (2012)](https://www.tandfonline.com/doi/full/10.1080/03610918.2011.588358).
Besides, the package contains several functions for simulation and visualization, as well as a real dataset taken from the literature.

## Installation

You can install the latest version from CRAN using:

``` r
install.packages( "bmixture" )
```

``` r
require( "bmixture" )
```

## Example 1: Finite mixture of Normal distributions using real world data

Here is a simple example to see the performance of the package for the Finite mixture of Normal distributions for the `galaxy` dataset:

``` r
data( galaxy )

# Runing bdmcmc algorithm for the galaxy dataset      
mcmc_sample = bmixnorm( data = galaxy )

summary( mcmc_sample ) 
plot( mcmc_sample )
print( mcmc_sample )
```

## Example 2: Finite mixture of Normal distributions using simulatoin data

Here is a simple example to see the performance of the package for the Finite mixture of Normal distributions using simulation data. First, we simulate data from the mixture of Normal with 3 components as follow:

``` r
n      = 500
mean   = c( 0  , 10 , 3   )
sd     = c( 1  , 1  , 1   )
weight = c( 0.3, 0.5, 0.2 )
    
data = rmixnorm( n = n, weight = weight, mean = mean, sd = sd )
   
# plot for simulation data      
hist( data, prob = TRUE, nclass = 30, col = "gray" )
    
x           = seq( -20, 20, 0.05 )
densmixnorm = dmixnorm( x, weight, mean, sd )
      
lines( x, densmixnorm, lwd = 2 )  
```

Now, we run the 'bdmcmc' algorithm for the above simulation data set 
     
``` r
bmixnorm.obj = bmixnorm( data, k = 3, iter = 1000 )
    
summary( bmixnorm.obj ) 
```



