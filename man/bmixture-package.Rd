\name{bmixture-package}
\alias{bmixture-package}

\alias{ bmixnorm_unknown_k }
\alias{ bmixnorm_fixed_k }

\docType{package}

\title{ Bayesian Estimation for Finite Mixture of Distributions }

\description{
The \code{R} package \pkg{bmixture} provides statistical tools for Bayesian estimation in finite mixture of distributions.       
The package implemented the improvements in the Bayesian literature, including Mohammadi and Salehi-Rad (2012) and Mohammadi et al. (2013).
Besides, the package contains several functions for simulation and visualization, as well as a real dataset taken from the literature.
}

\section{How to cite this package}{
Whenever using this package, please cite as
%\preformatted{
Mohammadi R. (2019). \pkg{bmixture}: Bayesian Estimation for Finite Mixture of 
Distributions, \code{R} package version 1.3, https://CRAN.R-project.org/package=bmixture
%}
}

\author{Reza Mohammadi <a.mohammadi@uva.nl>}

\references{
Mohammadi, A., Salehi-Rad, M. R., and Wit, E. C. (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics}, 28(2):683-700

Mohammadi, A, and Salehi-Rad, M. R. (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation}, 41(3):419-435

Stephens, M. (2000) Bayesian analysis of mixture models with an unknown number of components-an alternative to reversible jump methods. \emph{Annals of statistics}, 28(1):40-74

Richardson, S. and Green, P. J. (1997) On Bayesian analysis of mixtures with an unknown number of components. \emph{Journal of the Royal Statistical Society: series B}, 59(4):731-792

Green, P. J. (1995) Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. \emph{Biometrika}, 82(4):711-732

Cappe, O., Christian P. R., and Tobias, R. (2003) Reversible jump, birth and death and more general continuous time Markov chain Monte Carlo samplers. \emph{Journal of the Royal Statistical Society: Series B}, 65(3):679-700

Wade, S. and Ghahramani, Z. (2018) Bayesian Cluster Analysis: Point Estimation and Credible Balls (with Discussion). \emph{Bayesian Analysis}, 13(2):559-626
}

\examples{
\dontrun{

library( bmixture )

data( galaxy )

# Runing bdmcmc algorithm for the galaxy dataset      
mcmc_sample = bmixnorm( data = galaxy )

summary( mcmc_sample ) 
plot( mcmc_sample )
print( mcmc_sample)

# simulating data from mixture of Normal with 3 components
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
     
# Runing bdmcmc algorithm for the above simulation data set      
bmixnorm.obj = bmixnorm( data, k = 3, iter = 1000 )
    
summary( bmixnorm.obj ) 
}
}

\keyword{package}
