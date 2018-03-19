\name{print.bmixt}
\alias{print.bmixt}

\title{Print function for \code{S3} class \code{"bmixt"} }

\description{ Prints the information about the output of function \code{\link{bmixt}}. }

\usage{\method{print}{bmixt}( x, ... )}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"bmixt"}, from function \code{\link{bmixt}}. }
  \item{\dots}{System reserved (no specific usage). }
}

\references{
Mohammadi, A., M. R. Salehi-Rad, and E. C. Wit (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics} 28(2):683-700

Mohammadi, A., and M. R. Salehi-Rad (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation} 41(3):419-435
}

\author{ Reza Mohammadi }

\seealso{\code{\link{bmixt}}}

\examples{
\dontrun{
# simulating data from mixture of Normal with 3 components
n      = 500
weight = c( 0.3, 0.5, 0.2 )
mean   = c( 0  , 10 , 3   )
sd     = c( 1  , 1  , 1   )
    
data = rmixnorm( n = n, weight = weight, mean = mean, sd = sd )

# plot for simulation data      
hist( data, prob = TRUE, nclass = 30, col = "gray" )
  
x           = seq( -20, 20, 0.05 )
densmixnorm = dmixnorm( x, weight, mean, sd )
      
lines( x, densmixnorm, lwd = 2 )  
    
# Runing bdmcmc algorithm for the above simulation data set      
bmixt.obj = bmixt( data, k = 3, iter = 1000 )
    
print( bmixt.obj ) 
}
}
