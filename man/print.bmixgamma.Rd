\name{print.bmixgamma}
\alias{print.bmixgamma}

\title{Print function for \code{S3} class \code{"bmixgamma"} }

\description{ Prints the information about the output of function \code{\link{bmixgamma}}. }

\usage{\method{print}{bmixgamma}( x, ... )}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"bmixgamma"}, from function \code{\link{bmixgamma}}. }
  \item{\dots}{System reserved (no specific usage). }
}

\references{
Mohammadi, A., M. R. Salehi-Rad, and E. C. Wit (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics} 28(2):683-700

Mohammadi, A., and M. R. Salehi-Rad (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation} 41(3):419-435
}

\author{ Reza Mohammadi }

\seealso{\code{\link{bmixgamma}}}

\examples{
\dontrun{
# simulating data from mixture of gamma with two components
n      = 500 # number of observations
weight = c( 0.6, 0.4 )
alpha  = c( 12 , 1   )
beta   = c( 3  , 2   )

data <- rmixgamma( n = n, weight = weight, alpha = alpha, beta = beta )
  
# plot for simulation data    
hist( data, prob = TRUE, nclass = 50, col = "gray" )
  
x     = seq( 0, 10, 0.05 )
truth = dmixgamma( x, weight, alpha, beta )
      
lines( x, truth, lwd = 2 )
  
# Runing bdmcmc algorithm for the above simulation data set      
bmixgamma.obj <- bmixgamma( data, iter = 500 )
    
print( bmixgamma.obj )
}
}
