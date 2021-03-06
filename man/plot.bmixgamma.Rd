\name{plot.bmixgamma}
\alias{plot.bmixgamma}

\title{ Plot function for \code{S3} class \code{"bmixgamma"} }

\description{ Visualizes the results for function \code{\link{bmixgamma}}. }

\usage{\method{plot}{bmixgamma}( x, ... )}

\arguments{
  \item{x    }{ An object of \code{S3} class \code{"bmixgamma"}, from function \code{\link{bmixgamma}}. }
  \item{\dots}{ System reserved (no specific usage). }
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{bmixgamma}} }

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
bmixgamma.obj <- bmixgamma( data )
    
plot( bmixgamma.obj )
}
}

\keyword{hplot}

