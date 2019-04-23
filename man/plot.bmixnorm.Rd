\name{plot.bmixnorm}
\alias{plot.bmixnorm}

\title{ Plot function for \code{S3} class \code{"bmixnorm"} }

\description{ Visualizes the results for function \code{\link{bmixnorm}}. }

\usage{\method{plot}{bmixnorm}( x, ... )}

\arguments{
  \item{x    }{ An object of \code{S3} class \code{"bmixnorm"}, from function \code{\link{bmixnorm}}. }
  \item{\dots}{ System reserved (no specific usage). }
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{bmixnorm}} }

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
bmixnorm.obj = bmixnorm( data, k = 3 )
    
plot( bmixnorm.obj ) 
}
}
