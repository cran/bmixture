\name{galaxy}
\alias{galaxy}
\docType{data}
\title{ Galaxy data }
\description{
  This dataset considers of 82 observatons of the velocities (in 1000 km/second) of distant galaxies diverging from our own, from six well-separated conic sections of the Corona Borealis.
  The dataset has been analyzed under a variety of mixture models; See e.g. Stephens (2000).
}
\usage{ data( galaxy ) }
\format{
  A data frame with 82 observations on the following variable.
  \describe{
    \item{\code{speed}}{a numeric vector giving the speed of galaxies (in 1000 km/second).}
  }
}

\references{
Stephens, M. (2000) Bayesian analysis of mixture models with an unknown number of components-an alternative to reversible jump methods. \emph{Annals of statistics}, 28(1):40-74, \doi{10.1214/aos/1016120364}
}

\examples{
data( galaxy )
  
hist( galaxy, prob = TRUE, xlim = c( 0, 40 ), ylim = c( 0, 0.3 ), nclass = 20, 
      col = "gray", border = "white" )

lines( density( galaxy ), col = "black", lwd = 2 ) 
}

\keyword{ datasets }
