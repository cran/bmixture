\name{mixnorm}
\alias{rmixnorm}
\alias{dmixnorm}

\title{ Mixture of Normal distribution }

\description{
Random generation and density function for the finite mixture of univariate Normal distribution.
}

\usage{ 
rmixnorm( n = 10, weight = 1, mean = 0, sd = 1 ) 

dmixnorm( x, weight = 1, mean = 0, sd = 1 ) 
}

\arguments{
  \item{n     }{ The number of samples required.}
  \item{x     }{ The vector of quantiles.}
  \item{weight}{ The vector of probability weights, with length equal to number of components (\eqn{k}). This is assumed to sum to 1; if not, it is normalized.}
  \item{mean  }{ The vector of means.}
  \item{sd    }{ The vector of standard deviations.}
}

\details{
Sampling from finite mixture of Normal distribution, with density:
  \deqn{Pr(x|\underline{w}, \underline{\mu}, \underline{\sigma}) = \sum_{i=1}^{k} w_{i} N(x|\mu_{i}, \sigma_{i}).}
}

\value{	Generated data as an vector with size \eqn{n}. }

\references{
Mohammadi, A., Salehi-Rad, M. R., and Wit, E. C. (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics}, 28(2):683-700

Mohammadi, A., and Salehi-Rad, M. R. (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation}, 41(3):419-435
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{rmixt}}, \code{\link{rmixgamma}} }

\examples{
\dontrun{
n      = 10000   
weight = c( 0.3, 0.5, 0.2 )
mean   = c( 0  , 10 , 3   )
sd     = c( 1  , 1  , 1   )
    
data = rmixnorm( n = n, weight = weight, mean = mean, sd = sd )
  
hist( data, prob = TRUE, nclass = 30, col = "gray" )
  
x           = seq( -20, 20, 0.05 )
densmixnorm = dmixnorm( x, weight, mean, sd )
      
lines( x, densmixnorm, lwd = 2 )
}
}

\keyword{ distribution }
