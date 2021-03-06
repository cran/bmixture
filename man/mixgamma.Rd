\name{mixgamma}
\alias{rmixgamma}
\alias{dmixgamma}

\title{ Mixture of Gamma distribution }

\description{
Random generation and density function for a finite mixture of Gamma distribution.
}

\usage{
rmixgamma( n = 10, weight = 1, alpha = 1, beta = 1 )

dmixgamma( x, weight = 1, alpha = 1, beta = 1 )
}

\arguments{
  \item{n     }{ number of observations. }
  \item{x     }{ vector of quantiles.}
  \item{weight}{ vector of probability weights, with length equal to number of components (\eqn{k}). This is assumed to sum to 1; if not, it is normalized.}
  \item{alpha }{ vector of non-negative parameters of the Gamma distribution.}
  \item{beta  }{ vector of non-negative parameters of the Gamma distribution.}
}

\details{
Sampling from finite mixture of Gamma distribution, with density:
  \deqn{Pr(x|\underline{w}, \underline{\alpha}, \underline{\beta}) = \sum_{i=1}^{k} w_{i} Gamma(x|\alpha_{i}, \beta_{i}),}
where
\deqn{Gamma(x|\alpha_{i}, \beta_{i})=\frac{(\beta_{i})^{\alpha_{i}}}{\Gamma(\alpha_{i})} x^{\alpha_{i}-1} e^{-\beta_{i}x}.}
}

\value{	Generated data as an vector with size \eqn{n}. }

\references{
Mohammadi, A., Salehi-Rad, M. R., and Wit, E. C. (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics}, 28(2):683-700, \doi{10.1007/s00180-012-0323-3}

Mohammadi, A., and Salehi-Rad, M. R. (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation}, 41(3):419-435, \doi{10.1080/03610918.2011.588358}
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{rgamma}}, \code{\link{rmixnorm}}, \code{\link{rmixt}} }

\examples{
\dontrun{
n      = 10000   
weight = c( 0.6  , 0.3  , 0.1   )
alpha  = c( 100  , 200  , 300   )
beta   = c( 100/3, 200/4, 300/5 )
    
data = rmixgamma( n = n, weight = weight, alpha = alpha, beta = beta )
  
hist( data, prob = TRUE, nclass = 30, col = "gray" )
  
x            = seq( -20, 20, 0.05 )
densmixgamma = dmixnorm( x, weight, alpha, beta )
      
lines( x, densmixgamma, lwd = 2 )
}
}

\keyword{ distribution }

