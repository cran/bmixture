\name{mixt}
\alias{rmixt}
\alias{dmixt}

\title{ Mixture of t-distribution }

\description{
Random generation and density function for a finite mixture of univariate t-distribution.
}

\usage{ 
rmixt( n = 10, weight = 1, df = 1, mean = 0, sd = 1 ) 

dmixt( x, weight = 1, df = 1, mean = 0, sd = 1 ) 
}

\arguments{
  \item{n     }{ number of observations.}
  \item{x     }{ vector of quantiles.}
  \item{weight}{ vector of probability weights, with length equal to number of components (\eqn{k}). This is assumed to sum to 1; if not, it is normalized.}
  \item{df    }{ vector of degrees of freedom (> 0, maybe non-integer). df = Inf is allowed. }
  \item{mean  }{ vector of means.}
  \item{sd    }{ vector of standard deviations.}
}

\details{
Sampling from finite mixture of t-distribution, with density:

  \deqn{Pr(x|\underline{w}, \underline{df}, \underline{\mu}, \underline{\sigma}) = \sum_{i=1}^{k} w_{i} t_{df}(x| \mu_{i}, \sigma_{i}),}
where 
  \deqn{ t_{df}(x| \mu, \sigma) = \frac{ \Gamma( \frac{df+1}{2} ) }{ \Gamma( \frac{df}{2} ) \sqrt{\pi df} \sigma } \left( 1 + \frac{1}{df} \left( \frac{x-\mu}{\sigma} \right) ^2 \right) ^{- \frac{df+1}{2} }.}

}

\value{	Generated data as an vector with size \eqn{n}. }

\references{
Mohammadi, A., Salehi-Rad, M. R., and Wit, E. C. (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics}, 28(2):683-700, \doi{10.1007/s00180-012-0323-3}

Mohammadi, A., and Salehi-Rad, M. R. (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation}, 41(3):419-435, \doi{10.1080/03610918.2011.588358}
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{rt}}, \code{\link{rmixnorm}}, \code{\link{rmixgamma}} }

\examples{
\dontrun{
n      = 10000   
weight = c( 0.3, 0.5, 0.2 )
df     = c( 4  , 4  , 4   )
mean   = c( 0  , 10 , 3   )
sd     = c( 1  , 1  , 1   )
    
data = rmixt( n = n, weight = weight, df = df, mean = mean, sd = sd )
  
hist( data, prob = TRUE, nclass = 30, col = "gray" )
  
x           = seq( -20, 20, 0.05 )
densmixt = dmixt( x, weight, df, mean, sd )
      
lines( x, densmixt, lwd = 2 )
}
}

\keyword{ distribution }

