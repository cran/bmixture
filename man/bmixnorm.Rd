\name{bmixnorm}
\alias{bmixnorm}

\title{	Sampling algorithm for mixture of Normal distributions }

\description{
	This function consists of several sampling algorithms for Bayesian estimation for finite a mixture of Normal distributions. 
}
 
\usage{
bmixnorm( data, k = "unknown", iter = 1000, burnin = iter / 2, lambda = 1, 
          k.start = NULL, mu.start = NULL, sig.start = NULL, pi.start = NULL, 
          k.max = 30, trace = TRUE )
}

\arguments{
	\item{data     }{ vector of data with size \code{n}.	}
	\item{k        }{ number of components of mixture distribution. It can take an integer values. }
	\item{iter     }{ number of iteration for the sampling algorithm.}
	\item{burnin   }{ number of burn-in iteration for the sampling algorithm.}
	\item{lambda   }{ For the case \code{k = "unknown"}, it is the parameter of the prior distribution of number of components \code{k}. }
	\item{k.start  }{ For the case \code{k = "unknown"}, initial value for number of components of mixture distribution.}
	\item{mu.start }{ Initial value for parameter of mixture distribution.}
	\item{sig.start}{ Initial value for parameter of mixture distribution.}
	\item{pi.start }{ Initial value for parameter of mixture distribution.}
	\item{k.max    }{ For the case \code{k = "unknown"}, maximum value for the number of components of mixture distribution.}
	\item{trace    }{ Logical: if TRUE (default), tracing information is printed.}
}

\details{
Sampling from finite mixture of Normal distribution, with density:

  \deqn{Pr(x|k, \underline{\pi}, \underline{\mu}, \underline{\sigma}) = \sum_{i=1}^{k} \pi_{i} N(x|\mu_{i}, \sigma^2_{i}),}
  
where \code{k} is the number of components of mixture distribution (as a defult we assume is \code{unknown}).
The prior distributions are defined as below

\deqn{ P(K=k) \propto \frac{\lambda^k}{k!}, \ \ \ k=1,...,k_{max},}
\deqn{ \pi_{i} | k \sim Dirichlet( 1,..., 1 ), }
\deqn{ \mu_{i} | k  \sim N( \epsilon, \kappa ), }
\deqn{ \sigma_i | k \sim IG( g, h ), }

where \code{IG} denotes an inverted gamma distribution. For more details see for more details see Stephens, M. (2000), \doi{10.1214/aos/1016120364}.
}

\value{
	An object with \code{S3} class \code{"bmixnorm"} is returned:

	\item{all_k      }{ a vector which includes the waiting times for all iterations. It is needed for monitoring the convergence of the BD-MCMC algorithm. }
	\item{all_weights}{ a vector which includes the waiting times for all iterations. It is needed for monitoring the convergence of the BD-MCMC algorithm. }
	\item{pi_sample  }{ a vector which includes the MCMC samples after burn-in from parameter \code{pi} of mixture distribution. }
	\item{mu_sample  }{ a vector which includes the MCMC samples after burn-in from parameter \code{mu} of mixture distribution. }
	\item{sig_sample }{ a vector which includes the MCMC samples after burn-in from parameter \code{sig} of mixture distribution. }
	\item{data       }{ original data. }
}

\references{
Stephens, M. (2000) Bayesian analysis of mixture models with an unknown number of components-an alternative to reversible jump methods. \emph{Annals of statistics}, 28(1):40-74, \doi{10.1214/aos/1016120364}

Richardson, S. and Green, P. J. (1997) On Bayesian analysis of mixtures with an unknown number of components. \emph{Journal of the Royal Statistical Society: series B}, 59(4):731-792, \doi{10.1111/1467-9868.00095}

Green, P. J. (1995) Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. \emph{Biometrika}, 82(4):711-732, \doi{10.1093/biomet/82.4.711}

Cappe, O., Christian P. R., and Tobias, R. (2003) Reversible jump, birth and death and more general continuous time Markov chain Monte Carlo samplers. \emph{Journal of the Royal Statistical Society: Series B}, 65(3):679-700

Mohammadi, A., Salehi-Rad, M. R., and Wit, E. C. (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics}, 28(2):683-700, \doi{10.1007/s00180-012-0323-3}

Mohammadi, A., and Salehi-Rad, M. R. (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation}, 41(3):419-435, \doi{10.1080/03610918.2011.588358}

Wade, S. and Ghahramani, Z. (2018) Bayesian Cluster Analysis: Point Estimation and Credible Balls (with Discussion). \emph{Bayesian Analysis}, 13(2):559-626, \doi{10.1214/17-BA1073}
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{bmixt}}, \code{\link{bmixgamma}}, \code{\link{rmixnorm}} }

\examples{
\dontrun{
data( galaxy )

set.seed( 70 )
# Runing bdmcmc algorithm for the galaxy dataset      
mcmc_sample = bmixnorm( data = galaxy )

summary( mcmc_sample ) 
plot( mcmc_sample )
print( mcmc_sample)

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
bmixnorm.obj = bmixnorm( data, k = 3, iter = 1000 )
    
summary( bmixnorm.obj ) 
}
}

\keyword{sampling algorithms}
\keyword{iteration}

