\name{bmixgamma}
\alias{bmixgamma}

\title{	Sampling algorithm for mixture of distributions }

\description{
	This function consists of several sampling algorithms for Bayesian estimation for finite mixture of Gamma distributions. 
}
 
\usage{
bmixgamma( data, k = "unknown", iter = 1000, burnin = iter / 2, lambda = 1, 
           mu = NULL, nu = NULL, kesi = NULL, tau = NULL, k.start = NULL, 
           alpha.start = NULL, beta.start = NULL, pi.start = NULL, 
           k_max = 30, trace = TRUE )
}

\arguments{
	\item{data       }{ The vector of data with size \code{n}.	}
	\item{k          }{ The number of components of mixture distribution. Defult is \code{"unknown"}. It can take an integer values. }
	\item{iter       }{ The number of iteration for the sampling algorithm.}
	\item{burnin     }{ The number of burn-in iteration for the sampling algorithm.}
	\item{lambda     }{ For the case \code{k = "unknown"}, it is the parameter of the prior distribution of number of components \code{k}. }
	\item{mu         }{ The parameter of alpha in mixture distribution.}
	\item{nu         }{ The parameter of alpha in mixture distribution.}
	\item{kesi       }{ The parameter of beta in mixture distribution.}
	\item{tau        }{ The parameter of beta in mixture distribution.}
	\item{k.start    }{ For the case \code{k = "unknown"}, initial value for number of components of mixture distribution.}
	\item{alpha.start}{ Initial value for parameter of mixture distribution.}
	\item{beta.start }{ Initial value for parameter of mixture distribution.}
	\item{pi.start   }{ Initial value for parameter of mixture distribution.}
	\item{k_max      }{ For the case \code{k = "unknown"}, maximum value for the number of components of mixture distribution.}
	\item{trace      }{ Logical: if TRUE (default), tracing information is printed.}
}

\details{
Sampling from finite mixture of Gamma distribution, with density:
  \deqn{Pr(x|k, \underline{\pi}, \underline{\alpha}, \underline{\beta}) = \sum_{i=1}^{k} \pi_{i} Gamma(x|\alpha_{i}, \beta_{i}),}
where \code{k} is the number of components of mixture distribution (as a defult we assume is \code{unknown}) and  
\deqn{Gamma(x|\alpha_{i}, \beta_{i})=\frac{(\beta_{i})^{\alpha_{i}}}{\Gamma(\alpha_{i})} x^{\alpha_{i}-1} e^{-\beta_{i}x}.}
The prior distributions are defined as below
\deqn{ P(K=k) \propto \frac{\lambda^k}{k!}, \ \ \ k=1,...,k_{max},}
\deqn{ \pi_{i} | k \sim Dirichlet( 1,..., 1 ),}
\deqn{ \alpha_{i} | k  \sim Gamma(\nu, \upsilon),}
\deqn{ \beta_i | k \sim G(\eta, \tau),}
for more details see Mohammadi et al. (2013).
}


\value{
	An object with \code{S3} class \code{"bmixgamma"} is returned:

	\item{all_k       }{ A vector which includes the waiting times for all iterations. It is needed for monitoring the convergence of the BD-MCMC algorithm. }
	\item{all_weights }{ A vector which includes the waiting times for all iterations. It is needed for monitoring the convergence of the BD-MCMC algorithm. }
	\item{pi_sample   }{ A vector which includes the MCMC samples after burn-in from parameter \code{pi} of mixture distribution. }
	\item{alpha_sample}{ A vector which includes the MCMC samples after burn-in from parameter \code{alpha} of mixture distribution. }
	\item{beta_sample }{ A vector which includes the MCMC samples after burn-in from parameter \code{beta} of mixture distribution. }
	\item{data        }{ original data. }
}

\references{
Mohammadi, A., Salehi-Rad, M. R., and Wit, E. C. (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics}, 28(2):683-700

Mohammadi, A. and Salehi-Rad, M. R. (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation}, 41(3):419-435

Stephens, M. (2000) Bayesian analysis of mixture models with an unknown number of components-an alternative to reversible jump methods. \emph{Annals of statistics}, 28(1):40-74

Richardson, S. and Green, P. J. (1997) On Bayesian analysis of mixtures with an unknown number of components. \emph{Journal of the Royal Statistical Society: series B}, 59(4):731-792

Green, P. J. (1995) Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. \emph{Biometrika}, 82(4):711-732

Cappe, O., Christian P. R., and Tobias, R. (2003) Reversible jump, birth and death and more general continuous time Markov chain Monte Carlo samplers. \emph{Journal of the Royal Statistical Society: Series B}, 65(3):679-700

Wade, S. and Ghahramani, Z. (2018) Bayesian Cluster Analysis: Point Estimation and Credible Balls (with Discussion). \emph{Bayesian Analysis}, 13(2):559-626
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{bmixnorm}}, \code{\link{bmixt}}, \code{\link{bmixgamma}} }

\examples{
\dontrun{
# simulating data from mixture of gamma with two components
n      = 1000 # number of observations
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
bmixgamma.obj <- bmixgamma( data, iter = 1000 )
    
summary( bmixgamma.obj )  
   
plot( bmixgamma.obj )        
}
}
