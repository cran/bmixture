\name{bmixt}
\alias{bmixt}

\title{	Sampling algorithm for mixture of distributions }

\description{
	This function consists of several sampling algorithms for Bayesian estimation for finite mixture of Normal distributions. 
}
 
\usage{
bmixt( data, k = "unknown", iter = 1000, burnin = iter / 2, lambda = 1, df = 1, 
       k.start = NULL, mu.start = NULL, sig.start = NULL, pi.start = NULL, 
       k_max = 30, trace = TRUE )
}

\arguments{
	\item{data}{ The vector of data with size \code{n}.	}
	\item{k}{ The number of components of mixture distribution. Defult is \code{"unknown"}. It can take an integer values. }
	\item{iter}{The number of iteration for the sampling algorithm.}
	\item{burnin}{The number of burn-in iteration for the sampling algorithm.}
	\item{lambda}{For the case \code{k = "unknown"} and it is rate for birth and parameter of prior distribution of number of components \code{k}.}
	\item{df}{ Degrees of freedom (> 0, maybe non-integer). df = Inf is allowed. }
	\item{k.start}{For the case \code{k = "unknown"}. Initial value for number of components of mixture distribution.}
	\item{mu.start}{Initial value for parameter of mixture distribution.}
	\item{sig.start}{Initial value for parameter of mixture distribution.}
	\item{pi.start}{Initial value for parameter of mixture distribution.}
	\item{k_max}{For the case \code{k = "unknown"}. Maximum value for the number of components of mixture distribution.}
	\item{trace}{Logical: if TRUE (default), tracing information is printed.}
}

\details{
Sampling from finite mixture of Gamma distribution, with density:

  \deqn{Pr(x|k, \underline{\pi}, \underline{\mu}, \underline{\sigma}) = \sum_{i=1}^{k} \pi_{i} N(x|\mu_{i}, \sigma_{i}),}

where \code{k} is the number of components of mixture distribution (as a defult we assume is \code{unknown}).
The prior distributions are defined as below

\deqn{ P(K=k) \propto \frac{\lambda^k}{k!}, \ \ \ k=1,...,k_{max},}
  
\deqn{ \pi_{i} | k \sim Dirichlet( 1,..., 1 ),}
  
\deqn{ \alpha_{i} | k  \sim Gamma(\nu, \upsilon),}
  
\deqn{ \beta_i | k \sim G(\eta, \tau),,}


for more details see Mohammadi et al. (2013) and Mohammadi and Salehi-Rad (2012).
}


\value{
	An object with \code{S3} class \code{"bmixt"} is returned:

	\item{all_k}{ A vector which includes the waiting times for all iterations. It is needed for monitoring the convergence of the BD-MCMC algorithm. }
	\item{all_weights}{ A vector which includes the waiting times for all iterations. It is needed for monitoring the convergence of the BD-MCMC algorithm. }
	\item{pi_sample}{ A vector which includes the MCMC samples after burn-in from parameter \code{pi} of mixture distribution. }
	\item{mu_sample}{ A vector which includes the MCMC samples after burn-in from parameter \code{mu} of mixture distribution. }
	\item{sig_sample}{ A vector which includes the MCMC samples after burn-in from parameter \code{sig} of mixture distribution. }
	\item{data}{ original data. }
}

\references{
Stephens, M. (2000) Bayesian analysis of mixture models with an unknown number of components-an alternative to reversible jump methods.  \emph{Annals of statistics} 28(1):40-74

Richardson, S. and P. J. Green (1997) On Bayesian analysis of mixtures with an unknown number of components. \emph{Journal of the Royal Statistical Society: series B} 59(4):731-792

Green, P. J.  (1995) Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. \emph{Biometrika} 82(4):711-732

Cappe, O., Christian P. R., and Tobias R. (2003) Reversible jump, birth and death and more general continuous time Markov chain Monte Carlo samplers. \emph{Journal of the Royal Statistical Society: Series B} 65(3):679-700

Mohammadi, A., M. R. Salehi-Rad, and E. C. Wit (2013) Using mixture of Gamma distributions for Bayesian analysis in an M/G/1 queue with optional second service. \emph{Computational Statistics} 28(2):683-700

Mohammadi, A. and M. R. Salehi-Rad (2012) Bayesian inference and prediction in an M/G/1 with optional second service. \emph{Communications in Statistics-Simulation and Computation} 41(3):419-435

Mohammadi, A. and E. Wit (2015) Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015) \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 
}

\author{ Abdolreza Mohammadi }

\seealso{\code{\link{rmixnorm}}}

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
    
summary( bmixt.obj ) 
}
}
