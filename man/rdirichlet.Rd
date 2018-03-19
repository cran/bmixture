\name{rdirichlet}
\alias{rdirichlet}

\title{ Random generation for the Dirichlet distribution }

\description{ Random generation from the Dirichlet distribution. }

\usage{ rdirichlet( n = 10, alpha = c( 1, 1 ) ) }

\arguments{
  \item{n}{ The number of samples required. }
  \item{alpha}{Vector of shape parameters.}
}
\details{ The Dirichlet distribution is the multidimensional generalization of the beta distribution. }

\value{ A matrix with \code{n} rows, each containing a single Dirichlet random deviate. }

\author{ Reza Mohammadi }
  
\seealso{ \code{\link[stats]{Beta}} }

\examples{
  draws = rdirichlet( n = 20, alpha = c( 1, 1, 1 ) )
}

\keyword{distribution}
