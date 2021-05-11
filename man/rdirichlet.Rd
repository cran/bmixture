\name{rdirichlet}
\alias{rdirichlet}

\title{ Random generation for the Dirichlet distribution }

\description{ Random generation from the Dirichlet distribution. }

\usage{ rdirichlet( n = 10, alpha = c( 1, 1 ) ) }

\arguments{
  \item{n    }{ number of observations. }
  \item{alpha}{ vector of shape parameters. }
}
\details{ The Dirichlet distribution is the multidimensional generalization of the beta distribution. }

\value{ A matrix with \code{n} rows, each containing a single Dirichlet random deviate. }

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}}
  
%\seealso{ \code{\link[stats]{Beta}} }
\seealso{ \code{\link{Beta}} }

\examples{
  draws = rdirichlet( n = 500, alpha = c( 1, 1, 1 ) )
  
  boxplot( draws )
}

\keyword{distribution}
