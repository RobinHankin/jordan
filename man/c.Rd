\name{c}
\alias{c}
\alias{c.jordan}
\alias{concatenate}
\alias{conc}
\alias{conc_pair}
\alias{c,jordan-method}
\title{Concatenation}
\description{
Combines its arguments to form a single jordan object.
}
\usage{
\S4method{c}{jordan}(x,...)
}
\arguments{\item{x,\dots}{Jordan objects}}
\details{
  Returns a concatenated jordan of the same type as its arguments.
  Argument checking is not performed.
}
\value{Returns a Jordan object of the appropriate type (coercion is not performed)}
\author{Robin K. S. Hankin}
\note{
  
  Names are inherited from the behaviour of \code{cbind()}, not
  \code{c()}.
}

\examples{
c(rqhm(),rqhm()*10)
}
