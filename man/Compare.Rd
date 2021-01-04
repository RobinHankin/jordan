\name{Compare-methods}
\docType{methods}
\alias{Compare-methods}
\alias{Compare,jordan,jordan-method}
\alias{Compare,jordan,numeric-method}
\alias{Compare,numeric,jordan-method}
\alias{Compare,jordan,jordan-method}
\alias{Compare,jordan,ANY-method}
\alias{Compare,ANY,jordan-method}
\alias{jordan_compare}
\alias{jordan_equal_jordan}
\alias{jordan_equal_single}
\alias{jordan_compare_jordan}
\alias{jordan_compare_single}
\alias{spin_equal_spin}
\title{Methods for compare S4 group}
\description{
  Methods for comparison (equal to, greater than, etc) of
  jordans.  Only equality makes sense.
}
\usage{jordan_compare_jordan(e1,e2)}
\arguments{\item{e1,e2}{Jordan objects}}
\value{Return a boolean}

\examples{
# rspin() > 0 # meaningless and returns an error
}
\keyword{methods}
\keyword{math}
