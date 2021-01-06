\name{valid}
\alias{valid}
\alias{validity}
\alias{valid_rsm}
\alias{valid_chm}
\alias{valid_qhm}
\alias{is_ok_rsm}
\alias{is_ok_chm}
\alias{is_ok_qhm}
\alias{is_ok_albert}
\alias{is_ok_rsm}
\alias{valid_albert}
\title{Validity methods}
\description{Validity methods, to check that objects are well-formed}
\usage{
valid_rsm(object)
valid_chm(object)
valid_qhm(object)
valid_albert(object)
is_ok_rsm(r)
is_ok_chm(r)
is_ok_qhm(r)
is_ok_albert(r)
is_ok_rsm(r)
}
\arguments{
  \item{object}{Putative jordan object}
  \item{r}{Integer, number of rows in putative jordan object}
}
\details{
Validity methods.  The \code{validity_foo()} functions test for an object
to be the right type, and the  \code{is_ok_foo()} functions test the
number of rows being appropriate for a jordan object of some type; these
functions return an error if not appropriate, or, for
\code{jordan_matrix} objects, the size of the matrix worked with.
}
\value{Return a Boolean}
\author{Robin K. S. Hankin}
\examples{
 is_ok_qhm(45)   # 5x5 Hermitian quaternionic matrices
#is_ok_qhm(46)   # FALSE
}