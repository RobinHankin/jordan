\name{jordan}
\alias{jordan}
\alias{is.rsm}
\alias{is.real_symmetric_matrix}
\alias{is.quaternion_herm_matrix}
\alias{quaternion_herm_matrix}
\alias{complex_herm_matrix}
\alias{is.complex_herm_matrix}
\alias{real_symmetric_matrix}
\alias{albert}
\alias{is.albert}
\alias{is.jordan}
\alias{as.spin}
\alias{is.spin}
\alias{spin}
\alias{numeric_compare_jordan}
\title{Create jordan objects}
\description{Creation methods for jordan objects}
\arguments{
  \item{M}{A matrix with columns representing independent entries in a
    matrix-based Jordan algebra}
  \item{a,V}{Scalar and vector components of a spin factor}
}
\details{
  The functions documented here are the creation methods for the five
  types of jordan algebra.
  \itemize{
    \item\code{quaternion_herm_matrix()}
    \item\code{complex_herm_matrix()}
    \item\code{real_symmetric_matrix()}
    \item\code{albert()}
    \item\code{spin()}
}

(to generate quick \dQuote{get you going} Jordan algebra objects, use
the \code{rrsm()} family of functions, documented at \code{random.Rd}).

}
\value{Return jordans or Boolean as appropriate}
\author{Robin K. S. Hankin}
\seealso{\code{\link{random}}}
\examples{

A <- real_symmetric_matrix(1:10)  # vector of length 1
as.1matrix(A)                     # in matrix form

complex_herm_matrix(cbind(1:25,2:26))
quaternion_herm_matrix(1:15)

albert(1:27)
spin(-6,cbind(1:12,12:1))


x <- rrsm() ; y <- rrsm() ; z <- rrsm()  # also works with the other Jordans

x*(y*z) - (x*y)*z         # Jordan algebra is not associative...
(x*y)*(x*x) - x*(y*(x*x)) # but satisfies the Jordan identity


}
