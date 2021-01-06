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
\description{
  The functions documented here are the creation methods for the five
  types of jordan algebra
}
\usage{
real_symmetric_matrix(M)
complex_herm_matrix(M)
albert(M)
quaternion_herm_matrix(M)
spin(a,V)
}
\arguments{
  \item{M}{A matrix with columns representing independent entries in a
    matrix-based Jordan algebra}
  \item{a,V}{Scalar and vector components of a spin factor}
}
\details{
Details here
}
\value{Return jordans or Boolean as appropriate}
\author{Robin K. S. Hankin}
\examples{
4+5
}
