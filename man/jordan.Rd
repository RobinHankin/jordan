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

}
\value{Return jordans or Boolean as appropriate}
\author{Robin K. S. Hankin}
\examples{

A <- real_symmetric_matrix(1:10)  # vector of length 1
as.1matrix(A)                     # in matrix form

complex_herm_matrix(cbind(1:25,2:26))
quaternion_herm_matrix(1:15)

albert(1:27)
spin(-6,cbind(1:12,12:1))

}
