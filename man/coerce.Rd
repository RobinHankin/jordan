\name{coerce}
\alias{coerce}
\alias{coercion}
\alias{as.1matrix}
\alias{as.jordan}
\alias{as.1matrix,real_symmetric_matrix-method}
\alias{as.1matrix,complex_herm_matrix-method}
\alias{as.1matrix,quaternion_herm_matrix-method}
\alias{as.1matrix,albert-method}
\alias{as.1matrix,spin-method}
\alias{vec_to_rsm1}
\alias{vec_to_chm1}
\alias{vec_to_qhm1}
\alias{vec_to_albert1}
\alias{rsm1_to_vec}
\alias{chm1_to_vec}
\alias{qhm1_to_vec}
\alias{albert1_to_vec}
\alias{as.real_symmetric_matrix}
\alias{as.complex_herm_matrix}
\alias{as.quaternion_herm_matrix}
\alias{as.albert}
\alias{numeric_to_real_symmetric_matrix}
\alias{numeric_to_complex_herm_matrix}
\alias{numeric_to_quaternion_herm_matrix}
\alias{numeric_to_albert}
\alias{as.list}
\alias{as.list,albert-method}
\alias{as.matrix,jordan-method}
\title{Coercion}
\description{Various coercions needed in the package}
\usage{
as.jordan(x,class)
vec_to_rsm1(x)
vec_to_chm1(x)
vec_to_qhm1(x)
vec_to_albert1(x)
rsm1_to_vec(M)
chm1_to_vec(M)
qhm1_to_vec(M)
albert1_to_vec(H)
as.real_symmetric_matrix(x,d,single=FALSE)
as.complex_herm_matrix(x,d,single=FALSE)
as.quaternion_herm_matrix(x,d,single=FALSE)
as.albert(x,single=FALSE)
numeric_to_real_symmetric_matrix(x,d)
numeric_to_complex_herm_matrix(x,d)
numeric_to_quaternion_herm_matrix(x,d)
numeric_to_albert(e1)
as.list(x,...)
}
\arguments{
  \item{x,e1}{Numeric vector of independent entries}
  \item{M,H}{A matrix}
  \item{d}{Dimensionality of algebra}
  \item{single}{Boolean, indicating whether a single value is tobe
    returned}
  \item{class}{Class of object}
  \item{...}{Further arguments, currently ignored}
}
\details{
The prefix \dQuote{1} means a single element, usually a matrix.
}
\value{Return a coerced value.}
\author{Robin K. S. Hankin}
\examples{


vec_to_chm1(1:16)  # Hermitian matrix

as.1matrix(rchm())

as.complex_herm_matrix(matrix(runif(75),ncol=3))

matrix1_to_jordan(3)

}
