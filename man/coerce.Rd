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
\alias{matrix1_to_jordan}
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
matrix1_to_jordan(x)
}
\arguments{
  \item{x,e1}{Numeric vector of independent entries}
  \item{M,H}{A matrix}
  \item{d}{Dimensionality of algebra}
  \item{single}{Boolean, indicating whether a single value is to be
  returned}
  \item{class}{Class of object}
  \item{...}{Further arguments, currently ignored}
}
\details{

  The numeral \dQuote{1} in a function name means that it operates on,
  or returns, a single element, usually a matrix.  Thus function
  \code{as.1matrix()} is used to convert a jordan object to a list of
  matrices.  Length one jordan objects are converted to a matrix.
  
  Functions \code{vec_to_rsm1()} et seq convert a numeric vector to a
  (symmetric, complex, quaternion, octonion) matrix, that is, elements
  of a matrix-based Jordan algebra.

  Functions \code{rsm1_to_vec()} convert a (symmetric, complex,
  quaternion, octonion) matrix to a numeric vector of independent
  components.  The upper triangular components are used; no checking for
  symmetry is performed (the lower triangular components, and non-real
  components of the diagonal, are discarded).

  Functions \code{as.real_symmetric_matrix()},
  \code{as.complex_herm_matrix()}, \code{as.quaternion_herm_matrix()}
  and \code{as.albert()} take a numeric matrix and return a
  (matrix-based) Jordan object.

  Functions \code{numeric_to_real_symmetric_matrix()} have not been
  coded up yet.

  Function \code{matrix1_to_jordan()} takes a matrix and returns a
  length-1 (matrix based) Jordan vector.  It uses the class of the
  entries (real, complex, quaternion, octonion) to decide which type of
  Jordan to return.
  
}
\value{Return a coerced value.}
\author{Robin K. S. Hankin}
\examples{


vec_to_chm1(1:16)  # Hermitian matrix

as.1matrix(rchm())

as.complex_herm_matrix(matrix(runif(75),ncol=3))

matrix1_to_jordan(cprod(matrix(rnorm(35),7,5)))
matrix1_to_jordan(matrix(c(1,1+1i,1-1i,3),2,2))
matrix1_to_jordan(Oil + matrix(1,3,3))


}
