\name{jordan-package}
\alias{jordan-package}
\docType{package}
\title{
\packageTitle{jordan}
}
\description{
\packageDescription{jordan}
}
\details{

  A \dfn{Jordan algebra} is a non-associative algebra over the reals
  with a multiplication that satisfies the following identities:

  \deqn{xy=yx}

  \deqn{(xy)(xx) = x(y(xx))}

  (the second identity is known as the Jordan identity).  In literature
  one usually indicates multiplication by juxtaposition but one
  sometimes sees \eqn{x\circ y}{x o y}.  Package idiom is to use an
  asterisk, as in \code{x*y}.  There are five types of Jordan algebras:

  \enumerate{
    \item Real symmetric matrices, class \code{real_symmetric_matrix},
    abbreviated in the package to \code{rsm}
    \item Complex Hermitian matrices, class \code{complex_herm_matrix},
    abbreviated to \code{chm}
    \item Quaternionic Hermitian matrices, class
    \code{quaternion_herm_matrix}, abbreviated to \code{qhm}
    \item Albert algebras, the space of \eqn{3\times 3}{3*3}
    octonionic matrices, class \code{albert}
    \item Spin factors, class \code{spin}
  }

  (of course, the first two are special cases of the next).  The
  \code{jordan} package provides functionality to manipulate jordan
  objects using natural \R idiom.

  Objects of all these classes are stored in dataframe (technically, a
  matrix) form with columns being elements of the jordan algebra.

  The first four classes are matrix-based in the sense that the
  algebraic objects are symmetric or Hermitian matrices (the \code{S4}
  class is \dQuote{\code{jordan_matrix}}).  The fifth class, spin
  factors, is not matrix based.

  One can extract the symmetric or Hermitian matrix from objects of
  class \code{jordan_matrix} using \code{as.list()}, which will return a
  list of symmetric or Hermitian matrices.  A function name preceded by
  a \dQuote{1} (for example \code{as.1matrix()} or \code{vec_to_qhm1()})
  means that it deals with a single (symmetric or Hermitian) matrix.

  Algebraically, the matrix form of \code{jordan_matrix} objects is
  redundant (for example, a \code{real_symmetric_matrix} of size
  \eqn{n\times n}{n*n} has only \eqn{n(n+1)/2} independent entries,
  corresponding to the upper triangular elements).
  
}
\author{
\packageAuthor{jordan}

Maintainer: \packageMaintainer{jordan}
}
\references{
K. McCrimmon 1978. \dQuote{Jordan algebras and their applications}.
\emph{Bulletin of the American Mathematical Society}, Volume 84, Number 4.
}
\keyword{ package }
\examples{

rrsm()     # Random Real Symmetric matrices
rchm()     # Random Complex Hermitian matrices
rqhm()     # Random Quaternionic Hermitian matrices
ralbert()  # Random Albert algebra
rspin()    # Random spin factor

x <- rqhm(n=1)  
y <- rqhm(n=1)  
z <- rqhm(n=1)  

x/1.2 + 0.3*x*y     # Arithmetic works as expected ...
x*(y*z) -(x*y)*z    # ... but '*' is not associative


## Verify the Jordan identity for  type 3 algebras:

LHS <- (x*y)*(x*x)
RHS <- x*(y*(x*x))

diff <- LHS-RHS  # zero to numerical precision

diff[1,drop=TRUE]  # result in matrix form

}
