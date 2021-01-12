\name{jordan-package}
\alias{jordan-package}
\docType{package}
\title{
\packageTitle{jordan}
}
\description{
\packageDescription{jordan}
}
\details{\loadmathjax

  A \dfn{Jordan algebra} is a nonassociative algebra over the reals with
  a multiplication that satisfies the following identities:

  \mjsdeqn{xy=yx}

  \mjsdeqn{(xy)(xx) = x(y(xx))}

  (the second identity is known as the Jordan identity).  In literature
  one usually indicates multiplication by juxtaposition but one
  sometimes sees \mjeqn{x\circ y}{x o y}.  Package idiom is to use an
  asterisk, as in \code{x*y}.  There are five types of Jordan algebras:

  \enumerate{
    \item Real symmetric matrices, class \code{real_symmetric_matrix},
    abbreviated in the package to \code{rsm}
    \item Complex Hermitian matrices, class \code{complex_herm_matrix},
    abbreviated to \code{chm}
    \item Quaternionic Hermitian matrices, class
    \code{quaternion_herm_matrix}, abbreviated to \code{qhm}
    \item Albert algebras, the space of \mjeqn{3\times 3}{3*3}
    octonionic matrices, class \code{albert}
    \item Spin factors, class \code{spin}
  }

  (of course, the first two are special cases of the next).  The
  \code{jordan} package provides functionality to manipulate jordan
  objects using natural R idiom.

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
  \mjeqn{n\times n}{n*n} has only \mjseqn{n(n+1)/2} independent entries,
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
x <- rqhm()  # random Hermitian Quaternion algebra
y <- rqhm()  # 

x

LHS <- (x*y)*(x*x)
RHS <- x*(y*(x*x))

diff <- LHS-RHS  # zero to numerical precision


diff[1,,drop=TRUE]


}