\name{id}
\alias{id}
\alias{identity}
\alias{is.id}
\alias{as.id}
\alias{is.identity}
\alias{as.identity}
\alias{rsm_id}
\alias{chm_id}
\alias{qhm_id}
\alias{spin_id}
\alias{albert_id}
\alias{one}
\alias{as.one}
\title{Multiplicative identities}
\description{
Multiplying a jordan object by the \dfn{identity} leaves it unchanged.
}
\usage{
as.identity(x)
rsm_id(n,d)
chm_id(n,d)
qhm_id(n,d)
albert_id(n)
spin_id(n=3,d=5)
}
\arguments{
  \item{n}{Length of vector to be created}
  \item{d}{Dimensionality}
  \item{x}{In function \code{as.identity()}, a jordan object.  Return
    value will be a jordan object of the same dimensionality but entries
    equal to the identity}
}
\details{

  The identity object in the matrix-based classes (\code{jordan_matrix})
  is simply the identity matrix.  Function \code{as.identity()} takes an
  object of any of the five types (\code{rsm}, \code{chm}, \code{qhm},
  \code{spin}, or \code{albert}) and returns a vector of the same length
  and type, but comprising identity elements.

  Class \code{spin} has identity
  \eqn{\left(1,\mathbf{0}\right)}{(1,0)}.

}
\value{A jordan object is returned.}
\author{Robin K. S. Hankin}
\examples{

x <- as.albert(matrix(sample(1:99,81,replace=TRUE),nrow=27))
I <- as.identity(x)
x == x*I  # should be TRUE


rsm_id(6,3)




}
