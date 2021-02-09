\name{r_to_n}
\alias{r_to_n}
\alias{r_to_n_rsm}
\alias{r_to_n_chm}
\alias{r_to_n_qhm}
\alias{r_to_n_albert}
\alias{n_to_r}
\alias{n_to_r_rsm}
\alias{n_to_r_chm}
\alias{n_to_r_qhm}
\alias{n_to_r_albert}
\title{Sizes of Matrix-based Jordan algebras}
\description{\loadmathjax
Given the number of rows in a (matrix-based) Jordan object, return the
size of the underlying associative matrix algebra
}
\usage{
r_to_n_rsm(r)
r_to_n_chm(r)
r_to_n_qhm(r)
r_to_n_albert(r=27)
n_to_r_rsm(n)
n_to_r_chm(n)
n_to_r_qhm(n)
n_to_r_albert(n=3)
}
\arguments{
  \item{n}{Integer, underlying associative algebra being matrices of
    size \mjeqn{n\times n}{n*n}}
  \item{r}{Integer, number of rows of independent representation of a
    matrix-based jordan object}
}
\details{
These functions are here for consistency, and the \code{albert} ones for
completeness.

For the record, they are:


\itemize{
  \item Real symmetric matrices, \code{rsm}, \mjseqn{r=n(n+1)/2},
  \mjeqn{n=(\sqrt{1+4r}-1)/2}{n=(sqrt(1+4r)-1)/2}
  \item Complex Hermitian matrices, \code{chm}, \mjseqn{r=n^2}, 
  \mjeqn{n=\sqrt{r}}{n=sqrt(r)}
  \item Quaternion Hermitian matrices, \code{qhm}, \mjseqn{r=n(2n-1)}, 
  \mjeqn{n=\sqrt{1+8r}/4}{n=sqrt(1+8r)/4}
  \item Albert algebras, \mjseqn{r=27}, \mjseqn{n=3}
}

}

\value{Return non-negative integers}
\author{Robin K. S. Hankin}
\note{
I have not been entirely consistent in my use of these functions.
}


\examples{
r_to_n_qhm(nrow(rqhm()))
}