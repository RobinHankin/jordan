\name{zero}
\alias{zero}
\alias{is_zero_jordan}
\alias{is.zero,jordan-method}
\alias{is.zero}

\title{The zero Jordan object}
\description{Package idiom for the zero Jordan object, and testing}
\usage{
is.zero(x)
is_zero_jordan(e1,e2=0)
}
\arguments{
  \item{x,e1}{Jordan object to test for zeroness}
  \item{e2}{Dummy numeric object to make the \code{Arith} method work}
}
\details{

One often wants to test a jordan object for being zero, and natural
idiom would be \code{rchm()==0}.  The helper function is
\code{is_zero_jordan()}, and the generic is \code{is.zero()}.

}
\value{Returns a Boolean}
\author{Robin K. S. Hankin}
\examples{
rrsm()*0 == 0
}
