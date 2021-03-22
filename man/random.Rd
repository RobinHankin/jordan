\name{random}
\alias{ralbert}
\alias{rrsm}
\alias{rchm}
\alias{rqhm}
\alias{rspin}
\title{Random Jordan objects}
\description{Random jordan objects with specified properties}
\usage{
ralbert(n=3)
rrsm(n=3, d=5)
rchm(n=3, d=5)
rqhm(n=3,d=5)
rspin(n=3, d=5)
}
\arguments{
  \item{n}{Length of random object returned}
  \item{d}{Dimensionality of random object returned}
}
\details{
  These functions give a quick \dQuote{get you going} random Jordan
  object to play with.
}
\value{Return a jordan object}
\author{Robin K. S. Hankin}
\examples{
rrsm()
ralbert()
rspin()
}
