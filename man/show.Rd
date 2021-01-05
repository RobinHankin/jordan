\name{show}
\alias{show}
\alias{spin_show}
\alias{albert_show}
\alias{spin_show}
\title{Print methods}
\description{Show methods, to display objects at the prompt}
\usage{
albert_show(x)
}
\arguments{
  \item{x}{Albert object}
}
\details{
The matrix-based algebras use a show method that modifies the row and
column names of the underlying matrix slightly.

Spin factors are displayed with the scalar component offset from the
vector component.

}
\value{Returns the argument}
\author{Robin K. S. Hankin}
\examples{
rspin()

rqhm()

rchm()
}