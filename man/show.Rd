\name{show}
\alias{show}
\alias{spin_show}
\alias{jordan_matrix_show}
\alias{albert_show}
\alias{spin_show}
\title{Print methods}
\description{Show methods, to display objects at the prompt}
\usage{
albert_show(x)
spin_show(x)
jordan_matrix_show(x)
}
\arguments{\item{x}{Jordan object}}
\details{
The matrix-based algebras use a show method that modifies the row and
column names of the underlying matrix slightly.

Spin factors are displayed with the scalar component offset from the
vector component.

Both are sensitive to the value of option \code{top_and_bottom}, a
two-element integer vector indicating the number of start lines and end
lines to print.

}
\value{Returns the argument}
\author{Robin K. S. Hankin}
\examples{
rspin()

rqhm()

rchm()
}