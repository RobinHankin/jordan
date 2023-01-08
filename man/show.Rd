\name{show}
\alias{show}
\alias{top_and_bottom}
\alias{head_and_tail}
\alias{spin_show}
\alias{description}
\alias{jordan_matrix_show}
\alias{albert_show}
\alias{spin_show}
\alias{description}
\title{Print methods}
\description{Show methods, to display objects at the prompt}
\usage{
albert_show(x)
spin_show(x)
jordan_matrix_show(x)
description(x,plural=FALSE)
}
\arguments{
  \item{x}{Jordan object}
  \item{plural}{Boolean, indicating whether plural form is to
    be given}
}
\details{

  The special algebras use a bespoke show method,
  \code{jordan_matrix_show()} or \code{spin_show()}.  If the number of
  elements is small, they display a concise representation and modify
  the row and column names of the underlying matrix slightly; spin
  factors are displayed with the scalar component offset from the vector
  component.

  Print methods for special algebras are sensitive to the value of
  option \code{head_and_tail}, a two-element integer vector indicating
  the number of start lines and end lines to print.

  Function \code{description()} gives a natural-language description of
  its argument, used in the print method.

}
\value{Returns the argument}
\author{Robin K. S. Hankin}
\examples{
rspin()

rqhm()

rchm()
}
