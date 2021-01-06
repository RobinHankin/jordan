\name{Arith}
\docType{methods}
\alias{Arith}
\alias{Arith-methods}
\alias{Arith,jordan,missing-method}
\alias{Arith,jordan,ANY-method}
\alias{Arith,ANY,jordan-method}
\alias{Arith,jordan,jordan-method}
\alias{Arith,jordan,numeric-method}
\alias{jordan_negative}
\alias{jordan_inverse}
\alias{jordan_arith_jordan}
\alias{jordan_arith_numeric}
\alias{numeric_arith_jordan}
\alias{harmonize_oo}
\alias{harmonize_on}
\alias{jordan_plus_jordan}
\alias{jordan_plus_numeric}
\alias{jordan_prod_jordan}
\alias{octjordan_prod_octjordan}
\alias{quaternion_prod_quaternion}
\alias{jordan_prod_numeric}
\alias{jordan_prod_jordan}
\alias{jordan_power_jordan}
\alias{jordan_power_singleinteger}
\alias{jordan_power_numeric}
\alias{albert_arith_albert}
\alias{albert_arith_numeric}
\alias{albert_inverse}
\alias{albert_power_albert}
\alias{albert_power_numeric}
\alias{albert_power_single_n}
\alias{albert_prod_albert}
\alias{chm_arith_chm}
\alias{chm_arith_numeric}
\alias{chm_inverse}
\alias{chm_power_numeric}
\alias{chm_prod_chm}
\alias{numeric_arith_albert}
\alias{numeric_arith_chm}
\alias{numeric_arith_qhm}
\alias{numeric_arith_rsm}
\alias{qhm_arith_numeric}
\alias{qhm_arith_qhm}
\alias{qhm_inverse}
\alias{qhm_power_numeric}
\alias{qhm_prod_qhm}
\alias{rsm_arith_numeric}
\alias{rsm_arith_rsm}
\alias{rsm_inverse}
\alias{rsm_power_numeric}
\alias{rsm_prod_rsm}
\alias{spin_plus_numeric}
\alias{spin_plus_spin}
\alias{spin_power_numeric}
\alias{spin_power_single_n}
\alias{spin_power_spin}
\alias{spin_prod_numeric}
\alias{spin_prod_spin}
\alias{spin_inverse}
\alias{spin_negative}
\alias{vec_albertprod_vec}
\alias{vec_chmprod_vec}
\alias{vec_qhmprod_vec}
\alias{vec_rsmprod_vec}
\title{Methods for Function Arith in package Jordan}
\description{
Methods for Arithmetic functions for jordans: \code{+},
\code{-}, \code{*}, \code{/}, \code{^}
}
\usage{
jordan_negative(z)
jordan_plus_jordan(e1,e2)
jordan_plus_numeric(e1,e2)
jordan_prod_numeric(e1,e2)
jordan_power_jordan(e1,e2)
albert_arith_albert(e1,e2)
albert_arith_numeric(e1,e2)
albert_inverse(e1)
albert_power_albert(...)
albert_power_numeric(e1,e2)
albert_power_single_n(e1,n)
albert_prod_albert(e1,e2)
chm_arith_chm(e1,e2)
chm_arith_numeric(e1,e2)
chm_inverse(e1)
chm_power_numeric(e1,e2)
chm_prod_chm(e1,e2)
numeric_arith_albert(e1,e2)
numeric_arith_chm(e1,e2)
numeric_arith_qhm(e1,e2)
numeric_arith_rsm(e1,e2)
qhm_arith_numeric(e1,e2)
qhm_arith_qhm(e1,e2)
qhm_inverse(x)
qhm_power_numeric(e1,e2)
qhm_prod_qhm(e1,e2)
rsm_arith_numeric(e1,e2)
rsm_arith_rsm(e1,e2)
rsm_inverse(e1)
rsm_power_numeric(e1,e2)
rsm_prod_rsm(e1,e2)
spin_plus_numeric(e1,e2)
spin_plus_spin(e1,e2)
spin_power_numeric(e1,e2)
spin_power_single_n(e1,n)
spin_power_spin(...)
spin_prod_numeric(e1,e2)
spin_prod_spin(e1,e2)
spin_inverse(...)
spin_negative(e1)
vec_albertprod_vec(x,y)
vec_chmprod_vec(x,y)
vec_qhmprod_vec(x,y)
vec_rsmprod_vec(x,y)
}
\arguments{
  \item{z,e1,e2}{Jordan objects or numeric vectors}
  \item{n}{Integer for powers}
  \item{...}{Further arguments (ignored)}
  \item{x,y}{Numeric vectors, Jordan objects in independent form}
}
\details{

  The package implements the \code{Arith} group of \code{S4} generics so
  that idiom like \code{A + B*C} works as expected with jordans.

  Functions like \code{jordan_inverse()} and \code{jordan_plus_jordan()}
  are low-level helper functions.  The only really interesting operation
  is multiplication; functions like \code{jordan_prod_jordan()}.

  Names are implemented and the rules are inherited (via
  \code{onion::harmonize_oo()} and \code{onion::harmonize_on()}) from
  \code{rbind()}.

}
\value{generally return jordans}
\author{Robin K. S. Hankin}
\examples{

a <- rspin()
a[1] <- a[2]*7
a

}
\keyword{methods}
\keyword{math}
