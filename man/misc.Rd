\name{misc}
\alias{misc}
\alias{harmonize_spin_numeric}
\alias{harmonize_spin_spin}
\alias{jordan_compare_numeric}
\alias{mymatrixpower}
\alias{mymatrixpower_onion}
\alias{quadraticform}
\alias{r1}
\alias{rn}
\alias{dim,spin-method}
\alias{length,jordan-method}
\alias{names,jordan-method}
\alias{names<-,jordan-method}
\alias{sum,jordan-method}
\title{Miscellaneous Jordan functionality}
\description{Miscellaneous Jordan functionality that should be
  documented somewhere}
\usage{
harmonize_spin_numeric(e1,e2)
harmonize_spin_spin(e1,e2)
}
\arguments{
  \item{e1,e2}{Objects to harmonize}
}
\details{
  Miscellaneous low-level helper functions.
  
  The harmonize functions \code{harmonize_spin_numeric()} and
  \code{harmonize_spin_spin()} work for spin objects for the
  matrix-based classes \code{emulator::harmonize_foo()} is used.

}
\value{These are mostly low-level helper functions; they not
  particularly user-friendly.  They generally return either
  numeric or Jordan objects.
}
\author{Robin K. S. Hankin}
