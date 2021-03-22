\name{extract}
\Rdversion{1.1}
\alias{extract}
\alias{index-class}
\alias{[,albert,index,missing,logical-method}
\alias{[,complex_herm_matrix,index,missing,logical-method}
\alias{[,jordan,index,ANY,ANY-method}
\alias{[,jordan,index,missing,ANY-method}
\alias{[,quaternion_herm_matrix,index,missing,logical-method}
\alias{[,real_symmetric_matrix,index,missing,logical-method}
\alias{[,spin,index,missing,ANY-method}
\alias{[,spin,missing,index,ANY-method}
\alias{[<-,albert,index,missing,albert-method}
\alias{[<-,complex_herm_matrix,index,ANY,ANY-method}
\alias{[<-,complex_herm_matrix,index,missing,complex_herm_matrix-method}
\alias{[<-,jordan_matrix,index,missing,numeric-method}
\alias{[<-,quaternion_herm_matrix,index,missing,quaternion_herm_matrix-method}
\alias{[<-,real_symmetric_matrix,index,missing,real_symmetric_matrix-method}
\alias{[<-,spin,index,index,ANY-method}
\alias{[<-,spin,index,missing,numeric-method}
\alias{[<-,spin,index,missing,spin-method}
\alias{[,spin,missing,missing,ANY-method}
\alias{[<-,spin,ANY,missing,ANY-method}
\alias{[<-,spin,missing,ANY,numeric-method}
\alias{[<-,spin,missing,missing,numeric-method}
\alias{[<-,spin,missing,missing,spin-method}
\title{Extract and replace methods for jordan objects}
\description{Extraction and replace methods for jordan objects should
  work as expected.

  Replace methods can take a jordan or a numeric, but  the numeric must
  be zero.
  
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "albert", i = "index", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "complex_herm_matrix", i = "index", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "jordan", i = "index", j = "ANY", drop = "ANY")}: ... }
    \item{[}{\code{signature(x = "jordan", i = "index", j = "missing", drop = "ANY")}: ... }
    \item{[}{\code{signature(x = "quaternion_herm_matrix", i = "index", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "real_symmetric_matrix", i = "index", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "spin", i = "index", j = "missing", drop = "ANY")}: ... }
    \item{[}{\code{signature(x = "spin", i = "missing", j = "index", drop = "ANY")}: ... }
    \item{[<-}{\code{signature(x = "albert", i = "index", j = "missing", value = "albert")}: ... }
    \item{[<-}{\code{signature(x = "complex_herm_matrix", i = "index", j = "ANY", value = "ANY")}: ... }
    \item{[<-}{\code{signature(x = "complex_herm_matrix", i = "index", j = "missing", value = "complex_herm_matrix")}: ... }
    \item{[<-}{\code{signature(x = "jordan_matrix", i = "index", j = "missing", value = "numeric")}: ... }
    \item{[<-}{\code{signature(x = "quaternion_herm_matrix", i = "index", j = "missing", value = "quaternion_herm_matrix")}: ... }
    \item{[<-}{\code{signature(x = "real_symmetric_matrix", i = "index", j = "missing", value = "real_symmetric_matrix")}: ... }
    \item{[<-}{\code{signature(x = "spin", i = "index", j = "index", value = "ANY")}: ... }
    \item{[<-}{\code{signature(x = "spin", i = "index", j = "missing", value = "numeric")}: ... }
    \item{[<-}{\code{signature(x = "spin", i = "index", j = "missing", value = "spin")}: ... }
	 }
}
\author{Robin K. S. Hankin}
\value{Generally return a jordan object of the same class as the first argument}
\examples{
showClass("index")  # taken from the Matrix package

a <- rspin(7)
a[2:4] <- 0
a[5:7] <- a[1]*10
a

}
\keyword{classes}
