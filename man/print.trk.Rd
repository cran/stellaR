\name{print.trk}
\alias{print.trk}
\alias{print.hb}
\alias{print.zahb}
\alias{print.iso}
\title{Print stellar track objects}
\description{
Format and print the contents of the objects of classes \code{trk},
\code{hb}, \code{zahb}, \code{iso}. 
 }
\usage{
\method{print}{trk}(x, \ldots)
\method{print}{hb}(x, \ldots)
\method{print}{zahb}(x, \ldots)
\method{print}{iso}(x, \ldots)
     }
\arguments{
\item{x}{An object of classes \code{trk},
  \code{hb}, \code{zahb}, or \code{iso}.}
\item{\dots}{further arguments passed to low level print functions.}
}

\details{
The contents of the variables \code{mass} (or \code{age}), \code{z}, \code{y},
\code{ml}, \code{alpha.enh} (and \code{massRGB} if appliable)
from classes \code{trk},
\code{hb}, \code{zahb}, or \code{iso}
are fomatted and printed. The function returns its argument invisibly (via invisible(x)).
}
\examples{

  trk <- getTrk(0.9, 0.002, 0.25, 1.7, 0)
  trk

}
\keyword{manip}
