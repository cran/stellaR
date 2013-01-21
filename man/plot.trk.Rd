\name{plot.trk}
\alias{plot.trk}
\alias{plot.hb}
\alias{plot.zahb}
\alias{plot.iso}
\title{Plot stellar track or isochrones objects}
\description{
Plot one evolutionary track or one isochrone extracted from objects of
classes \code{trk}, 
\code{hb}, \code{zahb}, \code{iso}. 
 }
\usage{
\method{plot}{trk}(x, add=FALSE, \dots)
\method{plot}{hb}(x, add=FALSE, \dots)
\method{plot}{zahb}(x, add=FALSE, \dots)
\method{plot}{iso}(x, add=FALSE, \dots)
     }
\arguments{
\item{x}{an object of classes \code{trk},
  \code{hb}, \code{zahb}, or \code{iso}.}
\item{add}{a logical value indicating whether the plot should be added to
an existing device.}
\item{\dots}{further arguments passed to \code{\link{plotAstro}}.}
}

\details{
The contents of the vectors \code{x$data$logTe} and \code{x$data$logL} are
plotted by a call to \code{\link{plotAstro}}. Sensible values for the
axes major ticks and labels are chosen.
}
\seealso{
  \code{\link{plotAstro}}.
}
\examples{

  \dontrun{
  trk <- getTrk(0.9, 0.002, 0.25, 1.7, 0)

  ### check return value from CDS
  if(!is.na(trk)[1]) plot(trk) }

}
\keyword{manip}
