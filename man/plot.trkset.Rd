\name{plot.trkset}
\alias{plot.trkset}
\alias{plot.hbset}
\alias{plot.isoset}
\title{Plot a set stellar objects}
\description{
Plot a set of evolutionary tracks from objects of
classes \code{trkset}, and
\code{hbset}, or a set of isochrones from object of class \code{isoset}.
 }
\usage{
\method{plot}{trkset}(x, add=FALSE, col=1, lty=1, xlim=NULL, ylim=NULL, \dots)
\method{plot}{hbset}(x, add=FALSE, col=1, lty=1, xlim=NULL, ylim=NULL, \dots)
\method{plot}{isoset}(x, add=FALSE, col=1, lty=1, xlim=NULL, ylim=NULL, \dots)
     }
\arguments{
\item{x}{an object of class \code{trkset},
  \code{hbset}}, or \code{isoset}.
\item{add}{a logical value indicating whether the plot should be added to
  an existing device.}
\item{col}{the colors for lines.  Multiple colors can be
          specified so that each track can be given its own color.  If
          there are fewer colors than tracks they are recycled in the
          standard fashion. }
\item{lty}{the line types.  Multiple types can be
          specified so that each track can be given its own type.  If
          there are fewer type than tracks they are recycled in the
          standard fashion.}
\item{xlim}{the x limits (x1, x2) of the plot.  
  The default value, \code{NULL}, indicates that the range of the
  finite values to be plotted should be used.
}
\item{ylim}{the y limits of the plot.}
\item{\dots}{further arguments passed to \code{\link{plotAstro}}.}
}

\details{

For each object \code{X} contained in the input of class \code{trkset},
\code{hbset}, or 
\code{isoset} the contents 
of the vectors \code{X$data$logTe} and \code{X$data$logL} are 
plotted by a call to \code{\link{plotAstro}}. Sensible values for the
axes ranges, major ticks and labels are chosen.
}
\seealso{
  \code{\link{plotAstro}}.
}

\examples{

  trkset <- getTrkSet(c(0.7, 0.9, 1.1), 0.002, 0.25, 1.7, 0)
  plot(trkset)

}
\keyword{hplot}
