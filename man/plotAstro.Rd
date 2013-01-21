\name{plotAstro}
\alias{plotAstro}
\title{Plot function for stellar objects}
\description{
Draw a scatter plot in the active graphics window.
 }
 \usage{

plotAstro(x, y, type="l", xlab="", ylab="", xi=(max(x)-min(x))/5,
          yi=(max(y)-min(y))/5, xmt=3, ymt=3, revX=FALSE, revY=FALSE,
          xlim=NULL, ylim=NULL, cex=1.0, cex.axis=1.3, cex.lab=1.5,
          add=FALSE, \ldots)
     }
\arguments{
\item{x, y}{the \code{x} and \code{y} arguments provide the x and y
     coordinates for 
  the plot.  }
\item{type}{	
  1-character string giving the type of plot desired.  The 
  following values are possible: \code{"p"} 
  for points, \code{"l"} for lines (default value), \code{"b"} for both
  points and lines, 
  \code{"c"} for empty points joined by lines, \code{"o"} for overplotted
  points and lines, \code{"s"} and \code{"S"} for stair steps and \code{"h"}
  for histogram-like vertical lines.  Finally, \code{"n"} does not
  produce any points or lines.}
\item{xlab}{a label for the x axis, default to empty string.}
\item{ylab}{a label for the y axis, default to empty string.}
\item{xi}{major ticks spacing for x axis.}
\item{yi}{major ticks spacing for y axis.}
\item{xmt}{number of minor ticks between two major ticks for x axis.}
\item{ymt}{number of minor ticks between two major ticks for y axis.}
\item{revX}{a logical value indicating whether the x axis should be
  reversed.}
\item{revY}{a logical value indicating whether the y axis should be
  reversed.}
\item{xlim}{the x limits (x1, x2) of the plot.  
  The default value, \code{NULL}, indicates that the range of the
  finite values to be plotted should be used.
}
\item{ylim}{the y limits of the plot.}
\item{cex}{the amount by which plotting text
  and symbols should be magnified relative to the default.
  }
\item{cex.axis}{the magnification to be used for axis annotation
  relative to the current setting of \code{cex}.}
\item{cex.lab}{the magnification to be used for x and y labels
  relative to the current setting of \code{cex}.}
\item{add}{a logical value indicating whether the plot should be added to
  an existing device.}
\item{\dots}{further arguments passed to low level plot functions.}
}

\seealso{
  \code{\link{plot.trk}}, \code{\link{plot.trkset}}.
}
\examples{

  \dontrun{
  trk <- getTrk(0.9, 0.002, 0.25, 1.7, 0)

  ### check return value from CDS
  if(!is.na(trk)[1]) plot(trk) }

}
\keyword{hplot}
