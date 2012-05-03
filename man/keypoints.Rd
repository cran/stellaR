\name{keypoints}
\alias{keypoints}
\alias{keypoints.trk}
\alias{keypoints.trkset}
\alias{keypoints.iso}
\alias{keypoints.isoset}
\title{Extract relevant evolutionary points from stellar tracks or isochrones}
\description{
Extract some relevant evolutionary points from stellar tracks or
stellar isochrones.
 }
\usage{
\method{keypoints}{trk}(x, \dots)
\method{keypoints}{trkset}(x, \dots)
\method{keypoints}{iso}(x, \dots)
\method{keypoints}{isoset}(x, \dots)
     }
\arguments{
\item{x}{an object of classes \code{trk},
  \code{trkset}, \code{iso}, or \code{isoset}.}
\item{\dots}{further arguments, not implemented.}
}
\value{
  A data frame containing the observations relative to the selected
  evolutionary phases. The columns relative to mass (or age for
  isochrones), 
  metallicity, initial
  He value, mixing-length, alpha-enhancement, and phase identifier are
  inserted in the data frame. 
  }
\details{
\code{keypoints.trk} extracts from the object \code{x$data} the rows
relative to the following evolutionary stages:
\enumerate{
\item{ZAMS: }{ Zero-Age Main-Sequence, defined as the point for which the
  central H abundance drops below 99\% of its initial value. }
\item{TO: }{ Turn-Off, defined as the point for which the effective
  temperature reaches its maximum value. If multiple lines satisfy the
  constraint, the values of all the rows are averaged. }
\item{BTO: }{ Brighter Turn-Off, defined as the point for which the
  effective temperature drops below the one of the TO minus 100 K. The
  points could not exist for low masses. See Chaboyer et
  al. (1996) for an analysis of the advantages of this evolutionary
point with respect to TO.}
\item{exHc: }{ Central H exhaustion, defined as the point for which the central
  H abundance is zero. For low masses the point could coincide with
  TO. This is the last point of the tracks with mass lower or equal to
  0.50 Msun.}
\item{Heflash: }{ Helium flash, the last point of the track for masses
  higher than 0.50 Msun.}
}

\code{keypoints.trkset} performs the selection described above for all the
set of objects.

The returned object contains the variable \code{id}, which labels the
evolutionary phases with the following coding: 1 = ZAMS, 2 = TO, 3 =
BTO, 4 = exHc, 5 = Heflash.

\code{keypoints.iso} extracts from the object \code{x$data} the rows
relative to the following evolutionary stages:
\enumerate{
\item{TO: }{ Turn-Off, defined as the point for which the effective
  temperature reaches its maximum value. If multiple lines satisfy the
  constraint, the values of all the rows are averaged. }
\item{BTO: }{ Brighter Turn-Off, defined as the point for which the
  effective temperature drops below the one of the TO minus 100 K. }
}

\code{keypoints.isoset} performs the selection described above for all the
set of objects.

The returned object contains the variable \code{id}, which labels the
evolutionary phases with the following coding: 1  = TO, 2 = BTO.

}
\references{B. Chaboyer, P. Demarque, P.J. Kernan, L.M. Krauss, and
  A. Sarajedini (1996). \emph{An accurate relative age estimator for
    globular clusters}. MNRAS, 283, 683-689.
}
\examples{

  trk <- getTrk(0.9, 0.002, 0.25, 1.7, 0)
  ### check return value from CDS
  if(!is.na(trk)[1]) keypoints(trk)

  is <- getIso(11, 0.002, 0.25, 1.7, 0)
  ### check return value from CDS
  if(!is.na(is)[1]) keypoints(is)
}
\keyword{manip}
