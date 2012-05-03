\name{getTrkSet}
\alias{getTrkSet}
\alias{getIsoSet}
\title{Import a set of data}
\description{
Import a set of stellar evolutionary tracks from pre-main sequence to He
flash, or a set of isochrones.
 }
\usage{
 getTrkSet(m, z, y, ml, afe,
             baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/")
 getIsoSet(age, z, y, ml, afe,
             baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/")
     }
\arguments{
\item{m}{a vector of masses, in unit of solar mass, of the models to
  get. Allowed values are 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60,
  0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10.}
\item{age}{a vector of ages, in Gyr, of the models to
  get. Allowed values are 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5,
  12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0.}
\item{z}{a vector of the initial metallicity of the models to
  get. Allowed values 
  are 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007,
    0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,
    0.008, 0.009, 0.01.}
\item{y}{a vector of the initial helium abundance of the models to
  get. Allowed 
  values are 0.25, 0.27, 0.33, 0.38, 0.42.}
\item{ml}{a vector of the mixing-length of the models to get. Allowed
  values are 1.7, 1.8, 1.9.}
\item{afe}{the alpha-enhanchment of the models to get. It can be one of
\code{afe = 0} for [alpha/Fe] = 0.0 or \code{afe = 1} for [alpha/Fe] =
0.3.}
\item{baseURL}{the URL of the base directory of the database from where get the models.}
}
\value{

  \code{getTrkSet} returns an object of class \code{trkset}, i.e. a
  list of objects of class \code{trk}.
  \code{getIsoSet} returns an object of class \code{isoset}, i.e. a
  list of objects of class \code{iso}.

  \code{\link{print}}, and \code{\link{plot}} methods are available for
  the two classes.

  If \code{baseURL} points to an invalid path or CDS ftp site is
  unavailable
  the function returns \code{NA} and issues a warning.
}
\details{

  As a sanity check, the function performs a consistency test on the parameters
  \code{z}, \code{y}, \code{ml}, \code{afe} by a call to
  \code{\link{testComposition}}.

  By default the function requires an Internet access for data
  gathering. Tracks and isochrones are downloaded from the on-line
  repository CDS
  (\url{http://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/540/A26}).
  It is also possible to gather data from a local download of the
  catalog, specifying as \code{baseURL} the root directory of the
  download. In this case the structure of the catalog and the names of
  files and directory must be preserved.
}
\seealso{
\code{\link{getZahb}}, \code{\link{getHbgrid}},
\code{\link{getHb}}, \code{\link{getIso}}.
}

\references{
  
M. Dell'Omodarme, G. Valle, S. Degl'Innocenti, and P.G. Prada Moroni
(2012). \emph{The Pisa Stellar Evolution Data Base for low-mass
  stars}. Astronomy and Astrophysics, 540, A26.
}
\examples{
  ### get two masses
  trkset <- getTrkSet(c(0.9, 1.0), 0.002, 0.25, 1.7, 0)

  ### get two masses at three metallicity, for a total of 6 objects
  trkset <- getTrkSet(c(0.9, 1.0), c(0.002, 0.006, 0.01), 0.25, 1.7, 0)

  ### get data from local directory /data
  \dontrun{set <- getTrkSet(0.9, 0.002, c(0.25, 0.33), 1.7, 0, baseURL="/data/")}

   ### get four isochrones
  isoset <- getIsoSet(c(9, 10, 12), 0.002, 0.25, 1.7, 0)
}
\keyword{manip}
