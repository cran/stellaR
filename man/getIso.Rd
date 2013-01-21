\name{getIso}
\alias{getIso}
\title{Import stellar isochrones data}
\description{
Import a stellar isochrone.
 }
\usage{
 getIso(age, z, y, ml, afe,
          baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/")
     }
\arguments{
\item{age}{the age (in Gyr) of the isochrone to
  get. Allowed values are 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5,
  12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0.}
\item{z}{the initial metallicity of the model to get. Allowed values
  are 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007,
    0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,
    0.008, 0.009, 0.01.}
\item{y}{the initial helium abundance of the model to get. Allowed
  values are 0.25, 0.27, 0.33, 0.38, 0.42.}
\item{ml}{the mixing-length of the model to get. Allowed values are
  1.7, 1.8, 1.9.}
\item{afe}{the alpha-enhanchment of the model to get. It can be one of
\code{afe = 0} for [alpha/Fe] = 0.0 or \code{afe = 1} for [alpha/Fe] =
0.3.}
\item{baseURL}{the URL of the base directory of the database from where get the models.}
}
\value{

  \code{getIso} returns an object of class \code{iso}.
  More specifically this is a list with elements:
  \item{age}{the age of the isochrone (in Gyr)}
  \item{z}{ the metallicity of the model}
  \item{y}{ the initial helium content of the model}
  \item{alpha.enh}{ the mixing-length value of the model}
  \item{data}{ the dataset of stellar isocrone, containing the following
    five variables:
   \describe{ 
  \item{logL}{logarithm of the luminosity (L in unit of solar luminosity)}
  \item{logTe}{logarithm of the effective temperature (Te in K)}
  \item{mass}{mass of the star (in unit of solar mass)}
  \item{radius}{radius of the star (in unit of solar radius)}
  \item{logg}{logarithm of the surface gravity (g in cm s^-2)}
}
}
  \code{\link{print}}, and \code{\link{plot}} methods are available for
  the class. 

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
\code{\link{getZahb}}, \code{\link{getHbgrid}}, \code{\link{getTrk}}, \code{\link{getTrkSet}}.
}

\references{
  
M. Dell'Omodarme, G. Valle, S. Degl'Innocenti, and P.G. Prada Moroni
(2012). \emph{The Pisa Stellar Evolution Data Base for low-mass
  stars}. Astronomy and Astrophysics, 540, A26.
}
\examples{
  \dontrun{
  iso <- getIso(12.0, 0.002, 0.25, 1.7, 0)

  ### get data from local directory /data
  iso <- getIso(12.0, 0.002, 0.25, 1.7, 0, baseURL="/data/")}
}
\keyword{manip}
