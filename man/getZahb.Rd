\name{getZahb}
\alias{getZahb}
\title{Import stellar track data for ZAHB }
\description{
Import a stellar ZAHB model.
 }
\usage{
 getZahb(z, y, ml, afe,
           baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/")
     }
\arguments{
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

  \code{getZahb} returns an object of class \code{zahb}
  More specifically it is a list with elements:
  \item{z}{ the metallicity of the model}
  \item{y}{ the initial helium content of the model}
  \item{ml}{ the mixing-length value of the model}
  \item{alpha.enh}{the alpha-enhanchment of the model. It can be one of
    \code{alpha.enh = 0} for [alpha/Fe] = 0.0 or \code{alpha.enh = 1} for
    [alpha/Fe] = 0.3.  }
  \item{data}{ the dataset of stellar track, containing the following
    variables:
    \describe{ 
      \item{mass}{mass of the star (in unit of solar mass)}
      \item{logTe}{logarithm of the effective temperature (Te in K)}
      \item{logL}{logarithm of the luminosity (L in unit of solar luminosity)}
    }
  }

  \code{\link{print}}, and \code{\link{plot}} methods are available for
  the two classes. 
}
\details{

  As a sanity check, the function performs a consistency test on the parameters
  \code{z}, \code{y}, \code{ml}, \code{afe} by a call to
  \code{\link{testComposition}}.

  By default the function requires an Internet access for data
  gathering. Tracks and isochrones are downloaded from the on-line
  repository CDS
  (\url{http://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/540/A26})
  It is also possible to gather data from a local download of the
  catalog, specifying as \code{baseURL} the root directory of the
  download. In this case the structure of the catalog and the names of
  files and directory must be preserved.
}
\seealso{
\code{\link{getTrk}}, \code{\link{getHbgrid}}, \code{\link{getIso}}, \code{\link{getTrkSet}}.
}

\references{
  
M. Dell'Omodarme, G. Valle, S. Degl'Innocenti, and P.G. Prada Moroni
(2012). \emph{The Pisa Stellar Evolution Data Base for low-mass
  stars}. Astronomy and Astrophysics, 540, A26.
}
\examples{

  zahb <- getZahb(0.002, 0.25, 1.7, 0)

  ### get data from local directory /data
  \dontrun{zahb <- getZahb(0.002, 0.25, 1.7, 0, baseURL="/data/")}
}
\keyword{manip}
