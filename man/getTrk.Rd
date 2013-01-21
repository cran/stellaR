\name{getTrk}
\alias{getTrk}
\alias{getHb}
\title{Import stellar track data}
\description{
Import a stellar evolutionary track from pre-main sequence to He flash
and from ZAHB to thermal pulses.
 }
\usage{
 getTrk(m, z, y, ml, afe,
          baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/")
 getHb(m, z, y, ml, afe,
         baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/")
     }
\arguments{
\item{m}{the mass, in unit of solar mass, of the model to
  get. Allowed values are 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60,
  0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10.}
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
\item{baseURL}{the URL of the base directory of the database from where
  get the models. See \code{Details} for further information.}
}
\value{

  \code{getTrk} returns an object of class \code{trk}; \code{getHb}
  returns an object of class \code{hb}.
  
  More specifically these objects are lists with common elements:
  \item{mass}{the mass of the model (in unit of solar mass)}
  \item{z}{ the metallicity of the model}
  \item{y}{ the initial helium content of the model}
  \item{ml}{ the mixing-length value of the model}
  \item{alpha.enh}{the alpha-enhanchment of the model. It can be one of
\code{alpha.enh = 0} for [alpha/Fe] = 0.0 or \code{alpha.enh = 1} for
  [alpha/Fe] = 0.3.  }
  \item{data}{ the dataset of stellar track, containing the following
    variables:
    \describe{ 
      \item{mod}{model number}
      \item{time}{logarithm of the evolutionary time (in yr)}
      \item{logL}{logarithm of the luminosity (L in unit of solar luminosity)}
      \item{logTe}{logarithm of the effective temperature (Te in K)}
      \item{mass}{mass of the star (in unit of solar mass)}
      \item{Hc}{central H/He abundance}
      \item{logTc}{logarithm of the central temperature (Tc in K)}
      \item{logRHOc}{logarithm of the central density (RHOc in g cm^-3)}
      \item{MHEc}{mass of the He core (in unit of solar mass)}
      \item{Lpp}{luminosity of pp chain (in unit of total luminosity L)}
      \item{LCNO}{luminosity of CNO cycle (in unit of total luminosity L)}
      \item{L3a}{luminosity of triple-alpha (in unit of total luminosity
  L)}
      \item{Lg}{gravitational luminosity (in unit of total luminosity L)}
      \item{radius}{radius of the star (in unit of solar radius)}
      \item{logg}{logarithm of the surface gravity (g in cm s^-2)} 
    }
    
    Only \code{time}, \code{logL}, and \code{logTe} are mandatory, while
  all other variables are optional.
  }

  \code{hb} class contains also the element:
  \item{massRGB}{the mass of the RGB progenitor}

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
\code{\link{getZahb}}, \code{\link{getHbgrid}}, \code{\link{getIso}}, \code{\link{getTrkSet}}.
}

\references{
  
M. Dell'Omodarme, G. Valle, S. Degl'Innocenti, and P.G. Prada Moroni
(2012). \emph{The Pisa Stellar Evolution Data Base for low-mass
  stars}. Astronomy and Astrophysics, 540, A26.
}
\examples{

  \dontrun{
  trk <- getTrk(0.9, 0.002, 0.25, 1.7, 0)

  ### get data from local directory /data
  trk <- getTrk(0.9, 0.002, 0.25, 1.7, 0, baseURL="/data/")

  ###  multi-panel plot of the various quantities with time
  track <- getTrk(0.80, 0.001, 0.25, 1.90, 0)
  if(!is.na(track)[1]) {
  trkdata <- within(stack(track$data), time <- rep(track$data$time,
      length.out=length(values)) )
  require(lattice)
  xyplot( values ~ time | ind, data=trkdata, type="l",
      scales=list(y=list(relation="free")))
  } }
}
\keyword{manip}
