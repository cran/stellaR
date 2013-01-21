\name{interpTrk}
\alias{interpTrk}
\title{Interpolate stellar tracks data}
\description{
Construct a set of interpolated tracks from the objects in the database.
 }
\usage{
 interpTrk(z, y, ml, afe, vmass=seq(0.30,1.10, by=0.05),
           baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/")
     }
\arguments{
\item{z}{the initial metallicity of the models to generate. Allowed values
  are between \code{0.0001} and \code{0.01}.}
\item{y}{the initial helium abundance of the models to generate. Allowed
  values are between \code{0.25} and \code{0.42}.}
\item{ml}{the mixing-length of the model to generate. Allowed values are
  between \code{1.7} and \code{1.9}.}
\item{afe}{the alpha-enhanchment of the model to get. It can be one of
\code{afe = 0} for [alpha/Fe] = 0.0 or \code{afe = 1} for [alpha/Fe] =
0.3.}
\item{vmass}{vector of masses of the stellar model to generate. Default
  to the whole set in the database. Allowed values are
  0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60,
  0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10. }
\item{baseURL}{the URL of the base directory of the database from where get the models.}
}
\value{

  The function returns an object of class \code{trkset} containing the
  interpolated tracks.
}
\details{
  A 3D interpolation is performed on the database to construct the
  required set of tracks. Interpolation on helium abundance and on
  mixing-length occur linearly on the variables, while the metallicity
  is log-transformed before the interpolation procedure. 

}
\examples{

  ### slow!
  \dontrun{ 
  iptrk <- interpTrk(0.002, 0.25, 1.74, 0)

  ### get data from local directory /data
  iptrk <- interpTrk(0.002, 0.25, 1.74, 0, baseURL="/data/")}
}
\keyword{math}
