\name{makeIso}
\alias{makeIso}
\title{Construct stellar isochrones from tracks}
\description{
The function computes by interpolation stellar isochrones of given
ages from evolutionary tracks. 
 }
\usage{
 makeIso(age, z=NULL, y=NULL, ml=NULL, afe=NULL, log=FALSE, linear=TRUE,
         tr=NULL, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/")
     }
\arguments{
\item{age}{the ages (in Gyr) of the isochrones to
  construct. A vector of ages is allowed.}
\item{z}{the initial metallicity of the isochrones to build.}
\item{y}{the initial helium abundance of the isochrones to build.}
\item{ml}{the mixing-length of the isochrones to build.}
\item{afe}{the alpha-enhanchment of the isochrone to build. It can be one of
\code{afe = 0} for [alpha/Fe] = 0.0 or \code{afe = 1} for [alpha/Fe] =
0.3.}
\item{log}{a logical value indicating whether the interpolation of the
  tracks should be performed on the logarithm of the evolutionary time.}
\item{linear}{a logical value indicating whether linear interpolation
  should be performed. Only the option \code{linear=TRUE} is currently
  available. } 
\item{tr}{a set of track to be used for isocrones construction. See
  \code{Details} for further information.}
\item{baseURL}{the URL of the base directory of the database from where get the models.}
}
\value{

  The function returns an object of class \code{isoset}.
}
\details{

  Isochrones are obtained by mean of a standard interpolation procedure on
  the set of tracks. Let \eqn{S(m)} be the set of tracks, parametrized
  by the value of the mass \eqn{m}.
  Let \eqn{t_i(m)} be the evolutionary
  time for the ith point on the track of mass \eqn{m}.
  Let be \eqn{k} the point on the track of lower mass of \eqn{S(m)} for
  which \eqn{t_k(m)} is greater of the time required for the isochrone.

  For each point \eqn{j >= k} on \eqn{S(m)}, an interpolation of mass,
  logarithm of the effective temperature and logarithm of the luminosity
  is performed among tracks. These points define the required isochrone.
  
  If a set of tracks is supplied by mean of the argument \code{tr}, the
  function interpolates among these tracks. The values of
  \code{z}, \code{y}, \code{ml}, \code{afe} are recovered from the
  supplied objects and a test of consistency is performed to assure that
  the tracks are homogeneous in these parameters.   

}
\references{
  
M. Dell'Omodarme, G. Valle, S. Degl'Innocenti, and P.G. Prada Moroni
(2012). \emph{The Pisa Stellar Evolution Data Base for low-mass
  stars}. Astronomy and Astrophysics, 540, A26.
}
\examples{

  isoset <- makeIso(c(11.2, 12.4), 0.002, 0.25, 1.7, 0)

  ### get data from local directory /data
  \dontrun{isoset <- makeIso(c(11.2, 12.4), 0.002, 0.25, 1.7, 0, baseURL="/data/")}
}
\keyword{math}
