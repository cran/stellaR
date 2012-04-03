\name{masshb}
\docType{data}
\alias{masshb}
\alias{masshbgrid}
\title{Mass of stars in HB}
\description{
  These data sets give the stellar mass in HB for the progenitors in RGB.  
}
\usage{
masshb
masshbgrid
}
\format{The two data frames contain 6 variables:
  \describe{
    \item{mass}{the RGB mass of the model.}
    \item{z}{the initial metallicity of the model.}
    \item{y}{the initial helium abundance of the model.}
    \item{ml}{the mixing-length of the model.}
    \item{mix}{the alpha-enhanchment of the model. It is coded
      \code{mix=AS09a0} for [alpha/Fe] = 0.0 and \code{mix=AS09a3} for
      [alpha/Fe] = 0.3.} 
    \item{masshb}{the HB mass of the model. The dataset
      \code{masshbgrid} contains multiple values for \code{masshb} for
      each single RGB progenitor mass.}
    }
}
\source{M. Dell'Omodarme, G. Valle, S. Degl'Innocenti, and P.G. Prada Moroni
(2012). \emph{The Pisa Stellar Evolution Data Base for low-mass
  stars}. Astronomy and Astrophysics, 540, A26.}

\keyword{datasets}
