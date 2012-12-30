\name{compositions}
\docType{data}
\alias{compositions}
\title{Compostions of computed models}
\description{
  The data sets gives the compositions of the computed models.  
}
\usage{
compositions
}
\format{A list containing 4 variables:
  \describe{
    \item{z}{the initial metallicities of the models.}
    \item{y}{a data frame with the initial helium abundances of the models.}
    \item{ml}{the mixing-lengths of the models.}
    \item{afe}{the alpha-enhanchments of the models. It is coded
      \code{0} for [alpha/Fe] = 0.0 and \code{1} for
      [alpha/Fe] = 0.3.} 
    }
}
\source{M. Dell'Omodarme, G. Valle, S. Degl'Innocenti, and P.G. Prada Moroni
(2012). \emph{The Pisa Stellar Evolution Data Base for low-mass
  stars}. Astronomy and Astrophysics, 540, A26.}

\keyword{datasets}
