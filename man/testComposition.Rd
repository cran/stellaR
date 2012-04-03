\name{testComposition}
\alias{testComposition}
\title{Check the existence of a record in the database}
\description{
Check that the given chemical composition (Z, Y, [alpha/Fe]) and the
mixing-length value exist in the stellar database.  
 }
 \usage{
testComposition(Z, Y, ML, AFE)
}
\arguments{
\item{Z}{the initial metallicity of the model.}
\item{Y}{the initial helium abundance of the model.}
\item{ML}{the mixing-length of the model.}
\item{AFE}{the alpha-enhanchment of the model. Allowed values are \code{AFE =
0} ([alpha/Fe] = 0.0) or \code{AFE = 1} ([alpha/Fe] = 0.3).}
}
\value{
  The function returns \code{TRUE} if the combination of the inputs
  exists in the database. 
}
\seealso{
\code{\link{showComposition}}.
}

\references{
  
M. Dell'Omodarme, G. Valle, S. Degl'Innocenti, and P.G. Prada Moroni
(2012). \emph{The Pisa Stellar Evolution Data Base for low-mass
  stars}. Astronomy and Astrophysics, 540, A26.
}
\examples{

  testComposition(0.002, 0.25, 1.7, 0)
}
\keyword{manip}
