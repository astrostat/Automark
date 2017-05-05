\name{Automark-package}
\alias{Automark-package}
\alias{Automark}
\docType{package}
\title{
Automark
}
\description{
This package implements a method for modeling photon counts collected form observation of variable-intensity astronomical sources. It aims to mark the abrupt changes in the corresponding wavelength distribution of the emission automatically. In the underlying methodology, change points are embedded into a marked Poisson process, where photon wavelengths are regarded as marks and both the Poisson intensity parameter and the distribution of the marks are allowed to change. The details is given in Wong et. al. (2016).
}
\author{
Raymond K. W. Wong, Vinay L. Kashyap, Thomas C. M. Lee and David A. van Dyk

Maintainer: Raymond K. W. Wong <raymondkww.dev@gmail.com>
}
\references{
Raymond K. W. Wong, Vinay L. Kashyap, Thomas C. M. Lee and David A. van Dyk (2016) \emph{Detecting Abrupt Changes in the Spectra of High-energy Astrophysical Sources}. The Annals of Applied Statistics, 10(2), 1107-1134.
}
\keyword{ package }
\examples{
library(Automark)           # load library
data(dat)                   # load example data
y1 <- apply(dat$y, 1, sum)  # aggregate data over time
\dontrun{
out <- spec(dat$x2, y1, dat$A, dat$delta.t, dat$delta.w, reps=ncol(dat$y))
plotspec(dat$x2, out)
}
}
