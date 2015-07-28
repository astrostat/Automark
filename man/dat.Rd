\name{dat}
\alias{dat}
\docType{data}
\title{
  Example dataset
}
\description{
  For demostration purpose, we provide an artificial dataset derived from a Chandra gratings observation of FK Com (ObsID 12298) with an extra background contamination.
}
\usage{data("dat")}
\references{
Raymond K. W. Wong, Vinay L. Kashyap, Thomas C. M. Lee and David A. van Dyk (2015). \emph{Detecting Abrupt Changes in the Spectra of High-energy Astrophysical Sources}. Submitted.
}
\examples{
library(Automark)           # load library
data(dat)                   # load example data
y1 <- apply(dat$y, 1, sum)  # aggregate data over time
\dontrun{
out <- spec(dat$x2, y1, dat$A, dat$delta.t, dat$delta.w, reps=ncol(dat$y))
plotspec(dat$x2, out)
}
}
\keyword{datasets}
