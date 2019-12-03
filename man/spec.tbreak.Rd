\name{spec.tbreak}
\alias{spec.tbreak}
\title{
  Detect change points in spectrum across time
}
\description{
  Detect change points in spectrum across time. Between the change points, the spectra are
  fitted by \code{\link{spec}}. The number of change points is determined via the minimum description length principle.
}
\usage{
spec.tbreak(x1, x2, y, A, delta.t, delta.w, display = TRUE,
opt = get.opt.spec(), cpus = 4, assign.emiss = NA,
emdl = TRUE, v1 = 1, v2 = 1, simple = TRUE,
max.B = 4, min.span = 5)
}
\arguments{
  \item{x1}{
    A vector containing time points (mid-points of time bins)
  }
  \item{x2}{
    A vector containing wavelengths (mid-points of wavelength bins)
  }
  \item{y}{
    A matrix containing the observed counts. Its dimension is c(number of wavelengths, number of time points)
  }
  \item{A}{
    A vector containing the effective areas
  }
  \item{delta.t}{
    Bin width of time
  }
  \item{delta.w}{
    Bin width of wavelength
  }
  \item{display}{
    If TRUE, progress will be shown.
  }
  \item{opt}{
    An option list returned by \code{\link{get.opt.spec}}. This is used in the homogeneous spectrum fitting between change points. If assign.emiss=NULL, nlambda will be automatically set to 10000 if it is set smaller than 10000. If this automatic change is not desired, one can set lam.ratio.len=-1. (lam.ratio.len is a useless parameter if assign.emiss=NULL)
  }
  \item{cpus}{
    Number of cores to use for parallel execution.
  }
  \item{assign.emiss}{
    Indices of wavelength bins that are allowed to have emission line. If NA, assign_emiss will be set to 1:length(y).
  }
  \item{emdl}{
    If TRUE, extended MDL, rather than MDL, is used.
  }
  \item{v1}{
    An ad-hoc parameter for tuning the strength of the penalty in fitting each time-homogenous spectrum between change points. That means, (v * penalty) will be the ultimate penalty used in the MDL / extended MDL. More penalty leads to a smoother spectrum and/or fewer emission lines. It is set to 1.0 by default.
  }
  \item{v2}{
    Similar to v1, but for tuning the strength of the penalty for change points across time. Larger v2 leads to fewer change points across time.
  }
  \item{simple}{
    If TRUE, a simplified output is returned. That is, results of
            the grid search (beta_cube) will be not be output.
  }
  \item{max.B}{
    The maximum number of time segments allowed across time.
  }
  \item{min.span}{
    The minimum number of time bins allowed between breaks.
  }
}
\value{
  \item{B}{
    Number of time segments in the output fit
  }
  \item{best.fit}{
    List of time-homogeneous fits (in the output fit) corresponding to different time segments.
    E.g., best.fit[[1]] corresponds to the first time segment.
  }
  \item{breaks}{
    Indices of breakpoints. E.g., c(1,10) represents the first time segment spans time bins 1,2,...,9 and the second one spans the rest.
  }
  \item{content}{
    Information of the fitting process (debug use)
  }
}
\references{
Raymond K. W. Wong, Vinay L. Kashyap, Thomas C. M. Lee and David A. van Dyk (2016). \emph{Detecting Abrupt Changes in the Spectra of High-energy Astrophysical Sources}. The Annals of Applied Statistics, 10(2), 1107-1134.
}
\seealso{
  \code{\link{spec}}, \code{\link{plotspec}}
}
\examples{
library(Automark)       # load library
data(dat)               # load example data
\dontrun{
out <- spec.tbreak(dat$x1, dat$x2, dat$y, dat$A, dat$delta.t, dat$delta.w, cpus=5)
plotspec(dat$x2, out$best.fit[[1]])
plotspec(dat$x2, out$best.fit[[2]], np=F, col="red")
}
}
