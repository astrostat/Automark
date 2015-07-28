\name{spec}
\alias{spec}
\title{
  Fit a spectrum with no time component / homogeneous across time
}
\description{
  Fit a spectrum with no time component / homogeneous across time. The corresponding spectrum is modeled nonparametrically using a mixture of a smooth radial basis expansion and a number
  of local deviations from the smooth term representing spectral emission lines.
  The corresponding tuning parameter is determined via the minimum description length principle.
}
\usage{
spec(x, y, A, delta.t, delta.w, reps = 1, display = TRUE,
opt = get.opt.spec(), assign.emiss = NA, emdl = TRUE,
simple = TRUE, v = 1) 
}
\arguments{
  \item{x}{
    A vector containing wavelengths, corresponding to the vector of observed
    counts y
  }
  \item{y}{
    A vector containing the observed counts or aggregated observed counts over homogeneous time segments
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
  \item{reps}{
    Number of time bins that you aggregate the counts over
  }
  \item{display}{
    If TRUE, progress bar is printed during fitting
  }
  \item{opt}{
    An option list returned by \code{\link{get.opt.spec}}
  }
  \item{assign.emiss}{
    Indices of wavelength bins that are allowed to have emission line.
    If NA, assign_emiss will be set to 1:length(y).
  }
  \item{emdl}{
    If TRUE, extended MDL, rather than MDL, is used.
  }
  \item{simple}{
    If TRUE, a simplified output is returned. That is, results of
            the grid search (beta_cube) will be not be output.
  }
  \item{v}{
    An ad-hoc parameter for tuning the strength of the penalty. That
            means, (v * penalty) will be the ultimate penalty used in the MDL /
            extended MDL. More penalty leads to a smoother spectrum and/or fewer
            emission lines. It is set to 1.0 by default.
  }
}
\value{
  \item{best.beta}{
    Regression coefficients of the output spectrum
  }
  \item{fitted.bright}{
    Linear predictor of the output spectrum; exp(fitted_bright) is the fitted brightness function.
  }
  \item{fitted.bright.ne}{
    Linear predictor without emission lines
  }
  \item{df}{
    df[1]: Degree of freedoms of the continuum (does not count the intercept)\\
    df[2]: Degree of freedoms of the emission lines
  }
  \item{out.ind}{
    Indices of the time points with emission lines
  }
  \item{mdl}{
    All MDLs in the grid search of the tuning parameter
  }
  \item{mdl.ind}{
    Index of the smallest MDL in the grid search of the tuning parameter
  }
  \item{mdl.pen}{
    All penalties in MDLs in the grid search of the tuning parameter
  }
  \item{beta.cube}{
    All fitted regression coefficients of all fits in the grid search of the tuning parameter
  }
}
\references{
Raymond K. W. Wong, Vinay L. Kashyap, Thomas C. M. Lee and David A. van Dyk (2015). \emph{Detecting Abrupt Changes in the Spectra of High-energy Astrophysical Sources}. Submitted.
}
\seealso{
  \code{\link{plotspec}}
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
