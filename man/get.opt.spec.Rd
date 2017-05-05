\name{get.opt.spec}
\alias{get.opt.spec}
\title{
  Option list
}
\description{
  Generate a list of options.
}
\usage{
get.opt.spec(p = 3, K = 30, alpha = 1, nlambda = 100, bs = "rb",
             lam.ratio.lims = c(1e-05, 1 - 1e-05), lam.ratio.len = 100,
             glmnet.lam.min = 0.001, bgamma = 1)
}
\arguments{
  \item{p}{
    Degree of polynomial used in the radial basis
}
  \item{K}{
    Number of knots
}
  \item{alpha}{
    The elastic net mixing parameter in \code{\link{glmnet}}
}
  \item{nlambda}{
    Number of lambda values tried in \code{\link{glmnet}}
}
  \item{bs}{
    Choice of smoothers; "rb" corresponds to radial basis
}
  \item{lam.ratio.lims}{
    Limits of the lambda ratios between the continuum and the emission lines.
}
  \item{lam.ratio.len}{
    Number of ratios
}
  \item{glmnet.lam.min}{
    lambda_min_ratio of \code{\link{glmnet}}
}
  \item{bgamma}{
    Tuning parameter in extended MDL
}
}
\value{
  An list of options.
}
\references{
Raymond K. W. Wong, Vinay L. Kashyap, Thomas C. M. Lee and David A. van Dyk (2016). \emph{Detecting Abrupt Changes in the Spectra of High-energy Astrophysical Sources}. The Annals of Applied Statistics, 10(2), 1107-1134. 
}

\seealso{
\code{\link{spec}}
}
