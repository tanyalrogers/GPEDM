\name{GPEDM-package}
\alias{GPEDM-package}
\alias{GPEDM}
\docType{package}
\title{
  Gaussian Process regression for Empirical Dynamic Modeling
}
\description{
  Functions for fiting hierarchical, separable length scale GP models with
  automatic relevance determination (ARD) for use in Empirical Dynamic Modeling (EDM) 
  and other applications.
}
\details{
  The main function is \code{\link{fitGP}} which is used to train the model and can also 
  produce predictions if desired. Use \code{\link{predict.GP}} to generate other or additional
  predictions from a fitted model, \code{\link{plot.GPpred}} to plot observed and predicted 
  values and \code{\link{getconditionals}} to obtain condictional reponses. Also available
  for use is the function \code{\link{makelags}} which can be used to create delay vectors.
  See the (yet to be created) vignette for more detailed instructions.
}
\author{
Original MATLAB code by Stephan Munch, adapted to R by Tanya Rogers.

Maintainer: Tanya Rogers <tanya.rogers@noaa.gov>
}
\references{
Munch, S. B., Poynor, V., and Arriaza, J. L. 2017. Circumventing structural uncertainty:
a Bayesian perspective on nonlinear forecasting for ecology. Ecological Complexity, 32:134.
}
\keyword{ package }
\seealso{
  \code{\link{fitGP}}, \code{\link{plot.GPpred}}, \code{\link{getconditionals}}, 
  \code{\link{makelags}}
}
\examples{
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.   
  }
}
