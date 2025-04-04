% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_iter.R
\name{makelags_iter}
\alias{makelags_iter}
\title{Make lag matrix for iterated prediction}
\usage{
makelags_iter(data, nfore, y, pop = NULL, E, tau, time, ...)
}
\arguments{
\item{data}{A data frame}

\item{nfore}{Number of steps to forecast}

\item{y}{Vector of column names}

\item{pop}{Column name for pop. Optional.}

\item{E}{Embedding dimension. Required.}

\item{tau}{Time delay. Required.}

\item{time}{Column name for time. Required.}

\item{...}{Additional arguments passed to makelags}
}
\value{
A data frame.
}
\description{
Creates a lag matrix that can be used for making iterated predictions.
Runs \code{makelags} with specified parameters and \code{forecast=TRUE},
then appends \code{nfore-1} empty rows for each pop with timesteps iterated by 1.
Requires use of a dataframe and provision of column names.
}
\examples{
RHfore2pop=makelags_iter(nfore=50, data=RickerHarvest, 
  y=c("CPUE_index","Catch"), 
  time="Time", tau=1, E=1, pop="Region")
}
\keyword{functions}
