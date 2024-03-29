% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{sv_fore}
\alias{sv_fore}
\title{The forecasting model garch basd on \code{stochvol} package}
\usage{
sv_fore(x, train_h, PI_level)
}
\arguments{
\item{x}{The input time series.}

\item{train_h}{The amount of future time steps to forecast.}

\item{PI_level}{This function don't use it in the calculation.
This is used to make this function consistent with the parameter setting of the top four models.}
}
\value{
A list including forecasting mean and sd.
}
\description{
The forecasting model garch basd on \code{stochvol} package
}
