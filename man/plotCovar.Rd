% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotCovar}
\alias{plotCovar}
\title{Plot covariate coefficients}
\usage{
plotCovar(
  object,
  covar = NA,
  tlab = NA,
  return_data = FALSE,
  myTheme = ggplot2::theme_bw(),
  pvalue = "shape"
)
}
\arguments{
\item{object}{An ALASCA object}

\item{covar}{Which covariable(s) to plot (default: \code{NA} which prints all)}

\item{tlab}{Alternative names for the covariables}

\item{return_data}{Set to \code{TRUE} to return data instead of plot}

\item{myTheme}{A ggplot2 theme to use, defaults to \code{ggplot2::theme_bw()}}
}
\value{
A ggplot2 objects\.
}
\description{
This function returns a plot of the regression coefficients for covariates that is not included in the ASCA model itself
}
