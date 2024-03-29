% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ALASCA.R
\name{ALASCA}
\alias{ALASCA}
\title{Create an ALASCA model}
\usage{
ALASCA(df, formula, effects = NULL, ...)
}
\arguments{
\item{df}{Data frame to be analyzed}

\item{formula}{Regression formula}

\item{effects}{Vector with effect terms. If \code{NULL}, ALASCA will guess (default)}

\item{scale_function}{Either a custom function or string to define scaling function: \code{sdall}, \code{sdref}, \code{sdt1}, \code{sdreft1}}

\item{separate_effects}{Boolean. Separate effects?}

\item{equal_baseline}{Set to \code{TRUE} to remove interaction between effects}

\item{validate}{Boolean. If \code{TRUE}, give estimates for robustness}

\item{n_validation_runs}{number of validation runs}

\item{validation_method}{Choose between \code{bootstrap} (default) and "jack-knife"}

\item{stratification_column}{The column to stratify participants by}

\item{save}{Save models and plots automatically (default: \code{FALSE})}

\item{filename}{File name to save model and plots (when \code{save = TRUE})}

\item{use_Rfast}{Boolean. Defaults to \code{TRUE}}

\item{p_adjust_method}{Method for correcting p values for multiple testing, see p.adjust.methods}

\item{participant_column}{String. Name of the column containing participant identification}

\item{n_validation_folds}{Partitions when validating}
}
\value{
An ALASCA object
}
\description{
\code{ALASCA} initializes an ALASCA model and returns an ALASCA object.
}
\details{
This function builds your ALASCA model. It needs, at least, a data frame and a formula. The effect matrices can be specified with \code{effects}, e.g., \code{c("time", "time+group+time:group", "group+time:group")}.
}
