% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{PLAC_TI}
\alias{PLAC_TI}
\title{C++ Function for Solving the PLAC Estimator.
(with time-invariant convariates only)}
\usage{
PLAC_TI(Z, X, W, Ind1, Ind2, Dn, b, h, K = 100L)
}
\arguments{
\item{Z}{matrix for all the covariates history.}

\item{X}{the response matrix (As, Xs, Ds).}

\item{W}{the ordered observed event times.}

\item{Ind1}{risk-set indicators.}

\item{Ind2}{truncation pair indicators.}

\item{Dn}{number of ties at each observed event time.}

\item{b}{initial values of the regression coefficients.}

\item{h}{initial values of the baseline hazard function.}

\item{K}{maximal iteration number, the default is \code{K = 100}.}
}
\value{
list of model fitting results for both conditional approach and
the PLAC estimator.
}
\description{
C++ Function for Solving the PLAC Estimator.
(with time-invariant convariates only)
}
