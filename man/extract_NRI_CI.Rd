% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rap.R
\name{extract_NRI_CI}
\alias{extract_NRI_CI}
\title{Extract NRI confidence intervals}
\usage{
extract_NRI_CI(results.boot, conf.level, n.boot, dp)
}
\arguments{
\item{results.boot}{The matrix of n.boot metrics from within CI.NRI}

\item{conf.level}{The confidence interval expressed between 0 & 1 (eg 95\%CI is conf.level = 0.95)}

\item{n.boot}{The number of bootstrapped samples}

\item{dp}{the number of decimal places to report the point estimate and confidence interval}
}
\value{
A two column matrix with the metric name and statistic with a confidence interval
}
\description{
Extract a confidence in interval from the bootstrapped results. Used by CI.NRI
}
