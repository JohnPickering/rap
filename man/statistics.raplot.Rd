% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rap.R
\name{statistics.raplot}
\alias{statistics.raplot}
\title{Statistical metrics}
\usage{
statistics.raplot(x1, x2, y, t = NULL, NRI_return = FALSE)
}
\arguments{
\item{x1}{Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the baseline model.  Must be between 0 & 1}

\item{x2}{Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1}

\item{y}{Binary of outcome of interest. Must be 0 or 1 (if fitted models are provided this is extracted from the fit which for an rms fit must have x = TRUE, y = TRUE).}

\item{t}{The risk threshold(s) for groups. eg t<-c(0,0.1,1) is a two group scenario with a threshold of 0.1 & t<-c(0,0.1,0.3,1)  is a three group scenario with thresholds at 0.1 and 0.3. Nb. If no t is provided it defaults to a single threshold at the prevalence of the cohort.}

\item{NRI_return}{Flag to return NRI metrics, default is FALSE.}
}
\value{
A matrix of metrics for use within CI.raplot
}
\description{
The function statistics.raplot calculates the reclassification metrics. Used by CI.raplot.
}
