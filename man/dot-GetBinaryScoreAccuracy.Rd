% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model_evalution.R
\name{.GetBinaryScoreAccuracy}
\alias{.GetBinaryScoreAccuracy}
\title{Get Binary Accuracy}
\usage{
.GetBinaryScoreAccuracy(x, y, th, weights = replicate(n = length(x), expr = 1))
}
\arguments{
\item{x}{numerical predictions}

\item{y}{logical ground truth}

\item{th}{numerical threshold to apply on x to get logical predictions}

\item{weights}{vector of numerical of length identical to x. 
It represent the weights of the samples
Default vector of ones}
}
\value{
accuracy for x for a threshold th with samples weighted by weights
}
\description{
Function that compute accuracy based on a threshold
}
