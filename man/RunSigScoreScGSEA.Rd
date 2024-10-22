% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Signature_score_function.R
\name{RunSigScoreScGSEA}
\alias{RunSigScoreScGSEA}
\title{Run the scGSEA signature score}
\usage{
RunSigScoreScGSEA(
  X,
  signature,
  ranks = NULL,
  to.scale = T,
  scale.params = NULL
)
}
\arguments{
\item{X}{matrix; The matrix. Columns = cells, rows = some features}

\item{signature}{list: signature as defined in GetSignature()}

\item{ranks}{rank-matrix as defined by matrixStats::colRanks.
If provided, the function does not compute the rankings from X}

\item{to.scale}{logical: Weather to scale signature score or not.}
}
\value{
list with three elements: 
score: data.frame; The signature score data.frame. Row = cells and columns = names(signature)
scale.mean: mean value to scale (if to.scale)
scale.sd: sd value to scale (if to.scale)
}
\description{
Function to run the scGSEA signature score
}
