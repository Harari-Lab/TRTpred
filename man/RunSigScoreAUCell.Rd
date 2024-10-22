% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Signature_score_function.R
\name{RunSigScoreAUCell}
\alias{RunSigScoreAUCell}
\title{Run the AUCell signature score}
\usage{
RunSigScoreAUCell(
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

\item{ranks}{rank-matrix as defined by AUCell::AUCell_buildRankings(). 
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
Function to run the AUCell signature score
}
