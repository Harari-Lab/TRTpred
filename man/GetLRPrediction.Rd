% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Logistic_regression.R
\name{GetLRPrediction}
\alias{GetLRPrediction}
\title{Get Prediction for the LR model}
\usage{
GetLRPrediction(
  x,
  path.folder,
  hyperparameters = DEFAULT.LR.HYPERPARAMS,
  model.file.name = NULL,
  prob.threshold = 0.5
)
}
\arguments{
\item{x}{data.frame/matrix. The input data with the dependent and independent variables.}

\item{path.folder}{character; Folder to retrieve the model}

\item{hyperparameters}{list; The LR hyperparmeters
Elements of list are 
  - "alpha" (L2-regularization): Vector of alpha values (0-1)
  - "lambda" (L1-regularization): Vector of lambda values (0-1)}

\item{model.file.name}{character; the file name of the model in the folder 
if null retrive it from "res_fit_"<alpha>"_"<lambda>".rds"
Default = NULL}

\item{prob.threshold}{numerical; The LR probability threshold
Default = 0.5}
}
\value{
data.frame; Rows = barcode names in x. Cols = "score" and "pred"
}
\description{
The function to get prediction from a LR model
}
