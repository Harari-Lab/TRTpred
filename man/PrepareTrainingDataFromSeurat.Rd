% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data_preparation.R
\name{PrepareTrainingDataFromSeurat}
\alias{PrepareTrainingDataFromSeurat}
\title{Prepare training data from seurat object}
\usage{
PrepareTrainingDataFromSeurat(
  SeuratObject,
  label,
  label.order = NULL,
  assay = "RNA",
  slot = "data",
  DEA.assay = NULL,
  DEA.slot = "data",
  DEA.pseudobulk.col = NULL,
  covariates = NULL,
  sample.weights.col = NULL
)
}
\arguments{
\item{SeuratObject}{Seurat object; The Seurat object (required)}

\item{label}{character; The y-label (i.e tumor reactivtiy) (required)}

\item{label.order}{character; the order of the y-label to transform them into factor}

\item{assay}{character; The Seurat assay (default = "RNA")}

\item{slot}{character; The assay slot (default = "slot")}

\item{DEA.assay}{character; The Seurat assay used for the Differential Expression Analysis 
Default = NULL = same as assay}

\item{DEA.slot}{character; The assay slot used for the Differential Expression Analysis 
Default = NULL = same as slot}

\item{DEA.pseudobulk.col}{character; The column names in Seurat used to construct pseudobluk for the pseudobluk DEA
Default = NULL}

\item{covariates}{characters; The covariates list
Default = NULL}

\item{sample.weights.col}{character; Which col should be used to normalize the performance
Default = NULL = none!}
}
\value{
list; 
- "data.input" 
- "data.output" 
- "DEA.data" 
- "sample.weights"
}
\description{
This function prepare training data from Seurat object
}
