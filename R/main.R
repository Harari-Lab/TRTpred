# TRTpred main
# Author: Rémy Pétremand
# Date: 07.05.2024
# Description: main functions for training and testing TRTpred
# Reference: https://doi.org/10.1038/s41587-024-02232-0

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

suppressMessages(require(Seurat))
suppressMessages(require(foreach))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))

# ------------------------------------------------------------------------------
# Global Parameters
# ------------------------------------------------------------------------------

DEFAULT.LR.HYPERPARAMS <- list("alpha" = 0, "lambda" = 0)
DEFAULT.SIGNATURE.HYPERPARAMS <- list("signature.lengths" = c(20), 
                                      "signature.side" = c("both"),
                                      "signature.rm.regex" = "none",
                                      "signature.methods" = c("AUCell"),
                                      "signature.selection.method" = c("logFC"))

P.ADJUST.METHODS <- c("fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "none")
LIMMA.METHODS <- c("limma_voom", "limma_trend")
EDGER.METHODS <- c("edgeR_LRT", "edgeR_QFL")
DESEQ.METHODS <- c("DESeq2_Wald", "DESeq2_LRT")
SEURAT.METHODS <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR")
DEA.METHODS <- c(SEURAT.METHODS, DESEQ.METHODS, LIMMA.METHODS, EDGER.METHODS)
FEATURE.TRANS.METHODS <- c("pca", "opls")
EVALUATION.METRICS <- c("mcc", "accuracy", "F1", "kappa", "auc", "sensitivity", "specificity", "PPV", "NPV")

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

#' Get Models prediction
#' 
#' Home function to apply TRTpred
#' 
#' @param seurat.obj Seurat Object; RNA-sequencing single-cell object (required)
#' @param seurat.assay character; The name of the seurat assay. 
#' Default = "RNA"
#' @param seurat.slot character; The name of the seurat slot
#' Default = "scale.data"
#' @param covariate.cols list of character; covariates
#' Default = NULL
#' @param y.label character; The y-label name. 
#' Default = NULL
#' 
#' @return data.frame; Cell-Wise TRT predictions 
#' Rows = barcode names in x. Cols = "score" and "pred"
#' 
#' @export
GetTRTpred <- function(
    seurat.obj, seurat.assay = "RNA", seurat.slot = "scale.data",
    covariate.cols = NULL, y.label = NULL){
  
  # Get input data: 
  data.input <- Seurat::GetAssayData(object = seurat.obj,
                                     assay = seurat.assay, 
                                     slot = seurat.slot)
  data.input <- data.frame(t(data.input), check.rows = F, check.names = F)
  
  # Get parameters: 
    
  data("best.model.hyperparameters", package = "TRTpred")
  
  method <- best.model.hyperparameters$method; 
  
  data.trans.method <- best.model.hyperparameters$data.trans.method
  if(data.trans.method == "none"){
    data.output <- NULL
    y.label <- NULL
  } else {
    data.output <- seurat.obj@meta.data[, unique(c(y.label, covariate.cols))]
  }

  PrepareData.list <- PrepareData(
    x.train = data.input, 
    y.train = data.output, 
    y.label = y.label, 
    method = method,
    data.trans.method = data.trans.method, 
    rm.corr.features = F, 
    DA.method = "none")
  
  if (method == "signature") {
    
    DEA.res <- best.model.hyperparameters$DEA
    best.model.hyperparameters$DEA <- NULL
    
    signature.list <- GetSignatureList(
      df = DEA.res,
      pval_col = "padj", 
      log2FC_col = "logFC", 
      pval_limit = 0.05, 
      log2FC_limits = c(0, 0), 
      gene.selection.methods = best.model.hyperparameters$signature.selection.method,
      lengths = best.model.hyperparameters$signature.lengths, 
      sides = best.model.hyperparameters$signature.side,
      remove_regex = best.model.hyperparameters$signature.rm.regex)
    
    prediction.res <- GetPredictionSignature(
      x = PrepareData.list$data.train, 
      path.folder = NULL, 
      DEA.file.name = NULL, 
      signature.list = signature.list,
      hyperparameters = best.model.hyperparameters, 
      score.threshold = best.model.hyperparameters$signature.x.threshold, 
      score.to.scale = !is.null(best.model.hyperparameters$signature.score.scale.params), 
      score.scale.params = best.model.hyperparameters$signature.score.scale.params)
    
  } else if (method == "LR") {
    prediction.res <- GetLRPrediction(
      x = PrepareData.list$data.train, 
      path.folder = path.folder, 
      model.file.name = LR.model.file.name, 
      prob.threshold = LR.prob.threhsold,
      hyperparameters = hyperparameters)
  }
  
  return(prediction.res)
}

#' Get clone-wise prediction
#' 
#' Get clone wise prediction from cell-wise TRT predictions. 
#' 
#' @param pred.df data.frame; The prediction data-frame with 3 required columns: 
#' - the score column (score.col)
#' - the prediction column (pred.col)
#' - The clone-id column (clone.col)
#' @param score.col character; The numerical score column name in pred.df
#' Default = "score"
#' @param pred.col character; The boolean pred column name in pred.df
#' Default = "pred"
#' @param clone.col character; The character clone-id column in pred.df
#' Default = "CloneId"
#' 
#' @return data.frame; Rows = clones id. Cols = "score.clone" and "pred.clone" and number of cells
#' 
#' @export
GetClonePred <- function(pred.df, score.col = "score", pred.col = "pred", clone.col = "CloneId"){
  
  pred.df$score <- pred.df[[score.col]]
  pred.df$pred <- pred.df[[pred.col]]
  
  pred.clone.df <- pred.df %>% 
    group_by(.dots = clone.col) %>% summarise(
      pred.clone = any(pred, na.rm = T),
      score.clone = max(score, na.rm = T),
      n.cells = n(), 
      n.cells.TRT = sum(pred, na.rm = T)
    ) %>%
    as.data.frame()
  
  return(pred.clone.df)
}

#' Get Models prediction
#' 
#' Home function to train and test a model (LR or signature )
#' 
#' @param data.input data.frame or matrix; The input matrix corresponding to the feature 
#' space only. Additional covariates are found in `y`. Rows = samples, Columns = features
#' @param data.output data.frame or matrix; The output matrix. 
#' Columns = c(y.label, y.covariates, y.sample)
#' @param y.label character; The column of `y` corresponding to the binary output
#' @param path.folder character; Folder to retieve the model / data / info
#' @param method character; The prediction method. Possible values are: "LR": 
#' Logistic regression or "signature" signature-score approach. 
#' @param hyperparams list: The hyper-parameter list. 
#' If method = "LR": Elements of the list are 
#'   - "alpha" (L2-regularization): Vector of alpha values (0-1).
#'   - "lambda" (L1-regularization): Vector of lambda values (0-1).
#' If method = "signature": Elements of the list are 
#'   - "signature.lengths": Vector of signature lengths.
#'   - "signature.side": Vector of signature sides ("up", "down", or "both").
#'   - "signature.rm.regex": Vector of rm.regex (regex to remove genes from signature).
#'   - "signature.methods": Vector of signature score methods.
#'   - "signature.selection.method": How to select genes from DEA res using the log-Fold-Change ("logFC") or p-values ("pval")
#' @param data.trans.method character: The data transformation method. 
#' Possible values are "pca" 
#' Default = "none" i.e. no transformation. 
#' @param DEA.file.name character; The name of the DEA file 
#' Default = "dea.res.rds"
#' @param signature.list list; The signature list
#' Default = NULL
#' @param signature.x.threshold numerical; The signature score threshold
#' Default = 0
#' @param LR.model.file.name character; The name of the LR file. 
#' If NULL, generate the file name generated by LRCrossValidation()
#' Default = NULL
#' @param LR.prob.threhsold numeric; The probability threshold for the logistic regression
#' Default = 0.5
#' 
#' @return data.frame; Rows = barcode names in x. Cols = "score" and "pred"
#' 
#' @export
GetModelPrediction <- function(data.input, data.output, y.label, 
                               path.folder, 
                               method, hyperparameters,
                               data.trans.method, 
                               DEA.file.name = "dea.res.rds", signature.list = NULL, signature.x.threshold = 0,
                               signature.score.scale = T, signature.score.scale.params = NULL,
                               LR.model.file.name = NULL, LR.prob.threhsold = NULL){
  
  PrepareData.list <- PrepareData(
    x.train = data.input, 
    y.train = data.output, 
    y.label = y.label, 
    method = method,
    data.trans.method = data.trans.method, 
    rm.corr.features = F, 
    DA.method = "none")
  
  if (method == "signature") {
    prediction.res <- GetPredictionSignature(
      x = PrepareData.list$data.train, 
      path.folder = path.folder, 
      DEA.file.name = DEA.file.name, 
      signature.list = signature.list,
      hyperparameters = hyperparameters, 
      score.threshold = signature.x.threshold, 
      score.to.scale = signature.score.scale, 
      score.scale.params = signature.score.scale.params)
    
  } else if (method == "LR") {
    prediction.res <- GetLRPrediction(
      x = PrepareData.list$data.train, 
      path.folder = path.folder, 
      model.file.name = LR.model.file.name, 
      prob.threshold = LR.prob.threhsold,
      hyperparameters = hyperparameters)
  }
  
  return(prediction.res)
}

#' Train Model function
#' 
#' Home function to train and test a model (LR or signature )
#' 
#' @param x.train data.frame or matrix; The input matrix corresponding to the feature 
#' space only. Additional covariates are found in `y`. Rows = samples, Columns = features
#' @param y.train data.frame or matrix; The output matrix. 
#' Columns = c(y.label, y.covariates, y.sample)
#' @param x.test data.frame or matrix; The input data like x.train but for testing
#' @param y.test data.frame or matrix; The output matrix like y.train but for testing
#' @param y.label character; The column of `y` corresponding to the binary output
#' @param y.covariates character vector; The columns of `y` corresponding to the covariates
#' @param y.sample character; The column of `y` corresponding to the samples for when method = "signature"
#' @param method character; The prediction method. Possible values are: "LR": 
#' Logistic regression or "signature" signature-score approach. 
#' @param hyperparams list: The hyper-parameter list. 
#' If method = "LR": Elements of the list are 
#'   - "alpha" (L2-regularization): Vector of alpha values (0-1).
#'   - "lambda" (L1-regularization): Vector of lambda values (0-1).
#' If method = "signature": Elements of the list are 
#'   - "signature.lengths": Vector of signature lengths.
#'   - "signature.side": Vector of signature sides ("up", "down", or "both").
#'   - "signature.rm.regex": Vector of rm.regex (regex to remove genes from signature).
#'   - "signature.methods": Vector of signature score methods.
#'   - "signature.selection.method": How to select genes from DEA res using the log-Fold-Change ("logFC") or p-values ("pval")
#' @param data.trans.method character: The data transforamtion method. 
#' Possible values are "pca" 
#' Default = "none" i.e. no transformation. 
#' @param pca.explained.var.threshold numeric: If not NULL, apply this threshold 
#' to select PCs explaining more variance than the threhsold
#' @param rm.corr.features logical: Do we remove the mutually correlated feature 
#' while keeping the ones that best correlate with the outcome? 
#' default = F
#' @param DA.method character; The discriminant analysis (DA) method. 
#' Possible values are "wilcox"
#' Default = "none" i.e. no discriminant analysis
#' @param DA.p.adjust.method character; The p-value adjust method for the 
#' discriminant analysis. 
#' Default = "fdr"
#' @param DA.p.val.threshold numeric; The p-val significant threshold. 
#' @param DA.event.per.variable, numerical; The number of event per variable threshold
#' A common EPV is 10 meaning that we need 10 event per variable. Here we use it as 
#' a way to select a fewer amount of variables for a fixed amount of samples. 
#' Default = NULL (no selection on EPV)
#' @param LR.prob.threhsold numeric; The probability threshold for the logitic regression
#' @param signature.col.aggregate character; Column in data.train@meta.data to combine cells into pseudobulk. 
#' If NULL, no pseudobulk integration
#' Default = NULL
#' @param DEA.method character; The Differential-expression-analylsis method (see RunDEA())
#' default = "wilcox" from the Seurat::FindMarker() function
#' @param DEA.data data.frame or matrix; The input matrix for the DEA. 
#' Rows = samples, Columns = features. If NULL, the DEA.data is x
#' Default = NULL
#' @param sample.weights data.frame (rows = observation, column = weight); The 
#' observation weights data.frame. 
#' Default = NULL i.e. 1 for each observation
#' @param folds.save.folder character; Folder to save the model information. 
#' If null doesn't save anything. 
#' default = NULL
#' 
#' @return data.frame summarizing the Nested-cross-validation
#' 
#' @export
TrainModel <- function(x.train, y.train, 
                       y.label, y.covariates = NULL, y.sample = NULL,
                       x.test = NULL, y.test = NULL,
                       method = c("LR", "signature"),
                       hyperparams = NULL,
                       data.trans.method = c("none", FEATURE.TRANS.METHODS), 
                       pca.explained.var.threshold = NULL,
                       rm.corr.features = F,
                       DA.method = c("none", DEA.METHODS),
                       DA.p.adjust.method = P.ADJUST.METHODS,
                       DA.p.val.threshold = 0.05,
                       DA.event.per.variable = NULL,
                       LR.prob.threhsold = NULL,
                       signature.col.aggregate = NULL,
                       DEA.method = DEA.METHODS, 
                       DEA.data = NULL,
                       sample.weights = NULL,
                       folds.save.folder = NULL,
                       save.fold.model = F, 
                       save.fold.pred = F, 
                       save.fold.score = F){
  
  # Setting parameters: 
  method <- match.arg(method) # default = "LR"
  data.trans.method <- match.arg(data.trans.method) # default = "none"
  DA.method <- match.arg(DA.method) # default = "none"
  DA.p.adjust.method <- match.arg(DA.p.adjust.method) # default = "fdr"
  DEA.method <- match.arg(DEA.method) # default = "wilcox"
  
  if (!is.null(x.test)){
    if (nrow(x.test) == 0){
      x.test <- NULL
      y.test <- NULL
    }
  }
  
  data.input.list <- PrepareData(
    x.train = x.train,
    y.train = y.train,
    y.label = y.label,
    x.test = x.test,
    y.test = y.test,
    DEA.data = DEA.data,
    method = method,
    data.trans.method = data.trans.method,
    pca.explained.var.threshold = pca.explained.var.threshold,
    rm.corr.features = rm.corr.features,
    rm.corr.features.threhsold = 0.8,
    DA.method = DA.method,
    y.sample = y.sample, 
    y.covariates = y.covariates, 
    y.aggregate = signature.col.aggregate,
    DA.event.per.variable = DA.event.per.variable, 
    DA.p.val.threshold = DA.p.val.threshold,
    DA.p.adjust.method = DA.p.adjust.method)
  
  
  if (method == "LR"){
    
    # Define the design of the model
    col.features <- colnames(data.input.list$data.train)
    col.features <- col.features[!(col.features %in% colnames(y.train))]
    design.str <- paste(y.label, "~", paste(c(col.features, y.covariates), collapse = " + "), sep = " ")
    
    # Run the training
    res.model <- LRCrossValidation(
      data.train = data.input.list$data.train, 
      data.test = data.input.list$data.test,
      design.str = design.str, 
      hyperparams = hyperparams, 
      model.prob.threshold = LR.prob.threhsold,
      sample.weights = sample.weights,
      folds.save.folder = folds.save.folder,
      save.model = save.fold.model, 
      save.pred = save.fold.pred, 
      save.score = save.fold.score
      )
    
  } else if (method == "signature"){
    # Get the results of the model for all hyperparmeters
    res.model <- SignatureCrossValidation(
      data.train = data.input.list$data.train, 
      data.test = data.input.list$data.test, 
      DEA.data = data.input.list$DEA.data,
      y.label = y.label, 
      y.sample = y.sample, 
      y.covariates = y.covariates, 
      DEA.method = DEA.method,
      signature.selection.method = hyperparams$signature.selection.method,
      signature.lengths = hyperparams$signature.lengths, 
      signature.sides = hyperparams$signature.side, 
      signature.rm.regex = hyperparams$signature.rm.regex, 
      signature.methods = hyperparams$signature.methods, 
      assay = "RNA", # in data.train assay "RNA" contains now the correct set of data
      slot = "data", # in data.train slot "data" contains now the correct set of data
      col.aggregate = signature.col.aggregate,
      sample.weights = sample.weights,
      folds.save.folder = folds.save.folder,
      save.model = save.fold.model, 
      save.pred = save.fold.pred, 
      save.score = save.fold.score)
  }
  
  res.model[["n.data.input"]] <- nrow(x.train)
  res.model[["p.data.input"]] <- ncol(x.train)
  
  return(res.model)
}
