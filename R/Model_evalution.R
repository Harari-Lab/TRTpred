# Model Evaluation
# Author: Rémy Pétremand
# Date: 07.05.2024
# Description: Function to evaluate the performance of models
# Reference: https://doi.org/10.1038/s41587-024-02232-0

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Global Parameters
# ------------------------------------------------------------------------------

EVALUATION.METRICS <- c("mcc", "accuracy", "F1", "kappa", "auc", "sensitivity", "specificity", "PPV", "NPV")

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

#' Get Binary Accuracy
#' 
#' Function that compute accuracy based on a threshold
#' 
#' @param x numerical predictions
#' @param y logical ground truth
#' @param th numerical threshold to apply on x to get logical predictions
#' @param weights vector of numerical of length identical to x. 
#' It represent the weights of the samples
#' Default vector of ones 
#' 
#' @return accuracy for x for a threshold th with samples weighted by weights
.GetBinaryScoreAccuracy <- function(x, y, th, weights = replicate(n = length(x), expr = 1)){
  pred <-  (x >= th) == y
  sum.T <- sum(weights[pred])
  sum.all <- sum(weights)
  return(sum.T/sum.all)
}

#' Get Binary Prediction
#' 
#' Function to get the prediction for some numerical score. In brief, a threshold 
#' is applied on a score vector to get binary prediction which are returned by 
#' the function
#' The threshold can be provided by the user (method = "threshold") or the 
#' threshold can be computed in a "greedy" manner by optimizing the accuracy. 
#' So there are two methods: 
#' - "threshold": unbiased method where `x.threshold` must be provided
#' - "greedy": Optimal threshold is computed to optimize maxium accuracy. Here, 
#' 'ground.truth' needs to be provided. 
#' 
#' @param x.train numerical vector; The training score vector
#' @param method character; The threhsold computation method. Possible values 
#' are "greedy" and "threshold".
#' Default = "greedy"
#' @param ground.truth logical vector; The logical ground truth needed for when method = "greedy"
#' @param x.threshold numerical; The threshold input for when method = "threshold"
#' @param x.test numerical vector; The testing score vector. If NULL, return no testing prediction
#' Default = NULL
#' @param weights numerical vector of length identical to x.train; The observation weights.
#' Default = NULL i.e. 1 for each observation
#' 
#' @return list("y.train.pred" = ..., "y.test.pred" = ...)
#' 
#' @export
BinarizeScorePrediction <- function(x.train, 
                                        method = c("greedy", "threshold"), 
                                        ground.truth = NULL,
                                        x.threshold = NULL,
                                        x.test = NULL, 
                                        weights = NULL){
  method <- match.arg(method)
  
  # Clean input if there is NAs in the the ground.truths
  if (any(is.na(ground.truth))){
    mask.keep <- !is.na(ground.truth)
    ground.truth <- ground.truth[mask.keep]
    x.train <- x.train[mask.keep]
    weights <- weights[mask.keep]
  }
  
  # sanity checks + control inputs
  if (method == "threshold" & is.null(x.threshold)){
    warning("BinarizeScorePrediction: Method is 'threshold' and x.threshold is not defined. Continue with x.threshold = 0")
    x.threshold <- 0
  } else if (method == "greedy" & is.null(ground.truth)){
    if (is.null(x.threshold)){
      x.threshold <- 0
    }
    method <- "threshold"
    warning(paste0("BinarizeScorePrediction: Method is 'greedy' and ground.truth is not defined. Continue with method = 'threshold' and x.threshold = ", x.threshold))
  }
  
  if (method == "greedy"){
    range_ <- range(x.train, na.rm = T)
    if ((range_[2]-range_[1]) > 0){
      if (is.null(weights)){
        weights <- replicate(n = length(x.train), expr = 1)
      }
      opt.res <- stats::optimize(
        f = .GetBinaryScoreAccuracy, 
        interval = range_, 
        x = x.train, 
        y = ground.truth, 
        weights = weights,
        maximum = T)
      x.threshold <- opt.res$maximum
    } else {
      x.threshold <- 0
    }
    y.train.pred <- x.train >= x.threshold
  } else if (method == "threshold"){
    y.train.pred <- x.train >= x.threshold
  }
  
  if (!is.null(x.test)){
    y.test.pred <- x.test >= x.threshold
  } else {
    y.test.pred <- NULL
  }
  
  return(list("y.train.pred" = y.train.pred, "y.test.pred" = y.test.pred, "x.threshold" = x.threshold))
}

#' Get binary predictive metric
#' 
#' Function to get TP, TN, FP & FN
#' 
#' @param ground_truth character vector of ground truths
#' @param preds character vector of predictions
#' 
#' @return list; list composed of the number of TP, TP, FP & FN
#' 
#' @export
GetMetricsBinary <- function(ground_truth, preds){
  
  ground_truth <- c(ground_truth)
  preds <- c(preds)
  
  if (length(ground_truth) != length(preds)){
    stop("ground_truth and pred not the same length")
  }
  
  tp <- sum( ground_truth  &  preds, na.rm = T)
  tn <- sum(!ground_truth  & !preds, na.rm = T)
  fp <- sum(!ground_truth  &  preds, na.rm = T)
  fn <- sum( ground_truth  & !preds, na.rm = T)
  
  return(list("TP" = tp, "FP" = fp, "TN" = tn, "FN" = fn))
}

#' Get binary accuracy metric
#' 
#' Function to get a list of binary accuracy metrics
#' 
#' @param metrics character; The accuracy metrics. can be one or multiple metrics 
#' in "accuracy", "mcc", "F1", "kappa"
#' @param ground_truth character vector of ground truths
#' @param preds character vector of predictions
#' @param tp numerical; True Positives
#' @param tn numerical; True Negatives
#' @param fp numerical; False Positives
#' @param fn numerical; False Negatives
#' 
#' @details 
#' The user can provide ground_truth and preds and the function will compute tp, 
#' tn, fp, fn using the function GetMetricsBinary(). 
#' The user can also only provide tp, tn, fp, fn directly
#' 
#' @return list; list of accuracy measurements
#' 
#' @export
GetAccuracyMetrics <- function(metrics = EVALUATION.METRICS,
                               ground_truth = NULL, preds = NULL, scores = NULL,
                               tp = NULL, fp = NULL, tn = NULL, fn = NULL){
  
  
  if (!all(metrics %in% EVALUATION.METRICS)){
    stop(paste0("GetAccuracyMetrics: Invalild metrics provided. Possible values are: ", paste(EVALUATION.METRICS, collapse = ", ")))
  }
  
  if (is.null(tp) & !is.null(preds) & !is.null(ground_truth)){
    res <- GetMetricsBinary(ground_truth = ground_truth, preds = preds)
    tp <- res$TP
    fp <- res$FP
    tn <- res$TN
    fn <- res$FN
  } 
  
  N <- (tp + fp + tn + fn)
  
  res <- list("TP" = tp, "FP" = fp, "TN" = tn, "FN" = fn)
  
  if ("accuracy" %in% metrics){
    res[["accuracy"]] <- (tp + tn) / N
  }
  if ("mcc" %in% metrics){ 
    # References:
    # https://en.wikipedia.org/wiki/Phi_coefficient
    # https://github.com/davidechicco/MCC/blob/master/bin/confusion_matrix_rates.r
    if (N > 0){
      sum1 <- tp+fp; sum2 <-tp+fn ; sum3 <-tn+fp ; sum4 <- tn+fn
      denominator <- as.double(sum1)*sum2*sum3*sum4 # as.double to avoid overflow error on large products
      if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {
        denominator <- 1
      }
      denominator <- sqrt(denominator)
      nominator <- (tp*tn)-(fp*fn)
      res[["mcc"]] <- nominator / denominator
    } else {
      res[["mcc"]] <- NA
    }
    
  }
  if ("F1" %in% metrics){
    res[["F1"]] <- (2*tp) / ((2*tp) + fp + fn)
  }
  if ("kappa" %in% metrics){
    res[["kappa"]] <- (tp + tn) / N
  }
  if ("sensitivity" %in% metrics){
    res[["sensitivity"]] <- tp / (tp + fn)
  }
  if ("specificity" %in% metrics){
    res[["specificity"]] <- tn / (tn + fp)
  }
  if ("PPV" %in% metrics){
    res[["PPV"]] <- tp / (tp + fp)
  }
  if ("NPV" %in% metrics){
    res[["NPV"]] <- tn / (tn + fn)
  }
  if ("auc" %in% metrics & !is.null(scores)){
    scores[is.infinite(scores)] <- NA
    if (!(all(is.na(ground_truth)) | all(is.na(scores)))){
      if (length(unique(ground_truth)) > 1){
        roc.res <- pROC::roc(response = ground_truth, 
                             predictor = scores, 
                             direction = "<", 
                             levels = c(F, T), 
                             auc = T)
        res[["auc"]] <- as.numeric(roc.res$auc)
      } else {
        res[["auc"]] <- as.numeric(NA)
      }
    } else{
      res[["auc"]] <- as.numeric(NA)
    }
  }
  
  return(res)
}

#' Get mcc metric
#' 
#' Function to get the mcc metric.
#' 
#' @param ground_truth character vector of ground truths
#' @param preds character vector of predictions
#' @param tp numerical; True Positives
#' @param tn numerical; True Negatives
#' @param fp numerical; False Positives
#' @param fn numerical; False Negatives
#' 
#' @details 
#' The user can provide ground_truth and preds and the function will compute tp, 
#' tn, fp, fn using the function GetMetricsBinary(). 
#' The user can also only provide tp, tn, fp, fn directly
#' 
#' @return numerical; The mcc
#' 
#' @export
GetMccBinary <- function(ground_truth = NULL, preds = NULL, 
                         tp = NULL, fp = NULL, tn = NULL, fn = NULL){
  
  if (is.null(tp)){
    res <- GetMetricsBinary(ground_truth = ground_truth, preds = preds)
    tp <- res$TP
    fp <- res$FP
    tn <- res$TN
    fn <- res$FN
  }
  
  # https://github.com/davidechicco/MCC/blob/master/bin/confusion_matrix_rates.r
  sum1 <- tp+fp; sum2 <-tp+fn ; sum3 <-tn+fp ; sum4 <- tn+fn
  denominator <- as.double(sum1)*sum2*sum3*sum4 # as.double to avoid overflow error on large products
  if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {
    denominator <- 1
  }
  denominator <- sqrt(denominator)
  nominator <- (tp*tn)-(fp*fn)
  res[["mcc"]] <- nominator / denominator
  
  return(res)
}
#' Get accuracy metric
#' 
#' Function to get the accuracy metric
#' 
#' @param ground_truth character vector of ground truths
#' @param preds character vector of predictions
#' @param tp numerical; True Positives
#' @param tn numerical; True Negatives
#' @param fp numerical; False Positives
#' @param fn numerical; False Negatives
#' 
#' @details 
#' The user can provide ground_truth and preds and the function will compute tp, 
#' tn, fp, fn using the function GetMetricsBinary(). 
#' The user can also only provide tp, tn, fp, fn directly
#' 
#' @return numerical; The accuracy
#' 
#' @export
GetAccurayBinary <- function(ground_truth = NULL, preds = NULL, 
                             tp = NULL, fp = NULL, tn = NULL, fn = NULL){
  
  if (is.null(tp)){
    res <- GetMetricsBinary(ground_truth = ground_truth, preds = preds)
    tp <- res$TP
    fp <- res$FP
    tn <- res$TN
    fn <- res$FN
  }
  
  res <- (tp + tn) / (tp + fp + tn + fn)
  
  return(res)
}

#' Get F1 metric
#' 
#' Function to get the F1 metric
#' 
#' @param ground_truth character vector of ground truths
#' @param preds character vector of predictions
#' @param tp numerical; True Positives
#' @param tn numerical; True Negatives
#' @param fp numerical; False Positives
#' @param fn numerical; False Negatives
#' 
#' @details 
#' The user can provide ground_truth and preds and the function will compute tp, 
#' tn, fp, fn using the function GetMetricsBinary(). 
#' The user can also only provide tp, tn, fp, fn directly
#' 
#' @return numerical; The F1
#' 
#' @export
GetF1Binary <- function(ground_truth = NULL, preds = NULL, 
                        tp = NULL, fp = NULL, tn = NULL, fn = NULL){
  
  if (is.null(tp)){
    res <- GetMetricsBinary(ground_truth = ground_truth, preds = preds)
    tp <- res$TP
    fp <- res$FP
    tn <- res$TN
    fn <- res$FN
  }
  
  res <- (2*tp) / ((2*tp) + fp + fn)
  
  return(res)
}

#' Get kappa metric
#' 
#' Function to get the kappa metric
#' 
#' @param ground_truth character vector of ground truths
#' @param preds character vector of predictions
#' @param tp numerical; True Positives
#' @param tn numerical; True Negatives
#' @param fp numerical; False Positives
#' @param fn numerical; False Negatives
#' 
#' @details 
#' The user can provide ground_truth and preds and the function will compute tp, 
#' tn, fp, fn using the function GetMetricsBinary(). 
#' The user can also only provide tp, tn, fp, fn directly
#' 
#' @return numerical; The kappa
#' 
#' @export
GetKappaBinary <- function(ground_truth = NULL, preds = NULL, 
                           tp = NULL, fp = NULL, tn = NULL, fn = NULL){
  
  if (is.null(tp)){
    res <- GetMetricsBinary(ground_truth = ground_truth, preds = preds)
    tp <- res$TP
    fp <- res$FP
    tn <- res$TN
    fn <- res$FN
  }
  
  res <- (tp + tn) / (tp + fp + tn + fn)
  
  return(res)
}

#' Get AUC metrics
#' 
#' Function to get the AUC metric
#' 
#' @param ground_truth vector of ground truths.
#' @param scores numerical vector of scores
#' 
#' @return numerical; The AUC
#' 
#' @export
GetAUC <- function(ground_truth, scores){
  
  ground.truth <- as.factor(ground.truth)
  
  roc.res <- pROC::roc(
    response = ground.truth, 
    predictor = scores, 
    direction = ">", 
    levels = levels(ground.truth), 
    auc = T)
  
  res <- as.numeric(roc.res$auc)
  
  return(res)
}

#' Get sensitivity metric
#' 
#' Function to get the sensitivity metric
#' 
#' @param ground_truth character vector of ground truths
#' @param preds character vector of predictions
#' @param tp numerical; True Positives
#' @param tn numerical; True Negatives
#' @param fp numerical; False Positives
#' @param fn numerical; False Negatives
#' 
#' @details 
#' The user can provide ground_truth and preds and the function will compute tp, 
#' tn, fp, fn using the function GetMetricsBinary(). 
#' The user can also only provide tp, tn, fp, fn directly
#' 
#' @return numerical; The sensitivity
#' 
#' @export
GetSensitivityBinary <- function(ground_truth = NULL, preds = NULL, 
                           tp = NULL, fp = NULL, tn = NULL, fn = NULL){
  
  if (is.null(tp)){
    res <- GetMetricsBinary(ground_truth = ground_truth, preds = preds)
    tp <- res$TP
    fp <- res$FP
    tn <- res$TN
    fn <- res$FN
  }
  
  res <- (tp) / (tp + fn)
  
  return(res)
}

#' Get specificity metric
#' 
#' Function to get the specificity metric
#' 
#' @param ground_truth character vector of ground truths
#' @param preds character vector of predictions
#' @param tp numerical; True Positives
#' @param tn numerical; True Negatives
#' @param fp numerical; False Positives
#' @param fn numerical; False Negatives
#' 
#' @details 
#' The user can provide ground_truth and preds and the function will compute tp, 
#' tn, fp, fn using the function GetMetricsBinary(). 
#' The user can also only provide tp, tn, fp, fn directly
#' 
#' @return numerical; The specificity
#' 
#' @export
GetSpecificityBinary <- function(ground_truth = NULL, preds = NULL, 
                                 tp = NULL, fp = NULL, tn = NULL, fn = NULL){
  
  if (is.null(tp)){
    res <- GetMetricsBinary(ground_truth = ground_truth, preds = preds)
    tp <- res$TP
    fp <- res$FP
    tn <- res$TN
    fn <- res$FN
  }
  
  res <- (tn) / (tn + fp)
  
  return(res)
}
