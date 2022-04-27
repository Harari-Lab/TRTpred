# Model evaluation functions
# Rémy Pétremand

EVALUATION.METRICS <- c("accuracy", "mcc", "F1", "kappa")

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
  
  if (length(ground_truth) != length(preds)){
    stop("ground_truth and pred not the same length")
  }
  
  tp <- sum( ground_truth  &  preds)
  tn <- sum(!ground_truth  & !preds)
  fp <- sum(!ground_truth  &  preds)
  fn <- sum( ground_truth  & !preds)
  
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
                               ground_truth = NULL, preds = NULL, 
                               tp = NULL, fp = NULL, tn = NULL, fn = NULL){
  
  
  if (!all(metrics %in% EVALUATION.METRICS)){
    stop(paste0("GetAccuracyMetrics: Invalild metrics provided. Possible values are: ", paste(EVALUATION.METRICS, collapse = ", ")))
  }
  
  if (is.null(tp)){
    res <- GetMetricsBinary(ground_truth = ground_truth, preds = preds)
    tp <- res$TP
    fp <- res$FP
    tn <- res$TN
    fn <- res$FN
  }
  
  N <- (tp + fp + tn + fn)
  
  res <- list()
  
  if ("accuracy" %in% metrics){
    res[["accuracy"]] <- (tp + tn) / N
  }
  if ("mcc" %in% metrics){ # https://en.wikipedia.org/wiki/Phi_coefficient
    S <- (tp + fn) / N
    P <- (tp + fp) / N
    nominator <-  ((tp/N) - (S * P))
    denominator <- sqrt((S * P) * (1 - S) * (1 - P))
    if (denominator == 0){
      res[["mcc"]] <- NA
    } else {
      res[["mcc"]] <- nominator / denominator
    }
    
  }
  if ("F1" %in% metrics){
    res[["F1"]] <- (2*tp) / ((2*tp) + fp + fn)
  }
  if ("kappa"%in% metrics){
    res[["kappa"]] <- (tp + tn) / N
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
  
  N <- tp + fp + tn + fn
  S <- (tp + fn) / N
  P <- (tp + fp) / N
  
  res <- ((tp/N) - (S * P)) / sqrt((S * P) * (1 - S) * (1 - P))
  
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


