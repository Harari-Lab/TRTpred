# Logistic Regression

suppressMessages(library(glmnet))
suppressMessages(library(glmnetUtils))
suppressMessages(library(stats))

EVALUATION.METRICS <- c("accuracy", "mcc", "F1", "kappa")

#' Logitic Regression Cross validation function
#' 
#' The function to perform logistic regression within a cross-validation frame-work
#' 
#' @param data.train data.frame/matrix. The training input data with the dependent and independent variables. 
#' @param data.test data.frame/matrix. The testing input data with the dependent and independent variables
#' Default = NULL i.e. not computed
#' @param design.str character; The design of the model 
#' @param hyperparams list; The hyper-parameter list with two elements:
#' - "alpha" (L2-regularization): Vector of alpha values (0-1)
#' - "lambda" (L1-regularization): Vector of lambda values (0-1)
#' @param model.prob.threshold numerical; The model probability threshold
#' Default = 0.5
#' 
#' @return data.frame; The dataframe summarizing the accuracy of the model for varying hyperparameters
#' 
#' @export
LRCrossValidation <- function(data.train, design.str, hyperparams, data.test = NULL, 
                              model.prob.threshold = 0.5){
  
  # Sanity checks:
  if (!("lambda" %in% names(hyperparams))){
    hyperparams$lambda <- c(0)
  }
  if (!("alpha" %in% names(hyperparams))){
    hyperparams$lambda <- c(0)
  }
  
  y.columns <- gsub(" *~.*", "", design.str)
  y.train <- data.train[, y.columns]
  y.train <- as.character(y.train) == levels(y.train)[1]
  if (!is.null(data.test)){
    y.test <- data.test[, y.columns]
    y.test <- as.character(y.test) == levels(y.test)[1]
  }
  
  
  # Formula doesn't work with "-"
  # So we replace them with a double point ".." beacause it is easy to find and replace later on. 
  colnames(data.train) <- gsub("[-]", "..", colnames(data.train))
  if (!is.null(data.test)){
    colnames(data.test) <- gsub("[-]", "..", colnames(data.test))
  }
  design.str <- gsub("[-]", "..", design.str)
  
  x.columns <- trimws(strsplit(gsub(".*~ *", "", design.str), "[+]")[[1]])
 
  # Create Empty containers
  res.inner <- data.frame()
  
  if (length(x.columns) > 1){
    for (lambda in hyperparams$lambda){ # lambda <- 5
      for(alpha in hyperparams$alpha){# alpha <- 0.8
        # Train the models
        # lambda for L1, alpha for L2
        res.fit <- glmnetUtils::glmnet(
          formula = as.formula(design.str),
          data = data.train,
          family = "binomial", 
          alpha = alpha, 
          lambda = lambda)
        
        # Get training accuracies:
        y.train.prob <- stats::predict(object = res.fit, newdata = data.train, type = "respons")[,1] # newx for simple glmnet
        y.train.pred <- y.train.prob < model.prob.threshold # https://stackoverflow.com/questions/60682714/why-did-my-factor-levels-get-reversed-in-logistic-regression-mnist-data
        
        res.train.metrics <- GetMetricsBinary(ground_truth = y.train, preds = y.train.pred)
        res.train.metrics.acc <- GetAccuracyMetrics(tp = res.train.metrics$TP, 
                                                    fp = res.train.metrics$FP, 
                                                    tn = res.train.metrics$TN, 
                                                    fn = res.train.metrics$FN)
        
        # Get testing accuracy
        if (!is.null(data.test)){
          y.test.prob  <- stats::predict(object = res.fit, newdata = data.test, type = "response")[,1]
          y.test.pred <- y.test.prob < model.prob.threshold
          
          res.test.metrics <- GetMetricsBinary(ground_truth = y.test, preds = y.test.pred)
          res.test.metrics.acc <- GetAccuracyMetrics(tp = res.test.metrics$TP, 
                                                     fp = res.test.metrics$FP, 
                                                     tn = res.test.metrics$TN, 
                                                     fn = res.test.metrics$FN)
        } else {
          res.test.metrics <- list(TP = NA, FP = NA, TN = NA, FN = NA)
          res.test.metrics.acc <- list("accuracy" = NA, "mcc" = NA, "F1" = NA, "kappa" = NA) 
        }
        
        # Save results in cv.df.outer.inner:
        res.info <- list(
          "lambda" = lambda, 
          "alpha" = alpha,
          "model.prob.threshold" = model.prob.threshold)
        
        names(res.train.metrics) <- paste0("train.", names(res.train.metrics))
        names(res.train.metrics.acc) <- paste0("train.", names(res.train.metrics.acc))
        
        names(res.test.metrics) <- paste0("test.", names(res.test.metrics))
        names(res.test.metrics.acc) <- paste0("test.", names(res.test.metrics.acc))
        
        res.info <- c(res.info, res.train.metrics, res.train.metrics.acc, res.test.metrics, res.test.metrics.acc)
        
        res.inner <- rbind(res.inner, res.info)
      } 
    } 
  } else {
    for (lambda in hyperparams$lambda){ # lambda <- 0
      for(alpha in hyperparams$alpha){
        # Save results in cv.df.outer.inner:
        res.info <- list(
          "lambda" = lambda, 
          "alpha" = alpha)
        
        res.train.metrics <- replicate(n = 4, NA)
        names(res.train.metrics) <- paste0("train.", c("TP", "TN", "FP", "FN"))
        res.train.metrics.acc <- replicate(n = length(EVALUATION.METRICS), NA)
        names(res.train.metrics.acc) <- paste0("train.", EVALUATION.METRICS)
        
        res.test.metrics <- replicate(n = 4, NA)
        names(res.test.metrics) <- paste0("test.", c("TP", "TN", "FP", "FN"))
        res.test.metrics.acc <- replicate(n = length(EVALUATION.METRICS), NA)
        names(res.test.metrics.acc) <- paste0("test.", EVALUATION.METRICS)
        
        res.info <- c(res.info, res.train.metrics, res.train.metrics.acc, res.test.metrics, res.test.metrics.acc)
        
        res.inner <- rbind(res.inner, res.info)
      }
    }
  }
  
  return(res.inner)
}
