# Logistic Regression functions
# Author: Rémy Pétremand
# Date: 07.05.2024
# Description: Function to train and test logistic regression models
# Reference: https://doi.org/10.1038/s41587-024-02232-0

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

suppressMessages(library(glmnet))
suppressMessages(library(glmnetUtils))
suppressMessages(library(stats))

# ------------------------------------------------------------------------------
# Global Parameters
# ------------------------------------------------------------------------------

EVALUATION.METRICS <- c("mcc", "accuracy", "F1", "kappa", "auc", "sensitivity", "specificity", "PPV", "NPV")

DEFAULT.LR.HYPERPARAMS <- list("alpha" = 0, "lambda" = 0)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

#' Get Prediction for the LR model
#' 
#' The function to get prediction from a LR model 
#' 
#' @param x data.frame/matrix. The input data with the dependent and independent variables. 
#' @param path.folder character; Folder to retrieve the model
#' @param hyperparameters list; The LR hyperparmeters
#' Elements of list are 
#'   - "alpha" (L2-regularization): Vector of alpha values (0-1)
#'   - "lambda" (L1-regularization): Vector of lambda values (0-1)
#' @param model.file.name character; the file name of the model in the folder 
#' if null retrive it from "res_fit_"<alpha>"_"<lambda>".rds"
#' Default = NULL
#' @param prob.threshold numerical; The LR probability threshold
#' Default = 0.5
#' 
#' @return data.frame; Rows = barcode names in x. Cols = "score" and "pred"
#' 
#' @export
GetLRPrediction <- function(x, path.folder,
                            hyperparameters = DEFAULT.LR.HYPERPARAMS,
                            model.file.name = NULL,
                            prob.threshold = 0.5){
  
  # 0. Check hyperparameters
  for (name_ in names(DEFAULT.LR.HYPERPARAMS)){
    if (!(name_ %in% names(hyperparameters))){
      hyperparameters[[name_]] <- DEFAULT.LR.HYPERPARAMS[[name_]]
    }
  }
  
  # 1. Get the LR model i.e. the glmnet object
  hyperparams_key <- paste(hyperparameters$alpha, hyperparameters$lambda, sep = "_")
  
  # Get the model.file.name
  if (is.null(model.file.name)){
    model.file.name <- paste0("res_fit_", hyperparams_key, ".rds")
  }
  # Get model.res
  model.path <- paste0(path.folder, model.file.name)
  if (file.exists(model.path)){
    res.fit <- readRDS(file = model.path)
  } else {
    stop("GetLRPrediction: Model path does not exist. Please check : ", model.path)
  }
  
  # 2. Get prediction

  colnames(x) <- gsub("[-]", "..", colnames(x))
  mask.missing <- !(rownames(res.fit$beta) %in% colnames(x))
  if (any(mask.missing)){
    colnames.2.add <- rownames(res.fit$beta)[mask.missing]
    
    warning("GetLRPrediction: Missing features in `x` \n", paste(colnames.2.add, collapse = ", "))
    
    data.2.add <- matrix(nrow = nrow(x), ncol = length(colnames.2.add), data = 0)
    data.2.add <- data.frame(data.2.add)
    rownames(data.2.add) <- rownames(x)
    colnames(data.2.add) <- colnames.2.add
    
    x <- cbind(x, data.2.add)
  }

  y.prob  <- stats::predict(object = res.fit, newdata = x, type = "response")[,1]
  y.pred <- y.prob < prob.threshold
  
  # 3. Get return data.frame
  LR.pred <- data.frame("score" = y.prob, "pred" = y.pred, row.names = rownames(x))
  
  return(LR.pred)
}

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
#' @param folds.save.folder character; Folder to save the model information. 
#' If null doesn't save anything. If not null, it will save the following data into the folder: 
#' - res.fit.rds: glmnetUtils::glmnet() results
#' - y.train.pred.rds: y.train.pred data.frame (rows = cell name, col = prediction)
#' - y.test.pred.rds: y.test.pred data.frame (rows = cell name, col = prediction)
#' default = NULL
#' @param sample.weights data.frame (rows = observation, column = weight); The 
#' observation weights data.frame. 
#' Default = NULL i.e. 1 for each observation
#' 
#' @return data.frame; The dataframe summarizing the accuracy of the model for varying hyperparameters
#' 
#' @export
LRCrossValidation <- function(data.train, design.str, hyperparams, data.test = NULL, 
                              model.prob.threshold = NULL, 
                              sample.weights = NULL,
                              folds.save.folder = NULL, 
                              save.model = F, save.pred = F, save.score = F){
  
  # Sanity check on the saving of data
  if (is.null(folds.save.folder)){
    if (any(c(save.model, save.pred, save.score))){
      warning("folds.save.folder is NULL but one or more of `save.model`, `save.pred`, `save.score` is set as TRUE: Nothing will be saved")
      save.model = F;save.pred = F;save.score = F
    }
  } else {
    if (dir.exists(folds.save.folder)){
      if (any(!c(save.model, save.pred, save.score))){
        warning("folds.save.folder exists but nothing will be saved since `save.model`, `save.pred`, `save.score` are set as FALSE")
      }
    }
  }
  
  # Sanity checks:
  for (name_ in names(DEFAULT.LR.HYPERPARAMS)){
    if (!(name_ %in% names(hyperparams))){
      hyperparams[[name_]] <- DEFAULT.LR.HYPERPARAMS[[name_]]
    }
  }
  
  if (is.null(model.prob.threshold)){
    binarize.method <- "greedy"
  } else {
    binarize.method <- "threshold"
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
  mask.start.num <- grepl(pattern = "^[0-9]", x = colnames(data.train))
  if (any(mask.start.num)){colnames(data.train)[mask.start.num] <- paste0("..", colnames(data.train)[mask.start.num])}
  if (!is.null(data.test)){
    colnames(data.test) <- gsub("[-]", "..", colnames(data.test))
    mask.start.num <- grepl(pattern = "^[0-9]", x = colnames(data.test))
    if (any(mask.start.num)){colnames(data.test)[mask.start.num] <- paste0("..", colnames(data.test)[mask.start.num])}
  }
  # rm the "-" in gene names and replace them by ".."
  design.str <- gsub("[-]", "..", design.str)
  # add ".." in front of gene names strarting with a number
  design.str <- gsub(pattern = " [+] ([0-9])", replacement = " + ..\\1", x = design.str)
  
  x.columns <- trimws(strsplit(gsub(".*~ *", "", design.str), "[+]")[[1]])
 
  # Create Empty containers
  res.inner <- data.frame()
  if (save.pred){
    y.train.pred.df <- data.frame(row.names = rownames(data.train))
    y.test.pred.df <- data.frame(row.names = rownames(data.test))
    beta.coef.df <- data.frame(row.names = x.columns)
  }
  if (save.score){
    y.train.score.df <- data.frame(row.names = rownames(data.train))
    y.test.score.df <- data.frame(row.names = rownames(data.test))
  }
  
  if (is.null(sample.weights)){
    weights <- replicate(n = nrow(data.train), expr = 1)
  } else {
    weights <- sample.weights[rownames(data.train),1]
  }
  
  if (length(x.columns) > 1){
    for (lambda in hyperparams$lambda){ # lambda <- 0.5
      for(alpha in hyperparams$alpha){ # alpha <- 0.1
        # <alpha>_<lambda>
        hyperparams_key <- paste(alpha, lambda, sep = "_")
        
        # Train the models
        # lambda for L1, alpha for L2
        res.fit <- glmnetUtils::glmnet(
          formula = as.formula(design.str),
          data = data.train,
          family = "binomial", 
          alpha = alpha, 
          weights = weights,
          lambda = lambda)
        
        # Get training probabilities
        # https://stackoverflow.com/questions/60682714/why-did-my-factor-levels-get-reversed-in-logistic-regression-mnist-data
        y.train.prob <- stats::predict(object = res.fit, newdata = data.train, type = "respons")[,1] # newx for simple glmnet
        y.train.prob <- 1 - y.train.prob
        
        # Get testting probabilities
        if (!is.null(data.test)){
          y.test.prob  <- stats::predict(object = res.fit, newdata = data.test, type = "response")[,1]
          y.test.prob <- 1 - y.test.prob
        } else {
          y.test.prob <- NULL
        }
        
        # Get the binary prediction (greedy or threshold)
        pred.res <- BinarizeScorePrediction(
          x.train = y.train.prob, 
          x.test = y.test.prob,
          method = binarize.method, 
          weights = NULL,
          ground.truth = y.train, 
          x.threshold = model.prob.threshold)
        
        # Get train & test predictions: 
        y.train.pred <- pred.res$y.train.pred
        y.test.pred <- pred.res$y.test.pred
        x.threshold <- pred.res$x.threshold
        
        # Get training accuracy
        res.train.metrics <- GetAccuracyMetrics(ground_truth = y.train, preds = y.train.pred, scores = y.train.prob, metrics = EVALUATION.METRICS)
        
        # Get testing accuracy
        if (!is.null(data.test)){
          res.test.metrics <- GetAccuracyMetrics(ground_truth = y.test, preds = y.test.pred, scores = y.test.prob, metrics = EVALUATION.METRICS)
        } else {
          res.test.metrics <- replicate(n = (4 + length(EVALUATION.METRICS)), NA)
          names(res.test.metrics) <- c("TP", "FP", "TN", "FN", EVALUATION.METRICS)
        }
        
        # Save results in cv.df.outer.inner:
        res.info <- list(
          "lambda" = lambda, 
          "alpha" = alpha,
          "model.prob.threshold" = x.threshold)
        
        names(res.train.metrics) <- paste0("train.", names(res.train.metrics))
        
        names(res.test.metrics) <- paste0("test.", names(res.test.metrics))
        
        res.info <- c(res.info, res.train.metrics, res.test.metrics)
        
        res.inner <- rbind(res.inner, res.info)
        
        
        # Save Model - WARNING VERY SPACE HEAVY
        if (save.model){
          # file name = "res_fit_"<alpha>"_"<lambda>".rds"
          saveRDS(res.fit, file = paste0(folds.save.folder, "res_fit_", hyperparams_key, ".rds"))
        }
        # Store Prediction - Less space heavy
        if (save.pred){
          y.train.pred.df.tmp <- data.frame(y.train.pred, row.names = rownames(data.train))
          colnames(y.train.pred.df.tmp) <- hyperparams_key
          y.train.pred.df <- cbind(y.train.pred.df, y.train.pred.df.tmp)
          
          if (!is.null(data.test)){
            y.test.pred.df.tmp <- data.frame(y.test.pred, row.names = rownames(data.test))
            colnames(y.test.pred.df.tmp) <- hyperparams_key
            y.test.pred.df <- cbind(y.test.pred.df, y.test.pred.df.tmp)
          }
          
          beta.coef.df.tmp <- data.frame(res.fit$beta)
          colnames(beta.coef.df.tmp) <- hyperparams_key
          beta.coef.df <- cbind(beta.coef.df, beta.coef.df.tmp)
        }
        # Store Scores - WARNING Space heavy
        if (save.score){
          y.train.score.df.tmp <- data.frame(y.train.prob, row.names = rownames(data.train))
          colnames(y.train.score.df.tmp) <- hyperparams_key
          y.train.score.df <- cbind(y.train.score.df, y.train.score.df.tmp)
          
          if (!is.null(data.test)){
            y.test.score.df.tmp <- data.frame(y.test.prob, row.names = rownames(data.test))
            colnames(y.test.score.df.tmp) <- hyperparams_key
            y.test.score.df <- cbind(y.test.score.df, y.test.score.df.tmp)
          }
        }
      } 
    } 
  } else {
    for (lambda in hyperparams$lambda){ # lambda <- 0
      for(alpha in hyperparams$alpha){
        # Save results in cv.df.outer.inner:
        res.info <- list(
          "lambda" = lambda, 
          "alpha" = alpha,
          "model.prob.threshold" = 0.5)
        
        res.train.metrics <- replicate(n = (4 + length(EVALUATION.METRICS)), NA)
        names(res.train.metrics) <- paste0("train.", c("TP", "TN", "FP", "FN", EVALUATION.METRICS))
        
        res.test.metrics <- replicate(n = (4 + length(EVALUATION.METRICS)), NA)
        names(res.test.metrics) <- paste0("test.", c("TP", "TN", "FP", "FN", EVALUATION.METRICS))
        
        res.info <- c(res.info, res.train.metrics, res.test.metrics)
        
        res.inner <- rbind(res.inner, res.info)
      }
    }
  }
  
  if (!is.null(folds.save.folder)){
    
  }
  
  # Save Prediction - Less space heavy
  if (save.pred){
    saveRDS(y.train.pred.df, file =  paste0(folds.save.folder, "y_train_pred.rds"))
    if (!is.null(data.test)){
      saveRDS(y.test.pred.df, file = paste0(folds.save.folder, "y_test_pred.rds"))
    }
    saveRDS(beta.coef.df, file = paste0(folds.save.folder, "LR_beta_coef.rds"))
  }
  # Save Scores - WARNING Space heavy
  if (save.score){
    saveRDS(y.train.score.df, file =  paste0(folds.save.folder, "y_train_score.rds"))
    if (!is.null(data.test)){
      saveRDS(y.test.score.df, file = paste0(folds.save.folder, "y_test_score.rds"))
    }
  }
  
  return(res.inner)
}
