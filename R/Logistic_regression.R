# Logistic Regression

suppressMessages(library(glmnet))
suppressMessages(library(glmnetUtils))
suppressMessages(library(stats))

EVALUATION.METRICS <- c("accuracy", "mcc", "F1", "kappa")

DEFAULT.LR.HYPERPARAMS <- list("alpha" = 0, "lambda" = 0)


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
                              model.prob.threshold = 0.5, folds.save.folder = NULL, 
                              sample.weights = NULL,
                              save.full.model = F){
  
  # Sanity checks:
  for (name_ in names(DEFAULT.LR.HYPERPARAMS)){
    if (!(name_ %in% names(hyperparams))){
      hyperparams[[name_]] <- DEFAULT.LR.HYPERPARAMS[[name_]]
    }
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
  if (!is.null(folds.save.folder)){
    y.train.pred.df <- data.frame(row.names = rownames(data.train))
    y.test.pred.df <- data.frame(row.names = rownames(data.test))
    beta.coef.df <- data.frame(row.names = x.columns)
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
        
        if (!is.null(folds.save.folder)){
          beta.coef.df.tmp <- data.frame(res.fit$beta)
          colnames(beta.coef.df.tmp) <- hyperparams_key
          beta.coef.df <- cbind(beta.coef.df, beta.coef.df.tmp)
          
          # DO NOT do this --> Way to space heavy!!! 
          if (save.full.model){
            # file name = "res_fit_"<alpha>"_"<lambda>".rds"
            saveRDS(res.fit, file = paste0(folds.save.folder, "res_fit_", hyperparams_key, ".rds"))
          }
        }
        
        # Get training accuracies:
        y.train.prob <- stats::predict(object = res.fit, newdata = data.train, type = "respons")[,1] # newx for simple glmnet
        y.train.pred <- y.train.prob < model.prob.threshold # https://stackoverflow.com/questions/60682714/why-did-my-factor-levels-get-reversed-in-logistic-regression-mnist-data
        
        res.train.metrics <- GetMetricsBinary(ground_truth = y.train, preds = y.train.pred)
        res.train.metrics.acc <- GetAccuracyMetrics(tp = res.train.metrics$TP, 
                                                    fp = res.train.metrics$FP, 
                                                    tn = res.train.metrics$TN, 
                                                    fn = res.train.metrics$FN)
        
        if (!is.null(folds.save.folder)){
          y.train.pred.df.tmp <- data.frame(y.train.pred, row.names = rownames(data.train))
          colnames(y.train.pred.df.tmp) <- hyperparams_key
          y.train.pred.df <- cbind(y.train.pred.df, y.train.pred.df.tmp)
          # saveRDS(y.train.pred.df, file = paste0(folds.save.folder, "y_train_pred_", hyperparams_key, ".rds"))
        }
        
        # Get testing accuracy
        if (!is.null(data.test)){
          y.test.prob  <- stats::predict(object = res.fit, newdata = data.test, type = "response")[,1]
          y.test.pred <- y.test.prob < model.prob.threshold
          
          if (!is.null(folds.save.folder)){
            y.test.pred.df.tmp <- data.frame(y.test.pred, row.names = rownames(data.test))
            colnames(y.test.pred.df.tmp) <- hyperparams_key
            y.test.pred.df <- cbind(y.test.pred.df, y.test.pred.df.tmp)
            # saveRDS(y.test.pred.df, file = paste0(folds.save.folder, "y_test_pred_", hyperparams_key, ".rds"))
          }
          
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
  
  if (!is.null(folds.save.folder)){
    saveRDS(y.train.pred.df, file =  paste0(folds.save.folder, "y_train_pred.rds"))
    if (!is.null(data.test)){
      saveRDS(y.test.pred.df, file = paste0(folds.save.folder, "y_test_pred.rds"))
    }
    saveRDS(beta.coef.df, file = paste0(folds.save.folder, "LR_beta_coef.rds"))
  }
  
  return(res.inner)
}
