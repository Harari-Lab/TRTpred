# Cross validation

suppressMessages(require(Seurat))
suppressMessages(require(splitTools))
suppressMessages(require(foreach))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))

DEFAULT.LR.HYPERPARAMS <- list("alpha" = 0, "lambda" = 0)
DEFAULT.SIGNATURE.HYPERPARAMS <- list("signature.lengths" = c(20), 
                                      "signature.side" = c("both"),
                                      "signature.rm.regex" = "none",
                                      "signature.methods" = c("AUCell"))
DA.METHODS <- c("wilcox")
P.ADJUST.METHODS <- c("fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "none")
LIMMA.METHODS <- c("limma_voom", "limma_trend")
EDGER.METHODS <- c("edgeR_LRT", "edgeR_QFL")
DESEQ.METHODS <- c("DESeq2_Wald", "DESeq2_LRT")
SEURAT.METHODS <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR")
DEA.METHODS <- c(SEURAT.METHODS, DESEQ.METHODS, LIMMA.METHODS, EDGER.METHODS)
FEATURE.TRANS.METHODS <- c("pca")
EVALUATION.METRICS <- c("accuracy", "mcc", "F1", "kappa")

#' Create Cross-validation folds
#' 
#' Function to create cross validation folds using the `splitTools::create_folds()` function
#' 
#' @param data.input data.frame; input data
#' @param data.output data.frame; output data
#' @param column.split vector of character or just one string; 
#' Columns in `data.output` that will be used to create "balanced" splitted data-sets
#' @param k.outer numerical: The number of outer folds
#' Default = 5
#' @param k.inner numerical: The number of inner folds. If 0, no inner fold is computed
#' Default = 0
#' @param seed numerical: The seed
#' 
#' @return list of list of indexes creating the folds (see ? splitTools::create_folds())
#' 
#' @export
CreateCVFolds <- function(data.input, data.output, column.split, k.outer = 5, k.inner = 0, seed = 1234){
  
  column.split <- c(column.split)
  if (length(column.split) > 1){
    data.output$split_group <- apply(data.output[, column.split], 1, paste, collapse = "_")
  } else {
    data.output$split_group <- data.output[[column.split]]
  }
  
  # Create the outer folds:
  folds.outer <- splitTools::create_folds(y = data.output$split_group, k = k.outer, m_rep = 1, seed = seed)
  
  # Create the container of the inner folds for each outer fold
  folds.outer.inner <- list()
  
  # Create container for summarizing the size of every fold
  fold.size.df <- data.frame()
  
  for (fold.out.name in names(folds.outer)){
    # cat(paste0("fold.out.name: ", fold.out.name, "\n"))
    
    # Get the indexes of outer fold
    fold.out <- folds.outer[[fold.out.name]]
    
    # Get the outer training and testing data
    data.output.train.outer <- data.output[fold.out,]
    data.output.test.outer  <- data.output[-fold.out,]
    data.input.train.outer <- data.input[fold.out,]
    data.input.test.outer  <- data.input[-fold.out,]
    
    # Get the distribution of labels:
    table.train.outer <- table(data.output.train.outer$split_group)
    table.test.outer <- table(data.output.test.outer$split_group)
    
    # Create the inner folds
    if (k.inner > 0){
      # Create the inner folds:
      fold.inner <- splitTools::create_folds(y = data.output.train.outer$split_group, k = k.inner, m_rep = 1, seed = seed)
      
      # Create empty container: 
      folds.outer.inner[[fold.out.name]] <- list()
      for (fold.in.name in names(fold.inner)) {
        # cat(paste0("\tfold.in.name: ", fold.in.name, "\n"))
        
        # Get the indexes of inner fold
        fold.in <- fold.inner[[fold.in.name]]
        
        # Save results in folds.outer.inner
        folds.outer.inner[[fold.out.name]][[fold.in.name]] <- fold.in
        
        # Get the inner training and testing data
        data.output.train.inner <- data.output.train.outer[fold.in,]
        data.output.test.inner  <- data.output.train.outer[-fold.in,]
        data.input.train.inner <- data.input.train.outer[fold.in,]
        data.input.test.inner  <- data.input.train.outer[-fold.in,]
        
        # Get the distribution of labels:
        table.train.inner <- table(data.output.train.inner$split_group)
        table.test.inner <- table(data.output.test.inner$split_group)
        
        fold.size.df.line <- list(
          "outer" = fold.out.name,
          "n_train.outer" = paste(names(table.train.outer), as.numeric(table.train.outer), sep = ": ", collapse = ", "),
          "n_test.outer" = paste(names(table.test.outer), as.numeric(table.test.outer), sep = ": ", collapse = ", "),
          "inner" = fold.in.name,
          "n_train.inner" = paste(names(table.train.inner), as.numeric(table.train.inner), sep = ": ", collapse = ", "),
          "n_test.inner" = paste(names(table.test.inner), as.numeric(table.test.inner), sep = ": ", collapse = ", ")
        )
        
        fold.size.df <- rbind(fold.size.df, fold.size.df.line)
      }
    } else {
      fold.size.df.line <- list(
        "outer" = fold.out.name,
        "n_train.outer" = paste(names(table.train.outer), as.numeric(table.train.outer), sep = ": ", collapse = ", "),
        "n_test.outer" = paste(names(table.test.outer), as.numeric(table.test.outer), sep = ": ", collapse = ", "),
        "inner" = "",
        "n_train.inner" = 0,
        "n_test.inner" = 0
      )
      
      fold.size.df <- rbind(fold.size.df, fold.size.df.line)
    }
  }
  
  return(list("folds.outer" = folds.outer, 
              "folds.outer.inner" = folds.outer.inner, 
              "size.df" = fold.size.df))
}


#' Nested Cross Validataion
#' 
#' Home Nested Cross validation function
#' 
#' @param x data.frame or matrix; The input matrix corresponding to the feature 
#' space only. additional covariates are found in `y`. Rows = samples, 
#' Columns = features
#' @param y data.frame or matrix; The output matrix. 
#' Columns = c(y.label, y.covariates, y.sample)
#' @param folds.outer list; outer-folds definition list. See CreateCVFolds()
#' @param folds.outer.inner list; outer/inner folds definition list. See CreateCVFolds()
#' @param y.label character; The column of `y` corresponding to the binary output
#' @param y.covariates character vector; The columns of `y` corresponding to the covariates
#' @param y.sample character; The column of `y` corresponding to the samples for when method = "signature"
#' @param method character; The prediction method. Possible values are: "LR": 
#' Logistic regression or "signature" signature-score approach. 
#' @param hyperparams list: The hyper-parameter list. 
#' if method = "LR": Elements of list are 
#'   - "alpha" (L2-regularization): Vector of alpha values (0-1)
#'   - "lambda" (L1-regularization): Vector of lambda values (0-1)
#' if method = "signature": Elements of list are 
#'   - "signature.lengths": Vector of signature lengths
#'   - "signature.side": Vector of siganture sides ("up", "down" or "both")
#'   - "signature.rm.regex": Vector of rm.regex (regex to remove genes from signature)
#'   - "signature.methods": Vector of signature score methods
#' @param data.trans.method character: The data transforamtion method. 
#' Possible values are "pca" 
#' Default = "none" i.e. no transformation. 
#' @param pca.explained.var.threshold numeric: If not NULL, apply this threshold 
#' to select PCs explaining more variance than the threhsold
#' @param DA.method character; The discriminant analysis (DA) method. 
#' Possible values are "wilcox"
#' Default = "none" i.e. no discriminant analysis
#' @param DA.p.adjust.method character; The p-value adjust method for the 
#' discriminant analysis. 
#' Default = "fdr"
#' @param DA.p.val.threshold numeric; The p-val significant threshold. 
#' @param LR.prob.threhsold numeric; The probability threshold for the logitic regression
#' @param DEA.method character; The Differential-expression-analylsis method (see RunDEA())
#' default = "wilcox" from the Seurat::FindMarker() function
#' 
#' @return data.frame summarizing the Nested-cross-validation
#' 
#' @export
NestedCrossValidation <- function(x, y, 
                                  folds.outer, folds.outer.inner, 
                                  y.label, y.covariates = NULL, y.sample = NULL,
                                  method = c("LR", "signature"),
                                  hyperparams = NULL,
                                  data.trans.method = c("none", FEATURE.TRANS.METHODS), 
                                  pca.explained.var.threshold = NULL,
                                  DA.method = c("none", DA.METHODS),
                                  DA.p.adjust.method = P.ADJUST.METHODS,
                                  DA.p.val.threshold = 0.05,
                                  LR.prob.threhsold = 0.5,
                                  DEA.method = DEA.METHODS
                                  ){
  
  
  method <- match.arg(method) # default = "LR"
  data.trans.method <- match.arg(data.trans.method) # default = "none"
  DA.method <- match.arg(DA.method) # default = "none"
  DA.p.adjust.method <- match.arg(DA.p.adjust.method) # default = "fdr"
  DEA.method <- match.arg(DEA.method) # default = "wilcox"
  
  names.outer <- names(folds.outer)
  names.inner <- names(folds.outer.inner[[names.outer[1]]])
  
  NCV.res <- 
    foreach::foreach(fold.out.name = names.outer, .combine='rbind') %:% # fold.out.name <- names.outer[1]
    foreach::foreach(fold.in.name = names.inner, .combine='rbind') %dopar% { # fold.in.name <- names.inner[1]
      
      # Get the indexes of outer fold
      fold.out <- folds.outer[[fold.out.name]]
      
      # Get the indexes of the inner fold
      fold.in <- folds.outer.inner[[fold.out.name]][[fold.in.name]]
      
      # Get the outer training and testing data
      x.train.outer <- x[fold.out,, drop=F]
      y.train.outer <- y[fold.out,, drop=F]
      x.test.outer  <- x[-fold.out,, drop=F]
      y.test.outer  <- y[-fold.out,, drop=F]
      
      # Get the inner training and testing data
      x.train.inner <- x.train.outer[fold.in,, drop=F]
      y.train.inner <- y.train.outer[fold.in,, drop=F]
      x.test.inner  <- x.test.outer[-fold.in,, drop=F]
      y.test.inner  <- y.test.outer[-fold.in,, drop=F]
      
      
      if (params.feature.trans$method != "none"){
        data.transformed <- FeatureTrans(
          x.train = x.train.inner, 
          x.test = x.test.inner, 
          method = data.trans.method,
          explained.var.threshold = pca.explained.var.threshold)
        
        x.train.inner <- data.transformed$x.train
        x.test.inner <- data.transformed$x.test
        rm(data.transformed)
      }
      
      # Disciminant analysis
      if (DA.method != "none"){
        discr.res <- DiscrAnalysisBinary(
          x = x.train.inner, 
          y = y.train.inner[, y.label], 
          method = DA.method, 
          p.adjust.method = DA.p.adjust.method, 
          p.val.threshold = DA.p.val.threshold
          )
        sig.features <- rownames(subset(discr.res, Significant))
        
        if (length(sig.features) > 0){
          x.train.inner <- x.train.inner[, sig.features, drop = F]
          x.test.inner <- x.test.inner[, sig.features, drop = F]
        } else {
          warning("Discriminant Analysis: No significant features were identified. We continue without Discriminant Analysis")
        }
      }
      
      if (method == "LR"){
        # Prepare data
        # For LR, the input (x) and output (y) are in the same data.frame
        design.str <- paste(y.label, "~", paste(c(colnames(x.train.inner), y.covariates), collapse = " + "), sep = " ")
        data.train.inner <- cbind(x.train.inner, y.train.inner)
        data.test.inner <- cbind(x.test.inner, y.test.inner)
        
        # Run the training
        res.model <- LRCrossValidation(
          data.train = data.train.inner, 
          data.test = data.test.inner,
          design.str = design.str, 
          hyperparams = hyperparams, 
          model.prob.threshold = LR.prob.threhsold)
        
      } else if (method == "signature"){
        # Create Seurat object
        data.train.inner <- Seurat::CreateSeuratObject(counts = t(x.train.inner), assay = "RNA", project = "NCV", meta.data = y.train.inner)
        data.test.inner <- Seurat::CreateSeuratObject(counts = t(x.test.inner), assay = "RNA", project = "NCV", meta.data = y.test.inner)
        
        # Get the results of the model for all hyperparmeters
        res.model <- SignatureCrossValidation(
          data.train = data.train.inner, 
          data.test = data.test.inner, 
          y.label = y.label, 
          y.sample = y.sample, 
          y.covariates = y.covariates, 
          DEA.method = DEA.method,
          signature.lengths = hyperparams$signature.lengths, 
          signature.sides = hyperparams$signature.side, 
          signature.rm.regex = hyperparams$signature.rm.regex, 
          signature.methods = hyperparams$signature.methods, 
          assay = "RNA", 
          slot = "counts")
      }
      
      
      res.model$outer <- fold.out.name
      res.model$inner = fold.in.name
      
      # Get number of summary of the dimentions to store
      table.train.outer <- table(y.train.outer[, y.label])
      table.test.outer <- table(y.test.outer[, y.label])
      
      table.train.inner <- table(y.train.inner[, y.label])
      table.test.inner <- table(y.test.inner[, y.label])
      
      dimension.list <- list(
        "n.train.outer"=nrow(x.train.outer),
        "p.train.outer"=ncol(x.train.outer),
        "table.train.outer"=paste(names(table.train.outer), as.numeric(table.train.outer), sep = ": ", collapse = ", "),
        "n.test.outer"=nrow(x.test.outer),
        "p.test.outer"=ncol(x.test.outer),
        "table.test.outer"=paste(names(table.test.outer), as.numeric(table.test.outer), sep = ": ", collapse = ", "),
        "n.train.inner"=nrow(x.train.inner),
        "p.train.inner"=ncol(x.train.inner),
        "table.train.inner"=paste(names(table.train.inner), as.numeric(table.train.inner), sep = ": ", collapse = ", "),
        "n.test.inner"=nrow(x.test.inner),
        "p.test.inner"=ncol(x.test.inner),
        "table.test.inner"=paste(names(table.test.inner), as.numeric(table.test.inner), sep = ": ", collapse = ", ")
      )
      for (names_ in names(dimension.list)){
        res.model[[names_]] <- dimension.list[[names_]]
      }
      
    res.model
    }
  
  return(NCV.res)
}

#'  Cross Validataion
#' 
#' Home Cross validation function.
#' This function can be used as a stand-alone Cross-validation function by providing the 
#' full `hyperparams` parameter. 
#' This function can also serve as the outer fold (second part) of a Nested-
#' cross-validation by providing the `hyperparams.best.per.fold` parameter. 
#' 
#' @param x data.frame or matrix; The input matrix corresponding to the feature 
#' space only. additional covariates are found in `y`. Rows = samples, 
#' Columns = features
#' @param y data.frame or matrix; The output matrix. 
#' Columns = c(y.label, y.covariates, y.sample)
#' @param folds.outer list; outer-folds definition list. See CreateCVFolds()
#' @param y.label character; The column of `y` corresponding to the binary output
#' @param y.covariates character vector; The columns of `y` corresponding to the covariates
#' @param y.sample character; The column of `y` corresponding to the samples for when method = "signature"
#' @param method character; The prediction method. Possible values are: "LR": 
#' Logistic regression or "signature" signature-score approach. 
#' @param hyperparams list: The hyper-parameter list. 
#' If not NULL: Stand-alone cross-validation. If NULL, second part of NCV.
#' if method = "LR": Elements of list are 
#'   - "alpha" (L2-regularization): Vector of alpha values (0-1)
#'   - "lambda" (L1-regularization): Vector of lambda values (0-1)
#' if method = "signature": Elements of list are 
#'   - "signature.lengths": Vector of signature lengths
#'   - "signature.side": Vector of siganture sides ("up", "down" or "both")
#'   - "signature.rm.regex": Vector of rm.regex (regex to remove genes from signature)
#'   - "signature.methods": Vector of signature score method
#' @param hyperparams.best.per.fold data.frame; The best partameters for each 
#' fold found by Nested-cross-validation. This parameter is the output of 
#' GetBestOuterHyperparams()
#' @param data.trans.method character: The data transforamtion method. 
#' Possible values are "pca" 
#' Default = "none" i.e. no transformation. 
#' @param pca.explained.var.threshold numeric: If not NULL, apply this threshold 
#' to select PCs explaining more variance than the threhsold
#' @param DA.method character; The discriminant analysis (DA) method. 
#' Possible values are "wilcox"
#' Default = "none" i.e. no discriminant analysis
#' @param DA.p.adjust.method character; The p-value adjust method for the 
#' discriminant analysis. 
#' Default = "fdr"
#' @param DA.p.val.threshold numeric; The p-val significant threshold. 
#' @param LR.prob.threhsold numeric; The probability threshold for the logitic regression
#' @param DEA.method character; The Differential-expression-analylsis method (see RunDEA())
#' default = "wilcox" from the Seurat::FindMarker() function
#' 
#' @return data.frame summarizing the Nested-cross-validation
#' 
#' @export
CrossValidation <- function(x, y, 
                            folds.outer, 
                            y.label, y.covariates = NULL, y.sample = NULL,
                            method = c("LR", "signature"),
                            hyperparams = NULL,
                            hyperparams.best.per.fold = NULL,
                            data.trans.method = c("none", FEATURE.TRANS.METHODS), 
                            pca.explained.var.threshold = NULL,
                            DA.method = c("none", DA.METHODS),
                            DA.p.adjust.method = P.ADJUST.METHODS,
                            DA.p.val.threshold = 0.05,
                            LR.prob.threhsold = 0.5,
                            DEA.method = DEA.METHODS){
  
  # Setting parameters: 
  method <- match.arg(method) # default = "LR"
  data.trans.method <- match.arg(data.trans.method) # default = "none"
  DA.method <- match.arg(DA.method) # default = "none"
  DA.p.adjust.method <- match.arg(DA.p.adjust.method) # default = "fdr"
  DEA.method <- match.arg(DEA.method) # default = "wilcox"
  
  # Get the names of the outer folds: 
  names.outer <- names(folds.outer)
  
  # Set the hyper-parameter list
  # This allows to call this function with hyperparams.best.per.fold or with hyperparams
  hyperparams.list <- list()
  if (!is.null(hyperparams.best.per.fold)){
    for (fold.out.name in names.outer){
      list.tmp <- list()
      for (hp in colnames(hyperparams.best.per.fold)){
        list.tmp[[hp]] <- c(hyperparams.best.per.fold[fold.out.name, hp])
      }
      hyperparams.list[[fold.out.name]] <- list.tmp
    }
    
    if (!is.null(hyperparams)){
      warning("Cannot have both `hyperparams` and `hyperparams.best.per.fold` defined. Continue with hyperparams.best.per.fold")
    }
  } else {
    if (is.null(hyperparams)){
      message("Both `hyperparams` and `hyperparams.best.per.fold` NOT defined.")
      message("Use deafault method hyperparameters")
      if (method == "LR"){
        hyperparams <- DEFAULT.LR.HYPERPARAMS
      } else if (method == "signature"){
        hyperparams <- DEFAULT.SIGNATURE.HYPERPARAMS
      }
      str.output <- ""
      for(i in names(hyperparams)){str.output <- paste(str.output, i, ":", paste(hyperparams[[i]], collapse = ", "), "\n")}
      message(str.output)
    } else {
      for (fold.out.name in names.outer){
        hyperparams.list[[fold.out.name]] <- hyperparams
      }
    }
  }
  
  CV.res <- 
    foreach::foreach(fold.out.name = names.outer, .combine='rbind') %dopar% { # fold.out.name <- names.outer[1]
      
      # Get the indexes of outer fold
      fold.out <- folds.outer[[fold.out.name]]
      
      # Get the outer training and testing data
      x.train.outer <- x[fold.out,, drop=F]
      y.train.outer <- y[fold.out,, drop=F]
      x.test.outer  <- x[-fold.out,, drop=F]
      y.test.outer  <- y[-fold.out,, drop=F]
      
      hyperparams_ <- hyperparams.list[[fold.out.name]]
      
      if (params.feature.trans$method != "none"){
        data.transformed <- FeatureTrans(
          x.train = x.train.outer, 
          x.test = x.test.outer, 
          method = data.trans.method,
          explained.var.threshold = pca.explained.var.threshold)
        
        x.train.outer <- data.transformed$x.train
        x.test.outer <- data.transformed$x.test
        rm(data.transformed)
      }
      
      # Disciminant analysis
      if (DA.method != "none"){
        discr.res <- DiscrAnalysisBinary(
          x = x.train.outer, 
          y = y.train.outer[, y.label], 
          method = DA.method, 
          p.adjust.method = DA.p.adjust.method, 
          p.val.threshold = DA.p.val.threshold
        )
        sig.features <- rownames(subset(discr.res, Significant))
        
        if (length(sig.features) > 0){
          x.train.outer <- x.train.outer[, sig.features, drop = F]
          x.test.outer <- x.test.outer[, sig.features, drop = F]
        } else {
          warning("Discriminant Analysis: No significant features were identified. We continue without Discriminant Analysis")
        }
      }
      
      if (method == "LR"){
        
        # Prepare data
        # For LR, the input (x) and output (y) are in the same data.frame
        design.str <- paste(y.label, "~", paste(c(colnames(x.train.outer), y.covariates), collapse = " + "), sep = " ")
        
        data.train.outer <- cbind(x.train.outer, y.train.outer)
        data.test.outer <- cbind(x.test.outer, y.test.outer)
        
        # Run the training
        res.model <- LRCrossValidation(
          data.train = data.train.outer, 
          data.test = data.test.outer,
          design.str = design.str, 
          hyperparams = hyperparams_, 
          model.prob.threshold = LR.prob.threhsold)
        
      } else if (method == "signature"){
        # Create Seurat object
        data.train.outer <- Seurat::CreateSeuratObject(counts = t(x.train.outer), assay = "RNA", project = "NCV", meta.data = y.train.outer)
        data.test.outer <- Seurat::CreateSeuratObject(counts = t(x.test.outer), assay = "RNA", project = "NCV", meta.data = y.test.outer)
        
        # Get the results of the model for all hyperparmeters
        res.model <- SignatureCrossValidation(
          data.train = data.train.outer, 
          data.test = data.test.outer, 
          y.label = y.label, 
          y.sample = y.sample, 
          y.covariates = y.covariates, 
          DEA.method = DEA.method,
          signature.lengths = hyperparams_$signature.lengths, 
          signature.sides = hyperparams_$signature.side, 
          signature.rm.regex = hyperparams_$signature.rm.regex, 
          signature.methods = hyperparams_$signature.methods, 
          assay = "RNA", 
          slot = "counts")
        
      }
      
      res.model$outer <- fold.out.name

      # Get number of summary of the dimentions to store
      table.train.outer <- table(y.train.outer[, y.label])
      table.test.outer <- table(y.test.outer[, y.label])
      
      dimension.list <- list(
        "n.train.outer"=nrow(x.train.outer),
        "p.train.outer"=ncol(x.train.outer),
        "table.train.outer"=paste(names(table.train.outer), as.numeric(table.train.outer), sep = ": ", collapse = ", "),
        "n.test.outer"=nrow(x.test.outer),
        "p.test.outer"=ncol(x.test.outer),
        "table.test.outer"=paste(names(table.test.outer), as.numeric(table.test.outer), sep = ": ", collapse = ", ")
      )
      for (names_ in names(dimension.list)){
        res.model[[names_]] <- dimension.list[[names_]]
      }
      
      res.model
    }
  
  return(CV.res)
}


#' Get the best hyper-parameters for Nested cross validation
#' 
#' Function to get the best hyper-parameters for the nested cross validation. 
#' It can compute all hyperparameters in a list. 
#' It can simultaneously get the best hyper-parameters for mcc, accuracy, F1 and kappa. 
#' 
#' @param NCV.res data.frame; Nested-Cross-Validation results as returned by NestedCrossValidation()
#' @param method characters; The method of the origin model. possible values are  
#' "LR", "signature"
#' @param outer.col character; The name of the outer-fold column in outer.col
#' 
#' @return data.frame; Data frame summarizing the best hyper-parameters for each outerfold
#' 
#' @export
GetBestOuterHyperparams <- function(NCV.res, 
                                    method = c("LR", "signature"),
                                    outer.col = "outer", 
                                    selection.metric = EVALUATION.METRICS){
  
  selection.metric <- match.arg(selection.metric)
  method <- match.arg(method)
  
  if (method == "LR"){
    hyperparms.col = names(DEFAULT.LR.HYPERPARAMS)
  } else if (method == "signature"){
    hyperparms.col = names(DEFAULT.SIGNATURE.HYPERPARAMS)
  }
  
  if (!(paste0("test.", selection.metric) %in% colnames(NCV.res))){
    stop("GetBestOuterHyperparams: The metric '", selection.metric, "' not in NCV.res")
  }
  
  NCV.res$accaury_ <- NCV.res[[paste0("test.", selection.metric)]]
  
  hyperparams.best <- data.frame(row.names = unique(NCV.res[[outer.col]]))
  for (hp in hyperparms.col){
    if (hp %in% colnames(NCV.res)){
      NCV.res$hp <- NCV.res[[hp]]
      hyperparams.best.tmp <- NCV.res %>%
        group_by_at(outer.col) %>%
        summarise(
          value = hp[which(accaury_ == max(accaury_, na.rm = T))[1]]
        )
      hyperparams.best.tmp <- data.frame(hyperparams.best.tmp)
      rownames(hyperparams.best.tmp) <- hyperparams.best.tmp$outer
      hyperparams.best.tmp$outer <- NULL
      colnames(hyperparams.best.tmp) <- hp
      hyperparams.best <- cbind(hyperparams.best, hyperparams.best.tmp)
    } else {
      warning(paste0("In GetBestOuterHyperparams: Missing hyper-parameter ", hp, 
                     ". Skip this hyper-parameter"))
    }
  }
  
  return(hyperparams.best)
}
