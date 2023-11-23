# Cross validation

suppressMessages(require(Seurat))
suppressMessages(require(foreach))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))

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





#' Create Cross-validation folds
#' 
#' Function to create cross validation folds
#' 
#' @param y vector; the vector to create the CV folds from
#' @param k numerical: The number of outer folds
#' Default = 5
#' @param y.leave.1.out vector; If not NULL, vector to create Leave-1-out category CV
#' @param sampling.method character; The sampling method to balance the folds. 
#' if "under", Undersampling is performed to balance the categories in every fold
#' if "upper", Uppersampling is performed to balance the categories in every fold
#' if "none", no balancing sampling is performed. 
#' Default = "under"
#' @param replicates numerical; The number of replicates. 
#' At each replicate, the seed is increased by 1
#' @param seed numerical: The seed Default = .Random.seed[1]
#' 
#' @return list containing the index of each sample belonging to the training and testing set for every fold
#' names of list = names of CV folds. values of list = list with two names 
#' - train: Training sample indexes
#' - test: Testing sample indexes
#' 
#' @export
CreateCVFolds <- function(y, y.leave.1.out = NULL, k = 5, replicates = 1,
                          sampling.method = c("under", "upper", "none"), 
                          seed = .Random.seed[1]){
  
  # Sanity checks:
  sampling.method <- match.arg(sampling.method)
  if (replicates < 1){
    replicates <- 1
  }
  
  # Get the mask of positive y values: 
  y <- as.factor(y)
  mask_y_pos <- y == levels(y)[1]
  
  # Prepare the split group depending upon wheather you want to do a classical 
  # k-fold CV or a leave-one-out based on values in a vector (i.e. y.leave.1.out)
  if (!is.null(y.leave.1.out)){
    fold_start_name <- "FoldL1O"
    
    y.split.group <- as.character(y.leave.1.out)
    y.split.group.classes <- unique(y.leave.1.out)
  } else if (!is.null(k)){
    fold_start_name <- "FoldK"
    
    idx.pos <- which(mask_y_pos)
    idx.neg <- which(!mask_y_pos)
    
    y.split.group <- replicate(n = length(y), NA) 
    
    set.seed(seed)
    idx.pos.k.groups <- sample(1:length(idx.pos)%%k)
    for (k_i in 0:(k-1)){ # k_i = 0
      y.split.group[idx.pos[which(idx.pos.k.groups == k_i)]] <- k_i
    }
    set.seed(seed)
    idx.neg.k.groups <- sample(1:length(idx.neg)%%k)
    for (k_i in 0:(k-1)){ # k_i = 0
      y.split.group[idx.neg[which(idx.neg.k.groups == k_i)]] <- k_i
    }
    y.split.group <- as.character(y.split.group)
    y.split.group.classes <- as.character(0:(k-1))
  } else {
    stop(" IN CreateCVFolds: k or y.leave.1.out must be defined")
  }
  
  
  fold.list <- list() # empty container
  for (rep_ in 1:replicates){ # rep_ <- replicates[1]
    for (y.split.cat in y.split.group.classes){ # y.split.cat <- y.split.group.classes[1]
      fold.outer.name <- paste0(fold_start_name, "_", y.split.cat, "_", rep_)
      
      idx.train.pos <- which(y.split.group != y.split.cat & mask_y_pos)
      idx.train.neg <- which(y.split.group != y.split.cat & !mask_y_pos)
      idx.test.pos <-  which(y.split.group == y.split.cat & mask_y_pos)
      idx.test.neg <-  which(y.split.group == y.split.cat & !mask_y_pos)
      
      if (sampling.method == "under"){
        train.length <- min(c(length(idx.train.pos), length(idx.train.neg)))
        test.length <-  min(c(length(idx.test.pos),  length(idx.test.neg)))
        
        set.seed(seed = seed)
        idx.train<- c(idx.train.pos[sample.int(n = length(idx.train.pos), size = train.length, replace = F)], 
                      idx.train.neg[sample.int(n = length(idx.train.neg), size = train.length, replace = F)])
        set.seed(seed = seed)
        idx.test <- c(idx.test.pos[sample.int(n = length(idx.test.pos), size = test.length, replace = F)], 
                      idx.test.neg[sample.int(n = length(idx.test.neg), size = test.length, replace = F)])
        
      } else if (sampling.method == "upper"){
        train.length <- max(c(length(idx.train.pos), length(idx.train.neg)))
        test.length <-  max(c(length(idx.test.pos),  length(idx.test.neg)))
        
        set.seed(seed = seed)
        idx.train <- c(idx.train.pos[sample.int(n = length(idx.train.pos), size = train.length, replace = T)], 
                       idx.train.neg[sample.int(n = length(idx.train.neg), size = train.length, replace = T)])
        set.seed(seed = seed)
        idx.test <- c(idx.test.pos[sample.int(n = length(idx.test.pos), size = test.length, replace = T)], 
                      idx.test.neg[sample.int(n = length(idx.test.neg), size = test.length, replace = T)])
      } else {
        idx.train <- c(idx.train.pos, idx.train.neg)
        idx.test <- c(idx.test.pos, idx.test.neg)
      }
      fold.list[[fold.outer.name]] <- list("train" = idx.train, "test" = idx.test)
      seed <- seed + 1 # change seed for the next replicate
    } 
  }
  
  return(fold.list)
}


#' Create Nested-Cross-validation folds
#' 
#' Function to create Nested cross validation folds using the CreateCVFolds in this package.
#' See ? CreateCVFolds() for more details
#' 
#' @param y vector; the vector to create the CV folds from
#' @param k.out numerical: The number of outer folds
#' Default = 5
#' @param k.in numerical: The number of inner folds
#' Default = 5
#' @param y.leave.1.out vector; If not NULL, vector to create Leave-1-out category CV
#' @param sampling.method character; The sampling method to balance the folds. 
#' if "under", Undersampling is performed to balance the categories in every fold
#' if "upper", Uppersampling is performed to balance the categories in every fold
#' if "none", no balancing sampling is performed. 
#' Default = "under"
#' @param replicates numerical; The number of replicates. 
#' At each replicate, the seed is increased by 1
#' @param seed numerical: The seed Default = .Random.seed[1]
#' 
#' @return list of list. first level = outer-folds, second level = inner folds. 
#' The content of the inner-fold list is the same as in CreateCVFolds()
#' 
#' @export
CreateNestedCVFolds <- function(y, y.leave.1.out = NULL, k.out = 5, k.in = 5,
                                replicates = 1,
                                sampling.method = c("under", "upper", "none"), 
                                seed = .Random.seed[1]){
  
  # Sanity checks:
  sampling.method <- match.arg(sampling.method)
  
  if (replicates < 1){
    replicates <- 1
  }
  
  # Create container for summarizing the size of every fold
  fold.size.df <- data.frame()
  
  # Create the outer folds:
  folds.outer <- CreateCVFolds(
    y = y, 
    k = k.out,
    y.leave.1.out = y.leave.1.out, 
    replicates = replicates, 
    sampling.method = sampling.method,
    seed = seed)
  
  # Create the container of the inner folds for each outer fold
  folds.outer.inner <- list()
  for (fold.out.name in names(folds.outer)){ # fold.out.name <- names(folds.outer)[1]
    
    # Get the indexes of outer fold
    fold.out <- folds.outer[[fold.out.name]]
    
    # Get the outer training and testing data
    y.train.outer <- y[fold.out$train]
    y.test.outer  <- y[fold.out$test]
    if (!is.null(y.leave.1.out)){
      y.leave.1.out.train.outer <- y.leave.1.out[fold.out$train]
      y.leave.1.out.test.outer  <- y.leave.1.out[fold.out$test]
    } else {
      y.leave.1.out.train.outer <- NULL
      y.leave.1.out.test.outer  <- NULL
    }
    
    # Get the distribution of labels:
    table.train.outer <- table(y.train.outer)
    table.test.outer <- table(y.test.outer)
    
    # Create the inner folds
    # if (k.in > 0){
    folds.inner <- CreateCVFolds(
      y = y.train.outer, 
      k = k.in,
      y.leave.1.out = y.leave.1.out.train.outer, 
      replicates = replicates, 
      sampling.method = sampling.method,
      seed = seed)
    
    # Create empty container: 
    folds.outer.inner[[fold.out.name]] <- folds.inner
    for (fold.in.name in names(folds.inner)) {
      
      # Get the indexes of inner fold
      fold.in <- folds.inner[[fold.in.name]]
      
      # Get the inner training and testing data
      y.train.inner <- y.train.outer[fold.in$train]
      y.test.inner  <- y.train.outer[fold.in$test]
      if (!is.null(y.leave.1.out)){
        y.leave.1.out.train.inner <- y.leave.1.out.train.outer[fold.in$train]
        y.leave.1.out.test.inner <- y.leave.1.out.train.outer[fold.in$test]
      }
      
      # Get the distribution of labels:
      table.train.inner <- table(y.train.inner)
      table.test.inner <- table(y.test.inner)
      
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
#' @param rm.corr.features logical: Do we remove the mutually correlated feature 
#' while keeping the ones that best correlate with the outcome? 
#' default = F
#' @param DA.method character; The discriminant analysis (DA) method. 
#' Possible values are "wilcox"
#' Default = "none" i.e. no discriminant analysis
#' @param DA.p.adjust.method character; The p-value adjust method for the 
#' discriminant analysis. 
#' Default = "fdr"
#' @param DA.event.per.variable, numerical; The number of event per variable threshold
#' A common EPV is 10 meaning that we need 10 event per variable. Here we use it as 
#' a way to select a fewer amount of variables for a fixed amount of samples. 
#' Default = NULL (no selection on EPV)
#' @param DA.p.val.threshold numeric; The p-val significant threshold. 
#' @param LR.prob.threhsold numeric; The probability threshold for the logitic regression
#' @param signature.col.aggregate character; Column in y to combine cells into pseudobulk. 
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
#' @param folds.save.folder character; Folder to save the fold model information. 
#' If null doesn't save anything. if not null, a folder called "NCV_", fold.out.name, "_", fold.in.name
#' will be created to store the model information in. 
#' default = NULL
#' 
#' @return data.frame summarizing the Nested-cross-validation
.NestedCrossValidationLoop <- function(x, y, 
                                  folds.outer, folds.outer.inner, 
                                  y.label, y.covariates = NULL, y.sample = NULL,
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
  
  method <- match.arg(method) # default = "LR"
  data.trans.method <- match.arg(data.trans.method) # default = "none"
  DA.method <- match.arg(DA.method) # default = "none"
  DA.p.adjust.method <- match.arg(DA.p.adjust.method) # default = "fdr"
  DEA.method <- match.arg(DEA.method) # default = "wilcox"
  
  if (!is.null(DEA.data)){
    if (any(dim(DEA.data) != dim(x))){
      warning("DEA.data not the same format as x. Lets try to harmonize the data")
      tryCatch(expr = {
        DEA.data <- DEA.data[rownames(x), colnames(x)]
      }, error = function(cond){
        stop("DEA.data should have the same columns and rows as x")
      })
    }
  }
  
  names.outer <- names(folds.outer)
  
  # fold.out.name <- "FoldK_4_1"; fold.in.name <- "FoldK_3_1"
  # fold.out.name <- "FoldL1O_OvCa210_1"; fold.in.name <- "FoldL1O_CrCp7_1"
  
  NCV.res <- 
    foreach::foreach(fold.out.name = names.outer, .combine='rbind') %:% # fold.out.name <- names.outer[1]
    foreach::foreach(fold.in.name = names(folds.outer.inner[[fold.out.name]]), .combine='rbind') %dopar% { # fold.in.name <- names(folds.outer.inner[[fold.out.name]])[1]
      message(paste0("fold.out.name: ", fold.out.name, "\tfold.in.name: ", fold.in.name))
      
      # Get the indexes of outer fold
      fold.out <- folds.outer[[fold.out.name]]
      
      # Get the indexes of the inner fold
      fold.in <- folds.outer.inner[[fold.out.name]][[fold.in.name]]
      
      # Get the outer training and testing data
      x.train.outer <- x[fold.out$train,, drop=F]
      y.train.outer <- y[fold.out$train,, drop=F]
      x.test.outer  <- x[fold.out$test,, drop=F] # we do not need it here
      y.test.outer  <- y[fold.out$test,, drop=F] # we do not need it here
      
      # Get the inner training and testing data
      x.train.inner <- x.train.outer[fold.in$train,, drop=F]
      y.train.inner <- y.train.outer[fold.in$train,, drop=F]
      x.test.inner  <- x.train.outer[fold.in$test,, drop=F]
      y.test.inner  <- y.train.outer[fold.in$test,, drop=F]
      
      # Get the inner and outer DEA data
      DEA.data.train.outer <- NULL
      DEA.data.train.inner <- NULL
      if (!is.null(DEA.data)){
        DEA.data.train.outer <- DEA.data[fold.out$train,, drop=F]
        DEA.data.train.inner <- DEA.data.train.outer[fold.in$train,, drop=F]
      }
      
      if (!is.null(folds.save.folder)){
        folds.save.folder.tmp <- paste0(folds.save.folder, "NCV_", fold.out.name, "_", fold.in.name, "/")
        if (!dir.exists(folds.save.folder.tmp)){
          system(command = paste("mkdir", "-p", folds.save.folder.tmp))
        }
      } else {
        folds.save.folder.tmp <- NULL
      }
      
      res.model <- TrainModel(
        x.train = x.train.inner,
        y.train = y.train.inner,
        x.test = x.test.inner,
        y.test = y.test.inner,
        y.label = y.label,
        y.covariates = y.covariates,
        y.sample = y.sample,
        method = method,
        hyperparams = hyperparams,
        data.trans.method = data.trans.method,
        pca.explained.var.threshold = pca.explained.var.threshold,
        rm.corr.features = rm.corr.features,
        DA.method = DA.method,
        DA.p.adjust.method = DA.p.adjust.method,
        DA.p.val.threshold = DA.p.val.threshold,
        DA.event.per.variable = DA.event.per.variable,
        LR.prob.threhsold = LR.prob.threhsold,
        signature.col.aggregate = signature.col.aggregate,
        DEA.method = DEA.method,
        DEA.data = DEA.data.train.inner,
        sample.weights = sample.weights,
        folds.save.folder = folds.save.folder.tmp, 
        save.fold.model = save.fold.model,
        save.fold.pred = save.fold.pred,
        save.fold.score = save.fold.score)
     
      # Add to the results additional feature dimension and properties
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
      
      message(paste0("Start : fold.out.name: ", fold.out.name, "\tfold.in.name: ", fold.in.name))
      
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
#' @param CV.folds list; Cross-Validation folds definition list. See CreateCVFolds()
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
#' @param signature.col.aggregate character; Column in y to combine cells into pseudobulk. 
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
#' @param folds.save.folder character; Folder to save the fold model information. 
#' If null doesn't save anything. if not null, a folder called "CV_", fold.out.name
#' will be created to store the model information in. 
#' default = NULL
#' 
#' @return data.frame summarizing the Nested-cross-validation
.CrossValidationLoop <- function(x, y, 
                                 CV.folds, 
                                 y.label, y.covariates = NULL, y.sample = NULL,
                                 method = c("LR", "signature"),
                                 hyperparams = NULL,
                                 hyperparams.best.per.fold = NULL,
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
  
  if (!is.null(DEA.data)){
    if (any(dim(DEA.data) != dim(x))){
      warning("DEA.data not the same format as x. Lets try to harmonize the data")
      tryCatch(expr = {
        DEA.data <- DEA.data[rownames(x), colnames(x)]
      }, error = function(cond){
        stop("DEA.data should have the same columns and rows as x")
      })
    }
  }
  
  # Get the names of the outer folds: 
  names.outer <- names(CV.folds)
  
  # Set the hyper-parameter list
  # This allows to call this function with hyperparams.best.per.fold or with hyperparams
  hyperparams.list <- list()
  if (!is.null(hyperparams.best.per.fold)){
    for (fold.out.name in names.outer){
      list.tmp <- list()
      for (hp in colnames(hyperparams.best.per.fold)){
        list.tmp[[hp]] <- c(hyperparams.best.per.fold[fold.out.name, hp])
        if (is.na(list.tmp[[hp]])){
          if (method == "LR"){
            list.tmp[[hp]] <- DEFAULT.LR.HYPERPARAMS[[hp]]
          } else if (method == "signature"){
            list.tmp[[hp]] <- DEFAULT.SIGNATURE.HYPERPARAMS[[hp]]
          }
        }
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
      message(paste0("fold.out.name: ", fold.out.name))
      
      # Get the indexes of outer fold
      fold.out <- CV.folds[[fold.out.name]]
      
      # Get the outer training and testing data
      x.train.outer <- x[fold.out$train,, drop=F]
      y.train.outer <- y[fold.out$train,, drop=F]
      x.test.outer  <- x[fold.out$test,, drop=F] 
      y.test.outer  <- y[fold.out$test,, drop=F]
      
      # Get the inner and outer DEA data
      DEA.data.train.outer <- NULL
      if (!is.null(DEA.data)){
        DEA.data.train.outer <- DEA.data[fold.out$train,, drop=F]
      }
      
      hyperparams_ <- hyperparams.list[[fold.out.name]]
      
      if (!is.null(folds.save.folder)){
        folds.save.folder.tmp <- paste0(folds.save.folder, "CV_", fold.out.name, "/")
        if (!dir.exists(folds.save.folder.tmp)){
          system(command = paste("mkdir", "-p", folds.save.folder.tmp))
        }
      } else {
        folds.save.folder.tmp <- NULL
      }
      
      res.model <- TrainModel(
        x.train = x.train.outer,
        y.train = y.train.outer,
        x.test = x.test.outer,
        y.test = y.test.outer,
        y.label = y.label,
        y.covariates = y.covariates,
        y.sample = y.sample,
        method = method,
        hyperparams = hyperparams_,
        data.trans.method = data.trans.method,
        pca.explained.var.threshold = pca.explained.var.threshold,
        rm.corr.features = rm.corr.features,
        DA.method = DA.method,
        DA.p.adjust.method = DA.p.adjust.method,
        DA.p.val.threshold = DA.p.val.threshold,
        DA.event.per.variable = DA.event.per.variable,
        LR.prob.threhsold = LR.prob.threhsold,
        signature.col.aggregate = signature.col.aggregate,
        DEA.method = DEA.method,
        DEA.data = DEA.data.train.outer,
        sample.weights = sample.weights,
        folds.save.folder = folds.save.folder.tmp, 
        save.fold.model = save.fold.model, 
        save.fold.pred = save.fold.pred, 
        save.fold.score = save.fold.score)
      
      
      message(paste0("Finished fold.out.name: ", fold.out.name))
      
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
#' @param CV.res data.frame; Nested-Cross-Validation results as returned by .CrossValidationLoop() or .NestedCrossValidationLoop()
#' @param method characters; The method of the origin model. possible values are  
#' "LR", "signature"
#' @param outer.col character; The name of the outer-fold column in outer.col
#' 
#' @return data.frame; Data frame summarizing the best hyper-parameters for each outerfold
#' 
#' @export
GetBestHyperparams <- function(CV.res, 
                               method = c("LR", "signature"),
                               CV.type = c("CV", "NCV"),
                               outer.col = "outer", 
                               selection.metric = EVALUATION.METRICS){
  
  method <- match.arg(method)
  CV.type <- match.arg(CV.type)
  selection.metric <- match.arg(selection.metric)
  
  
  if (method == "LR"){
    hyperparms.col = names(DEFAULT.LR.HYPERPARAMS)
  } else if (method == "signature"){
    hyperparms.col = names(DEFAULT.SIGNATURE.HYPERPARAMS)
  }
  
  if (!(paste0("test.", selection.metric) %in% colnames(CV.res))){
    stop("GetBestOuterHyperparams: The metric '", selection.metric, "' not in CV.res")
  }
  
  CV.res$accuracy_ <- CV.res[[paste0("test.", selection.metric)]]
  
  if (CV.type == "CV"){
    group.columns <- c(hyperparms.col)
  } else if (CV.type == "NCV"){
    group.columns <- c(outer.col, hyperparms.col)
  }
  
  CV.res.avg <- CV.res %>% 
    filter(!is.na(accuracy_)) %>%
    group_by(across(all_of(group.columns))) %>% 
    summarise(acc_avg = mean(accuracy_),
              acc_max = max(accuracy_),
              acc_min = min(accuracy_),
              .groups = "keep")
  CV.res.avg <- data.frame(CV.res.avg)
  
  if (CV.type == "CV"){
    hyperparams.best <- CV.res.avg[CV.res.avg$acc_avg == max(CV.res.avg$acc_avg),, drop = F]
    if (nrow(hyperparams.best) > 1){
      warning(": hyperparams.best has multiple values possible. Continue with the first one \n", hyperparams.best)
      hyperparams.best <- hyperparams.best[1, , drop=F]
    }
  } else if (CV.type == "NCV"){
    hyperparams.best <- CV.res.avg %>% 
      group_by(across(all_of(outer.col))) %>% 
      filter(acc_avg == max(acc_avg))
    hyperparams.best <- data.frame(hyperparams.best)
    hyperparams.best <- hyperparams.best[!duplicated(hyperparams.best$outer), ]
    rownames(hyperparams.best) <- hyperparams.best$outer
    hyperparams.best$outer <- NULL
    
    possible.outer <- unique(CV.res$outer)
    hyperparams.best <- hyperparams.best[possible.outer, ]
  }

  return(hyperparams.best)
}


NestedCrossValidation <- function(x, y, 
                                  CV.k.in = 5, CV.k.out = 5, CV.y.leave.1.out.column = NULL,
                                  CV.replicates = 1,
                                  CV.sampling.method = c("under", "upper", "none"),
                                  y.label, y.covariates = NULL, y.sample = NULL,
                                  method = c("LR", "signature"),
                                  metric.to.optimize = EVALUATION.METRICS,
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
                                  file.NCV.inner = NULL, 
                                  file.NCV.outer = NULL,
                                  seed = 1234){
  
  CV.sampling.method <- match.arg(CV.sampling.method)
  method <- match.arg(method)
  metric.to.optimize <- match.arg(metric.to.optimize)
  data.trans.method <- match.arg(data.trans.method)
  DA.method <- match.arg(DA.method)
  DA.p.adjust.method <- match.arg(DA.p.adjust.method)
  DEA.method <- match.arg(DEA.method)
  
  cat("Nested-Cross-Validation - Create folds: \n")
  
  # Define the NCV CV.y.leave.1.out parameter if it exists:
  if (!is.null(CV.y.leave.1.out.column)){
    CV.y.leave.1.out <- y[, CV.y.leave.1.out.column]
  } else {
    CV.y.leave.1.out <- NULL
  }
  # Get folds:
  CV_folds <- CreateNestedCVFolds(
    y = y[, y.label],
    y.leave.1.out = CV.y.leave.1.out, 
    k.out = CV.k.out, 
    k.in = CV.k.in,
    replicates = CV.replicates,
    sampling.method = CV.sampling.method, 
    seed = seed)
  
  cat("Nested-Cross-Validation - Inner Loop: \n")
  
  if (!is.null(file.NCV.inner)){
    cat(paste0("Reading Inner NCV results form ", file.NCV.inner, "\n"))
    NCV.inner.res <- readRDS(file = file.NCV.inner)
  } else {
    # Perform the NSC
    NCV.inner.res <- .NestedCrossValidationLoop(
      x = x, 
      y = y,
      y.label = y.label,
      y.sample = y.sample,
      y.covariates = y.covariates,
      folds.outer = CV_folds$folds.outer,
      folds.outer.inner = CV_folds$folds.outer.inner,
      method = method,
      hyperparams = hyperparams,
      data.trans.method = data.trans.method,
      pca.explained.var.threshold = pca.explained.var.threshold,
      rm.corr.features = rm.corr.features,
      DA.method = DA.method, 
      DA.p.adjust.method = DA.p.adjust.method,
      DA.p.val.threshold = DA.p.val.threshold, 
      DA.event.per.variable = DA.event.per.variable,
      LR.prob.threhsold = LR.prob.threhsold,
      signature.col.aggregate = signature.col.aggregate,
      DEA.method = DEA.method,
      DEA.data = DEA.data,
      sample.weights = sample.weights,
      folds.save.folder = folds.save.folder, 
      save.fold.model = F, 
      save.fold.pred = !is.null(folds.save.folder),
      save.fold.score = F
    )
  }
  
  
  cat("nested-Cross-Validation - Outer Loop: \n")
  
  cat("Get best hyperparameters: \n")
  
  # Get the best hyperparameters: 
  hyperparams.best <- GetBestHyperparams(
    CV.res = NCV.inner.res,
    method = method,
    CV.type = "NCV",
    selection.metric = metric.to.optimize)
  
  cat("Run outer loops (CV): \n")
  
  if (!is.null(file.NCV.outer)){
    # Read the NCV results
    cat(paste0("Reading Outer NCV results form ", file.NCV.outer, "\n"))
    NCV.outer.res <- readRDS(file = file.NCV.outer)
  } else {
    NCV.outer.res <- .CrossValidationLoop(
      x = x,
      y = y,
      y.label = y.label,
      y.sample = y.sample,
      y.covariates = y.covariates,
      CV.folds = CV_folds$folds.outer,
      method = method,
      hyperparams = NULL,
      hyperparams.best.per.fold = hyperparams.best,
      data.trans.method = data.trans.method,
      pca.explained.var.threshold = pca.explained.var.threshold,
      rm.corr.features = rm.corr.features,
      DA.method = DA.method, 
      DA.p.adjust.method = DA.p.adjust.method,
      DA.p.val.threshold = DA.p.val.threshold, 
      DA.event.per.variable = DA.event.per.variable,
      LR.prob.threhsold = LR.prob.threhsold,
      signature.col.aggregate = signature.col.aggregate,
      DEA.method = DEA.method,
      DEA.data = DEA.data,
      sample.weights = sample.weights,
      folds.save.folder = folds.save.folder, 
      save.fold.model = !is.null(folds.save.folder), 
      save.fold.pred = !is.null(folds.save.folder),
      save.fold.score = !is.null(folds.save.folder)
    )
  }
  
  return(list("NCV.outer.res" = NCV.outer.res, "NCV.inner.res" = NCV.inner.res))
}


CrossValidation <- function(x, y, 
                            CV.k = 5, 
                            CV.y.leave.1.out.column = NULL,
                            CV.replicates = 1,
                            CV.sampling.method = c("under", "upper", "none"),
                            y.label, y.covariates = NULL, y.sample = NULL,
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
                            file.CV = NULL,
                            seed = 1234){
  
  
  CV.sampling.method <- match.arg(CV.sampling.method)
  method <- match.arg(method)
  data.trans.method <- match.arg(data.trans.method)
  DA.method <- match.arg(DA.method)
  DA.p.adjust.method <- match.arg(DA.p.adjust.method)
  DEA.method <- match.arg(DEA.method)
  
  cat("Cross-Validation - Create folds: \n")
  
  # Define the NCV CV.y.leave.1.out parameter if it exists:
  if (!is.null(CV.y.leave.1.out.column)){
    CV.y.leave.1.out <- y[, CV.y.leave.1.out.column]
  } else {
    CV.y.leave.1.out <- NULL
  }
  # Get folds:
  CV.folds <- CreateCVFolds(
    y = y[, y.label],
    y.leave.1.out = CV.y.leave.1.out, 
    k = CV.k, 
    replicates = CV.replicates,
    sampling.method = CV.sampling.method, 
    seed = seed)
  
  cat("Cross-Validation - Loop: \n")
  
  if (!is.null(file.CV)){
    cat(paste0("Reading CV results form ", file.CV, "\n"))
    CV.res <- readRDS(file = file.CV)
  } else {
    
    CV.res <- .CrossValidationLoop(
      x = x, 
      y = y,
      y.label = y.label,
      y.sample = y.sample,
      y.covariates = y.covariates,
      CV.folds = CV.folds,
      method = method,
      hyperparams = hyperparams,
      hyperparams.best.per.fold = NULL,
      data.trans.method = data.trans.method,
      pca.explained.var.threshold = pca.explained.var.threshold,
      rm.corr.features = rm.corr.features,
      DA.method = DA.method, 
      DA.p.adjust.method = DA.p.adjust.method,
      DA.p.val.threshold = DA.p.val.threshold, 
      DA.event.per.variable = DA.event.per.variable,
      LR.prob.threhsold = LR.prob.threhsold,
      signature.col.aggregate = signature.col.aggregate,
      DEA.method = DEA.method, 
      DEA.data = DEA.data,
      sample.weights = sample.weights,
      folds.save.folder = folds.save.folder, 
      save.fold.model = F, # It saves as many model as combination of hyperparmeters 
      save.fold.pred = !is.null(folds.save.folder), # Only saves one data.frame
      save.fold.score = !is.null(folds.save.folder) # Only saves one data.frame
    )
  }
  
  return(CV.res)
}
