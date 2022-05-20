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



#' Train Model function
#' 
#' Home function to train and test a model (LR or signature )
#' 
#' @param x.train data.frame or matrix; The input matrix corresponding to the feature 
#' space only. Additional covariates are found in `y`. Rows = samples, Columns = features
#' @param y.train data.frame or matrix; The output matrix. 
#' Columns = c(y.label, y.covariates, y.sample)
#' @param  x.test data.frame or matrix; The input data like x.train but for testing
#' @param  y.test data.frame or matrix; The output matrix like y.train but for testing
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
#' @param save.DEA.folder character; Folder to save the DEA information in
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
                       DA.method = c("none", DA.METHODS),
                       DA.p.adjust.method = P.ADJUST.METHODS,
                       DA.p.val.threshold = 0.05,
                       LR.prob.threhsold = 0.5,
                       DEA.method = DEA.METHODS, 
                       save.DEA.folder){
  
  # Setting parameters: 
  method <- match.arg(method) # default = "LR"
  data.trans.method <- match.arg(data.trans.method) # default = "none"
  DA.method <- match.arg(DA.method) # default = "none"
  DA.p.adjust.method <- match.arg(DA.p.adjust.method) # default = "fdr"
  DEA.method <- match.arg(DEA.method) # default = "wilcox"

  if (nrow(x.test) == 0){
    x.test <- NULL
    y.test <- NULL
  }
  
  # Feature transformation
  if (data.trans.method != "none"){
    data.transformed <- FeatureTrans(
      x.train = x.train, 
      x.test = x.test, 
      method = data.trans.method,
      explained.var.threshold = pca.explained.var.threshold)
    
    x.train <- data.transformed$x.train
    x.test <- data.transformed$x.test
    
    rm(data.transformed)
  }
  
  # Discriminant analysis
  if (DA.method != "none"){
    discr.res <- DiscrAnalysisBinary(
      x = x.train, 
      y = y.train[, y.label], 
      method = DA.method, 
      p.adjust.method = DA.p.adjust.method, 
      p.val.threshold = DA.p.val.threshold
    )
    sig.features <- rownames(subset(discr.res, Significant))
    
    if (length(sig.features) > 0){
      x.train <- x.train[, sig.features, drop = F]
      if (!is.null(x.test)){
        x.test <- x.test[, sig.features, drop = F]
      }
    } else {
      warning("Discriminant Analysis: No significant features were identified. We continue without Discriminant Analysis")
    }
  }
  
  if (method == "LR"){
    
    # Prepare data
    # For LR, the input (x) and output (y) are in the same data.frame
    design.str <- paste(y.label, "~", paste(c(colnames(x.train), y.covariates), collapse = " + "), sep = " ")
    
    data.train <- cbind(x.train, y.train)
    data.test <- cbind(x.test, y.test)
    
    # Run the training
    res.model <- LRCrossValidation(
      data.train = data.train, 
      data.test = data.test,
      design.str = design.str, 
      hyperparams = hyperparams, 
      model.prob.threshold = LR.prob.threhsold)
    
  } else if (method == "signature"){
    # Create Seurat object
    data.train <- Seurat::CreateSeuratObject(counts = t(x.train), assay = "RNA", project = "NCV", meta.data = y.train)
    if(!is.null(x.test)){
      data.test  <- Seurat::CreateSeuratObject(counts = t(x.test), assay = "RNA", project = "NCV", meta.data = y.test)
    } else {
      data.test <- NULL
    }
    
    # Get the results of the model for all hyperparmeters
    res.model <- SignatureCrossValidation(
      data.train = data.train, 
      data.test = data.test, 
      y.label = y.label, 
      y.sample = y.sample, 
      y.covariates = y.covariates, 
      DEA.method = DEA.method,
      signature.lengths = hyperparams$signature.lengths, 
      signature.sides = hyperparams$signature.side, 
      signature.rm.regex = hyperparams$signature.rm.regex, 
      signature.methods = hyperparams$signature.methods, 
      assay = "RNA", 
      slot = "counts",
      save.DEA.folder = save.DEA.folder)
  }
  
  return(res.model)
}


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
        DA.method = DA.method,
        DA.p.adjust.method = DA.p.adjust.method,
        DA.p.val.threshold = DA.p.val.threshold,
        LR.prob.threhsold = LR.prob.threhsold,
        DEA.method = DEA.method,
        save.DEA.folder = NULL)
      
      # Add to the results additional feature dimention and properties
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
    foreach::foreach(fold.out.name = names.outer, .combine='rbind') %dopar% { # fold.out.name <- names.outer[11]
      message(paste0("fold.out.name: ", fold.out.name))
      
      # Get the indexes of outer fold
      fold.out <- folds.outer[[fold.out.name]]
      
      # Get the outer training and testing data
      x.train.outer <- x[fold.out$train,, drop=F]
      y.train.outer <- y[fold.out$train,, drop=F]
      x.test.outer  <- x[fold.out$test,, drop=F] 
      y.test.outer  <- y[fold.out$test,, drop=F]
      
      hyperparams_ <- hyperparams.list[[fold.out.name]]
      
      res.model <- TrainModel(
        x.train = x.train.outer, # [, 1:1000],
        y.train = y.train.outer,
        x.test = x.test.outer, # [, 1:1000],
        y.test = y.test.outer,
        y.label = y.label,
        y.covariates = y.covariates,
        y.sample = y.sample,
        method = method,
        hyperparams = hyperparams_,
        data.trans.method = data.trans.method,
        pca.explained.var.threshold = pca.explained.var.threshold,
        DA.method = DA.method,
        DA.p.adjust.method = DA.p.adjust.method,
        DA.p.val.threshold = DA.p.val.threshold,
        LR.prob.threhsold = LR.prob.threhsold,
        DEA.method = DEA.method,
        save.DEA.folder = NULL)
      
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
