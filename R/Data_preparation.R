# Function to prepare data for machine learning models

suppressMessages(require(stats))
suppressMessages(require(ropls))

FEATURE.TRANS.METHODS <- c("pca", "opls")
DA.METHODS <- c("wilcox")
P.ADJUST.METHODS <- c("fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "none")

#' Prepare data for machine learning
#' 
#' Function which compiles the different steps to prepare the input data for some model. 
#' The steps are:
#' 1. Data Transformation (data.trans.method)
#' 2. Remove Correlated Features (rm.corr.features)
#' 3. Discriminant analysis (DA.method)
#' 
#' Note that every step can be skipped by setting the corresponding parameters 
#' to "none", False, "none" respectively. 
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
#' @param DEA.data data.frame or matrix; The input matrix for the DEA. 
#' Rows = samples, Columns = features. If NULL, the DEA.data is x
#' Default = NULL
#' @param data.trans.method character: The data transforamtion method. 
#' Possible values are "pca" 
#' Default = "none" i.e. no transformation. 
#' @param pca.explained.var.threshold numeric: If not NULL, apply this threshold 
#' to select PCs explaining more variance than the threhsold
#' @param rm.corr.features logical: Do we remove the mutually correlated feature 
#' while keeping the ones that best correlate with the outcome? 
#' default = F
#' @param rm.corr.features.threhsold numerical: The correlation threshold of rm.corr.features
#' Default = 0.8
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
#' 
#' @return data.frame summarizing the Nested-cross-validation
#' 
#' @export
PrepareData <- function(
    x.train, y.train, y.label, 
    x.test = NULL, y.test = NULL,
    DEA.data = NULL,
    method = c("LR", "signature"),
    data.trans.method = c("none", FEATURE.TRANS.METHODS), 
    pca.explained.var.threshold = NULL,
    rm.corr.features = F, 
    rm.corr.features.threhsold = 0.8,
    DA.method = c("none", DA.METHODS), 
    DA.p.adjust.method = P.ADJUST.METHODS, 
    DA.p.val.threshold = 0.05, 
    DA.event.per.variable = NULL
    ){
  
  # TODO DEA.data implementation?  
  
  method <- match.arg(method)
  data.trans.method <- match.arg(data.trans.method)
  DA.method <- match.arg(DA.method)
  DA.p.adjust.method <- match.arg(DA.p.adjust.method)
  
  # Feature transformation
  if (data.trans.method != "none"){
    data.transformed <- FeatureTrans(
      x.train = x.train, 
      x.test = x.test,
      y.train.lab = y.train[, y.label],
      method = data.trans.method,
      explained.var.threshold = pca.explained.var.threshold)
    
    x.train <- data.transformed$x.train
    x.test <- data.transformed$x.test
    
    rm(data.transformed)
  }
  
  # Remove Mutually correlated features
  if (rm.corr.features){
    corr.res <- RemoveMutuallyCorrFeatures(
      x = x.train,  
      y = y.train[, y.label], 
      threshold = rm.corr.features.threhsold)
    
    keep.features <- rownames(subset(corr.res, keep))
    if (length(keep.features) > 0){
      x.train <- x.train[, keep.features, drop = F]
      if (!is.null(x.test)){
        x.test <- x.test[, keep.features, drop = F]
      }
    } else {
      warning("Remove Mutually Correlated features: No mutuallly corr features were identified.")
    }
  }
  
  # Discriminant analysis
  if (DA.method != "none"){
    discr.res <- DiscrAnalysisBinary(
      x = x.train, 
      y = y.train[, y.label], 
      method = DA.method, 
      p.adjust.method = DA.p.adjust.method, 
      p.val.threshold = DA.p.val.threshold,
      event.per.variable = DA.event.per.variable
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
  
  if (method == "signature"){
    # Create Seurat object
    # data.train:
    data.train <- Seurat::CreateSeuratObject(counts = t(x.train), assay = "RNA", project = "NCV", meta.data = y.train)
    # data.test:
    if(!is.null(x.test)){
      data.test <- Seurat::CreateSeuratObject(counts = t(x.test), assay = "RNA", project = "NCV", meta.data = y.test)
    } else {
      data.test <- NULL
    }
    # DEA.data:
    if (!is.null(DEA.data)){
      DEA.data <- Seurat::CreateSeuratObject(counts = t(DEA.data), assay = "RNA", project = "NCV", meta.data = y.train)
    } 
  } else if (method == "LR"){
    # Prepare data
    # For LR, the input (x) and output (y) are in the same data.frame
    data.train <- cbind(x.train, y.train)
    if (!is.null(x.test)){
      data.test <- cbind(x.test, y.test)
    }
    if (!is.null(DEA.data)){
      DEA.data <- cbind(DEA.data, y.train)
    } else {
      data.test <- NULL
    }
  } else {
    data.train <- cbind(x.train, y.train)
    data.train <- cbind(x.train, y.train)
    if (!is.null(x.test)){
      data.test <- cbind(x.test, y.test)
    } else {
      data.test <- NULL
    }
    if (!is.null(DEA.data)){
      DEA.data <- cbind(DEA.data, y.train)
    }
  }
  
  
  # Get return list
  return.list <- list()
  return.list[["data.train"]] <- data.train
  return.list[["data.test"]] <- data.test
  return.list[["DEA.data"]] <- DEA.data
  
  return(return.list)
}


#' Feature space transformation function
#' 
#' Function to transforme the feature space into another feature space. Example 
#' with PCA
#' 
#' @param x.train data.frame or matrix; The training input data. Train the PCA on it
#' @param x.train data.frame or matrix; The testing input data. Apply (predict) the trained PCA on it
#' 
#' @param method character; The data transformation method.  
#' possible values are : pca
#' Default = pca
#' @param explained.var.threshold numerical; If method = "pca". Threshold to 
#' select the PCs with explained variable above the threshold
#' Default = NULL (no threshold)
#' 
#' @return list of traingin and testing data.frames of transformed data
#' 
#' @export
FeatureTrans <- function(x.train, x.test = NULL, 
                         y.train.lab = NULL,
                         method = FEATURE.TRANS.METHODS, 
                         explained.var.threshold = NULL){
  
  method <- match.arg(method)
  
  if (method == "pca"){
    pca.res <- stats::prcomp(x.train, center = F, scale. = F)
    x.train.res <- data.frame(pca.res$x, check.rows = F, check.names = F)
    if (!is.null(x.test)){
      x.test.res <- stats::predict(object = pca.res, newdata = x.test)
    } else {
      x.test.res <- NULL
    }
    
    if (!is.null(explained.var.threshold)){
      pca.res <- summary(pca.res)
      pca.pcs.sig <- names(pca.res$importance["Proportion of Variance", pca.res$importance["Proportion of Variance", ] >= explained.var.threshold])
      x.train.res <- x.train.res[, pca.pcs.sig]
      if (!is.null(x.test)){
        x.test.res <- x.test.res[, pca.pcs.sig]
      }
    }
  } else if (method == "opls"){
    y.train.lab <- as.factor(y.train.lab)
    y.train.lab <- factor(y.train.lab, levels = c(levels(y.train.lab)[2], levels(y.train.lab)[1]))
    
    opls.res <- ropls::opls(
      x = x.train, 
      y = y.train.lab, 
      scaleC = 'standard', # c("none", "center", "pareto", "standard")
      log10L = FALSE,
      predI = 1,
      orthoI = 1,
      fig.pdfC = "none", # no figure is generated
      info.txtC = "none") # no save file name is provided)
    
    opls.res.df <- cbind(opls.res@weightMN, opls.res@weightStarMN, opls.res@orthoWeightMN,
                         opls.res@vipVn, opls.res@orthoVipVn, 
                         opls.res@loadingMN, opls.res@orthoLoadingMN, 
                         opls.res@coefficientMN)
    colnames(opls.res.df) <- c("weightMN", "weightStarMN", "orthoWeightMN",
                               "vipVn", "orthoVipVn", 
                               "loadingMN", "orthoLoadingMN", 
                               "coefficientMN") 
    # Add to this data.frame the variables that were removed due to having a 
    # variance < 2.2e-16 in the full or partial (cross-validation) dataset
    ZeroVarVi <-  data.frame(row.names = names(opls.res@xZeroVarVi))
    if (nrow(ZeroVarVi) > 0){
      for (col_ in colnames(opls.res.df)){ZeroVarVi[[col_]] <- NA}
      opls.res.df <- rbind(opls.res.df, ZeroVarVi[, colnames(opls.res.df)])
    }
    
  }
  
  return(list("x.train" = x.train.res, "x.test" = x.test.res))
}


#' Discriminant analysis function
#' 
#' @description 
#' This function compute the discriminant analysis of some feature space when 
#' comparing two classes in some categorical variable.
#' The discriminant analysis methods implemented are "wilcoxon rank sum test" 
#' 
#' @param x data.frame; The input feature space data
#' @param y categorical vector; The binary output
#' @param method character; The discriminant analysis method. Possible values are 
#' wilcox. 
#' @param p.adjust.method, character; The name of the method to adjust the p-value.
#' possible values are methods in `stats::p.adjust()` i.e. "bonferroni", "fdr", 
#' "holm", "hochberg", "hommel", "BH", "BY", "none".
#' Default = "fdr"
#' @param p.val.threshold, numerical; The p-val significant threshold 
#' Default = 0.05
#' @param event.per.variable, numerical; The number of event per variable threshold
#' A common EPV is 10 meaning that we need 10 event per variable. Here we use it as 
#' a way to select a fewer amount of variables for a fixed amount of samples. 
#' Default = NULL (no selection on EPV)
#' 
#' @return data.frame. The results form the discriminant analysis. 
#'  - row.names = features
#'  - column names = "p_val" (numerical), "p_val_adj" (numerical) & "Significant" (logical)
#' 
#' @export
DiscrAnalysisBinary <- function(x, y, method = DA.METHODS,
                                p.adjust.method = P.ADJUST.METHODS,
                                p.val.threshold = 0.05, 
                                event.per.variable = NULL){
  
  method <- match.arg(method) # default = wilcox
  p.adjust.method <- match.arg(p.adjust.method) #default = "fdr"
  
  y <- as.factor(y)
  classes <- levels(y)
  
  x <- x[, colSums(x) > 0]
  
  if (method == "wilcox"){
    # separate the x between the two classes
    x.T <- x[y == classes[1], ]
    x.F <- x[y == classes[2], ]
    
    # Get the raw p-values
    test.res <- sapply(X = colnames(x), function(col){
      WRS.res <- stats::wilcox.test(x = x.T[, col], y = x.F[, col], 
                                    alternative = "two.sided", 
                                    correct = F, exact = F, paired = F)
      # list("PC" = x, "statistic" = WRS.res$statistic, "p.value" = WRS.res$p.value)
      WRS.res$p.value
    })
    
    test.res <- data.frame(test.res)
    colnames(test.res) <- c("p_val")
    
    # correct the p-value  
    if (p.adjust.method != "none"){
      test.res$p_val_adj <- stats::p.adjust(p = test.res$p_val, method = p.adjust.method)
    } else {
      test.res$p_val_adj <- test.res$p_val
    }
    
    if (!is.null(event.per.variable)){
      # Get the number of features
      num.feature <- ceiling(nrow(x)/event.per.variable)
      # Identify features that pass the EPV threhsold
      # here we sort the test.res from most to least significant and take the num.feature first feature
      test.res <- test.res[order(test.res$p_val_adj), ]
      test.res$EPV_selection <- F
      test.res$EPV_selection[1:num.feature] <- T
      # Get which are significant
      test.res$Significant <- (test.res$p_val_adj <= p.val.threshold) & test.res$EPV_selection
    } else {
      # Get which are significant
      test.res$Significant <- test.res$p_val_adj <= p.val.threshold
    }
  }
  
  return(test.res)
}

#' Remove mutually correlated features
#' 
#' @description 
#' This function outputs a data.frame indicating which features are mutually 
#' correlated and which one to remove.
#' 
#' @param x data.frame; The input feature space data
#' @param y categorical vector; The binary output
#' @param threshold numerical; Significant correlaltion threhsolld
#' default = 0.8
#' 
#' @return data.frame. The results form the discriminant analysis. 
#'  - row.names = features
#'  - column names = "mutually_corr" (logical), "keep" (logical) & "abs_y_corr" (numerical)
#' 
#' @export
RemoveMutuallyCorrFeatures <- function(x, y, threshold = 0.8){
  
  y <- as.factor(y)
  
  corr.res <- data.frame(row.names = colnames(x))
  corr.res$mutually_corr <- F
  corr.res$keep <- T
  corr.res$abs_y_corr <- abs(cor(x = x, y = (y == levels(y)[1]), method = "pearson"))
  
  corr.mat <- cor(x, method = "pearson")
  corr.mat.bool <- abs(corr.mat) >= threshold
  corr.features <- rownames(corr.mat.bool)[rowSums(corr.mat.bool) > 1]
  
  if (length(corr.features) > 0){
    corr.res[corr.features, ]$mutually_corr <- T
    
    # feat.keep <- c()
    feat.rm <- c()
    for (corr.feat in corr.features){
      mutually.corr.feat <- rownames(corr.mat.bool)[corr.mat.bool[corr.feat, ]]
      best.corr.feat <- names(which(corr.res$abs_y_corr[mutually.corr.feat, 1] == 
                                      max(corr.res$abs_y_corr[mutually.corr.feat, 1])))[1]
      # feat.keep <- c(feat.keep, best.corr.feat)
      feat.rm <- c(feat.rm, mutually.corr.feat[mutually.corr.feat != best.corr.feat])
    }
    corr.res[feat.rm, ]$keep <- F
  }
  
  return(corr.res)
}



