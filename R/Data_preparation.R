# Data Prepation functions
# Author: Rémy Pétremand
# Date: 07.05.2024
# Description: Function to prepare data for machine learning models
# Reference: https://doi.org/10.1038/s41587-024-02232-0

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

suppressMessages(require(stats))
suppressMessages(require(ropls))

# ------------------------------------------------------------------------------
# Global Parameters
# ------------------------------------------------------------------------------

FEATURE.TRANS.METHODS <- c("pca", "opls")
# DA.METHODS <- c("wilcox")
P.ADJUST.METHODS <- c("none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY")
LIMMA.METHODS <- c("limma_voom", "limma_trend")
EDGER.METHODS <- c("edgeR_LRT", "edgeR_QFL")
DESEQ.METHODS <- c("DESeq2_Wald", "DESeq2_LRT")
SEURAT.METHODS <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR")
DEA.METHODS <- c(SEURAT.METHODS, DESEQ.METHODS, EDGER.METHODS, LIMMA.METHODS)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

#' Prepare Data for Machine Learning
#' 
#' This function compiles the various steps to prepare input data for a model. 
#' The steps are:
#' 1. Data Transformation (data.trans.method)
#' 2. Remove Correlated Features (rm.corr.features)
#' 3. Discriminant Analysis (DA.method)
#' 
#' Note that each step can be skipped by setting the corresponding parameters 
#' to "none", FALSE, and "none" respectively. 
#' 
#' @param x.train data.frame or matrix; The input matrix corresponding to the feature 
#' space only. Additional covariates are found in `y`. Rows represent samples, Columns represent features.
#' @param y.train data.frame or matrix; The output matrix. 
#' Columns = c(y.label, y.covariates, y.sample)
#' @param x.test data.frame or matrix; The input data similar to x.train but for testing.
#' @param y.test data.frame or matrix; The output matrix similar to y.train but for testing.
#' @param y.label character; The column of `y` corresponding to the binary output.
#' @param y.covariates character vector; The columns of `y` corresponding to the covariates.
#' @param y.sample character; The column of `y` corresponding to the samples when method = "signature".
#' @param method character; The prediction method. Possible values are: "LR" (Logistic regression) or "signature" (signature-score approach).
#' @param DEA.data data.frame or matrix; The input matrix for the DEA. 
#' Rows represent samples, Columns represent features. If NULL, DEA.data is set to x.
#' Default = NULL.
#' @param data.trans.method character; The data transformation method. 
#' Possible values are "pca". Default = "none" (i.e., no transformation).
#' @param pca.explained.var.threshold numeric; If not NULL, apply this threshold 
#' to select PCs explaining more variance than the threshold.
#' @param rm.corr.features logical; Should mutually correlated features be removed 
#' while keeping the ones that best correlate with the outcome? 
#' Default = FALSE.
#' @param rm.corr.features.threshold numeric; The correlation threshold for rm.corr.features.
#' Default = 0.8.
#' @param DA.method character; The discriminant analysis (DA) method. 
#' Possible values are "wilcox". Default = "none" (i.e., no discriminant analysis).
#' @param DA.p.adjust.method character; The p-value adjustment method for the 
#' discriminant analysis. Default = "fdr".
#' @param DA.p.val.threshold numeric; The p-value significance threshold.
#' @param DA.event.per.variable numeric; The number of events per variable threshold.
#' A common EPV is 10, meaning that we need 10 events per variable. Here, it is used 
#' to select fewer variables for a fixed number of samples. 
#' Default = NULL (no selection based on EPV).
#' 
#' @return data.frame summarizing the Nested-cross-validation.
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
    DA.method = c("none", DEA.METHODS),
    y.sample = NULL, 
    y.covariates = NULL, 
    y.aggregate = NULL, 
    DA.p.adjust.method = P.ADJUST.METHODS, 
    DA.p.val.threshold = 0.05, 
    DA.event.per.variable = NULL
    ){
  
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
    
    if (!is.null(DEA.data)){
      DEA.data <- NULL
      warning(paste0("main:PrepareData() - data.trans.method = ", data.trans.method, " and DEA is not null. DEA data is set as NULL and will be y.train in the future steps"))
    }
    
    rm(data.transformed)
  }
  
  # Discriminant analysis
  if (DA.method != "none"){
    
    if (is.null(DEA.data)){
      DEA.data <- x.train
    }
    
    discr.res <- DiscrAnalysisBinary(
      x = DEA.data,
      y = y.train,
      y.label = y.label,
      y.sample = y.sample,
      y.covariates = y.covariates,
      y.aggregate = y.aggregate,
      method = DA.method, 
      p.adjust.method = DA.p.adjust.method,
      p.val.threshold = DA.p.val.threshold,
      event.per.variable = DA.event.per.variable
    )
    sig.features <- subset(discr.res, Significant)$genes
    
    if (length(sig.features) > 0){
      x.train <- x.train[, sig.features, drop = F]
      if (!is.null(x.test)){
        x.test <- x.test[, sig.features, drop = F]
      }
      if (!is.null(DEA.data)){
        DEA.data <- DEA.data[, sig.features, drop = F]
      }
    } else {
      warning("Discriminant Analysis: No significant features were identified. We continue without Discriminant Analysis")
    }
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
      if (!is.null(DEA.data)){
        DEA.data <- DEA.data[, keep.features, drop = F]
      }
    } else {
      warning("Remove Mutually Correlated features: No mutuallly corr features were identified.")
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
    } else {
      data.test <- NULL
    } 
    if (!is.null(DEA.data)){
      DEA.data <- cbind(DEA.data, y.train)
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

#' Feature Space Transformation Function
#' 
#' Function to transform the feature space into another feature space, e.g., with PCA.
#' 
#' @param x.train data.frame or matrix; The training input data. Train the PCA on this data.
#' @param x.test data.frame or matrix; The testing input data. Apply (predict) the trained PCA on this data.
#' @param method character; The data transformation method. Possible values are: "pca".
#' Default = "pca".
#' @param explained.var.threshold numeric; If method = "pca", this threshold is used to 
#' select the PCs with explained variance above the threshold.
#' Default = NULL (no threshold).
#' 
#' @return list; A list of training and testing data.frames of transformed data.
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

#' Discriminant Analysis Function
#' 
#' Computes discriminant analysis of a feature space for comparing two classes in a categorical variable.
#' The implemented discriminant analysis method is Wilcox rank sum test.
#' 
#' @param x data.frame; The input feature space data.
#' @param y categorical vector; The binary output.
#' @param method character; The discriminant analysis method. Possible value is "wilcox".
#' @param p.adjust.method character; The name of the method to adjust the p-value.
#' Possible values are methods in `stats::p.adjust()`, such as "bonferroni", "fdr" (deafault), 
#' "holm", "hochberg", "hommel", "BH", "BY", "none".
#' @param p.val.threshold numeric; The p-value significance threshold. Default = 0.05.
#' @param event.per.variable numeric; The number of events per variable threshold.
#' A common EPV is 10, meaning that we need 10 events per variable. Here, it is used 
#' to select fewer variables for a fixed number of samples. 
#' Default = NULL (no selection based on EPV).
#' 
#' @return data.frame; The results from the discriminant analysis.
#'   - row names = features
#'   - column names = "p_val" (numeric), "p_val_adj" (numeric), and "Significant" (logical)
#' 
#' @export
DiscrAnalysisBinary <- function(x, y, y.label, 
                                y.sample = NULL, 
                                y.covariates = NULL, 
                                y.aggregate = NULL,
                                method = DEA.METHODS,
                                p.adjust.method = P.ADJUST.METHODS,
                                p.val.threshold = 0.05, 
                                event.per.variable = NULL){
  
  method <- match.arg(method) # default = wilcox
  p.adjust.method <- match.arg(p.adjust.method) #default = "fdr"
  
  # Create Seurat object
  x <- Seurat::CreateSeuratObject(counts = t(x), assay = "RNA", project = "NCV", meta.data = y)
  
  # Run DEA
  dea.res <- RunDEA(
    object = x, 
    col.DE_group = y.label, 
    col.sample = y.sample, 
    col.covariate = y.covariates, 
    assay = "RNA", slot = "data", 
    method = method, 
    col.aggregate = y.aggregate)
  
  if (p.adjust.method != "none"){
    dea.res$p_val_adj <- stats::p.adjust(p = dea.res$pval, method = p.adjust.method)
  } else {
    dea.res$p_val_adj <- dea.res$padj # The default is Bonferroni here
  }
  
  if (!is.null(event.per.variable)){
    # Get the number of features
    num.feature <- ceiling(nrow(x)/event.per.variable)
    # Identify features that pass the EPV threhsold
    # here we sort the test.res from most to least significant and take the num.feature first feature
    dea.res <- dea.res[order(dea.res$p_val_adj), ]
    dea.res$EPV_selection <- F
    dea.res$EPV_selection[1:num.feature] <- T
    # Get which are significant
    dea.res$Significant <- (dea.res$p_val_adj <= p.val.threshold) & dea.res$EPV_selection
  } else {
    # Get which are significant
    dea.res$Significant <- dea.res$p_val_adj <= p.val.threshold
  }
  
  # y <- as.factor(y)
  # classes <- levels(y)
  # 
  # x <- x[, colSums(x) > 0]
  # 
  # if (method == "wilcox"){
  #   # separate the x between the two classes
  #   x.T <- x[y == classes[1], ]
  #   x.F <- x[y == classes[2], ]
  #   
  #   # Get the raw p-values
  #   test.res <- sapply(X = colnames(x), function(col){
  #     WRS.res <- stats::wilcox.test(x = x.T[, col], y = x.F[, col], 
  #                                   alternative = "two.sided", 
  #                                   correct = F, exact = F, paired = F)
  #     # list("PC" = x, "statistic" = WRS.res$statistic, "p.value" = WRS.res$p.value)
  #     WRS.res$p.value
  #   })
  #   
  #   test.res <- data.frame(test.res)
  #   colnames(test.res) <- c("p_val")
  #   
  #   # correct the p-value  
  #   if (p.adjust.method != "none"){
  #     test.res$p_val_adj <- stats::p.adjust(p = test.res$p_val, method = p.adjust.method)
  #   } else {
  #     test.res$p_val_adj <- test.res$p_val
  #   }
  #   
  #   if (!is.null(event.per.variable)){
  #     # Get the number of features
  #     num.feature <- ceiling(nrow(x)/event.per.variable)
  #     # Identify features that pass the EPV threhsold
  #     # here we sort the test.res from most to least significant and take the num.feature first feature
  #     test.res <- test.res[order(test.res$p_val_adj), ]
  #     test.res$EPV_selection <- F
  #     test.res$EPV_selection[1:num.feature] <- T
  #     # Get which are significant
  #     test.res$Significant <- (test.res$p_val_adj <= p.val.threshold) & test.res$EPV_selection
  #   } else {
  #     # Get which are significant
  #     test.res$Significant <- test.res$p_val_adj <= p.val.threshold
  #   }
  # }
  
  return(dea.res)
}

#' Remove mutually correlated features
#' 
#' @description 
#' This function outputs a data.frame indicating which features are mutually 
#' correlated and which one to remove.
#' 
#' @param x data.frame; The input feature space data
#' @param y categorical vector; The binary output
#' @param threshold numerical; Significant correlation threshold
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

#' Prepare training data from seurat object
#' 
#' @description 
#' This function prepare training data from Seurat object
#' 
#' @param SeuratObject Seurat object; The Seurat object (required)
#' @param label character; The y-label (i.e tumor reactivtiy) (required)
#' @param label.order character; the order of the y-label to transform them into factor
#' @param assay character; The Seurat assay (default = "RNA")
#' @param slot character; The assay slot (default = "slot")
#' @param DEA.assay character; The Seurat assay used for the Differential Expression Analysis 
#' Default = NULL = same as assay
#' @param DEA.slot character; The assay slot used for the Differential Expression Analysis 
#' Default = NULL = same as slot
#' @param DEA.pseudobulk.col character; The column names in Seurat used to construct pseudobluk for the pseudobluk DEA
#' Default = NULL
#' @param covariates characters; The covariates list
#' Default = NULL
#' @param sample.weights.col character; Which col should be used to normalize the performance
#' Default = NULL = none!
#' 
#' @return list; 
#' - "data.input" 
#' - "data.output" 
#' - "DEA.data" 
#' - "sample.weights" 
#' 
#' @export
PrepareTrainingDataFromSeurat <- function(SeuratObject, label, label.order = NULL,
                                          assay = "RNA", slot = "data", 
                                          DEA.assay = NULL, DEA.slot = "data", DEA.pseudobulk.col = NULL,
                                          covariates = NULL,
                                          sample.weights.col = NULL){
  # ----------------------------------------------------------------------------
  # sanity checks:
  # ----------------------------------------------------------------------------
  if (!(assay %in% names(SeuratObject@assays))){
    stop("assay not in SeuratObject's assays")
  }
  if (is.null(attr(SeuratObject@assays[[assay]], slot))){
    stop("slot not in SeuratObject's slots in the assay selected")
  }
  if (!(label %in% colnames(SeuratObject@meta.data))){
    stop("label not in SeuratObject's meta.data")
  }
  if (!is.null(covariates)){
    if (!all(covariates %in% colnames(SeuratObject@meta.data))){
      stop("some covariates are not in SeuratObject's meta.data")
    }
  }
  if (!is.null(DEA.assay)){
    if (!(DEA.assay %in% names(SeuratObject@assays))){
      stop("DEA.assay not in SeuratObject's assays")
    }
    if (!is.null(DEA.slot)){
      if (is.null(attr(SeuratObject@assays[[DEA.assay]], DEA.slot))){
        stop("slot not in SeuratObject's slots in the DEA.assay selected")
      }
    }
  }
  if (!is.null(DEA.pseudobulk.col)){
    if (!(DEA.pseudobulk.col %in% colnames(SeuratObject@meta.data))){
      stop("DEA.pseudobulk.col not in SeuratObject's meta.data")
    }
  }
  
  # ----------------------------------------------------------------------------
  # Set Label as factor without NAs
  # ----------------------------------------------------------------------------
  
  SeuratObject@meta.data$select <- !is.na(SeuratObject@meta.data[, label])
  SeuratObject <- subset(SeuratObject, subset = select)
  SeuratObject@meta.data$select <- NULL
  
  if (is.null(label.order)){label.order <- unique(SeuratObject@meta.data[, label])}
  SeuratObject@meta.data[, label] <- factor(SeuratObject@meta.data[, label], levels = label.order)
  
  # ----------------------------------------------------------------------------
  # Get input and output data
  # ----------------------------------------------------------------------------
  
  # Get data.input
  data.input <- Seurat::GetAssayData(object = SeuratObject, assay = assay, slot = slot)
  data.input <- data.frame(t(data.input), check.rows = F, check.names = F)
  
  # Get data.output
  data.output <- SeuratObject@meta.data
  
  # Get DEA data
  if (!is.null(DEA.assay)){
    DEA.data <- Seurat::GetAssayData(object = SeuratObject, assay = DEA.assay, slot = DEA.slot)
    DEA.data <- data.frame(t(DEA.data), check.rows = F, check.names = F)
    DEA.data <- DEA.data[rownames(data.input), colnames(data.input)]
    
    if (any(dim(DEA.data) != dim(data.input))){
      warning("DEA.data not the same format as x. Lets try to harmonize the data")
      tryCatch(expr = {
        DEA.data <- DEA.data[rownames(data.input), colnames(data.input)]
      }, error = function(cond){
        stop("DEA.data should have the same columns and rows as x")
      })
    }
    
  } else {
    DEA.data <- NULL
  }
  
  # Get Sample weight
  if (!is.null(sample.weights.col)){
    if (sample.weights.col %in% colnames(SeuratObject@meta.data)){
      # Get group.weight.df
      group.weight.df <- data.frame(table(SeuratObject@meta.data[, sample.weights.col]))
      rownames(group.weight.df) <- group.weight.df[, 1]
      colnames(group.weight.df) <- c(sample.weights.col, "Count")
      # Get weights
      group.weight.df$weight <- sum(group.weight.df$Count)/(group.weight.df$Count * nrow(group.weight.df))
      # Set weights per observation
      SeuratObject@meta.data$weight <- group.weight.df[SeuratObject@meta.data[, sample.weights.col], ]$weight
      
      # Get final sample.weights (rows = observation, columns = weight)
      sample.weights <- SeuratObject@meta.data[, "weight", drop = F]
      
      if (sum(sample.weights$weight) != ncol(SeuratObject)){
        warning("The sum of weights need to equal to the number of observations")
      }
    } else {
      warning("sample.weights.col: '", sample.weights.col, "' is not present in colnames(SeuratObject@meta.data) ",
              "We continue the script without taking this into account")
      sample.weights <- NULL
    } 
  } else {
    sample.weights <- NULL
  }
  
  # ----------------------------------------------------------------------------
  # Return
  # ----------------------------------------------------------------------------
  
  return(list("data.input" = data.input, "data.output" = data.output, "DEA.data" = DEA.data, "sample.weights" = sample.weights))
}
