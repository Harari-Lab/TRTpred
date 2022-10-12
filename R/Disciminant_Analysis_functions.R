# Discriminant analysis functions

DA.METHODS <- c("wilcox")
P.ADJUST.METHODS <- c("fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "none")

suppressMessages(require(stats))

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
