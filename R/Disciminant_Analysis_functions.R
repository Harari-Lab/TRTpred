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
#' 
#' @return data.frame. The results form the discriminant analysis. 
#'  - row.names = features
#'  - column names = "p_val" (numerical), "p_val_adj" (numerical) & "Significant" (logical)
#' 
#' @export
DiscrAnalysisBinary <- function(x, y, method = DA.METHODS,
                                p.adjust.method = P.ADJUST.METHODS,
                                p.val.threshold = 0.05){
  
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
      WRS.res <- stats::wilcox.test(x = x.T[, col], y = x.F[, col], correct = F)
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
    
    # get which are significant
    test.res$Significant <- test.res$p_val_adj <= p.val.threshold
  }
  
  return(test.res)
}
