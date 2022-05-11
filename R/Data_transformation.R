# Function to transforme datta

suppressMessages(require(stats))

FEATURE.TRANS.METHODS <- c("pca")

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
      pca.pcs.sig <- names(pca.res$importance["Proportion of Variance", 
                                              pca.res$importance["Proportion of Variance", ] >= explained.var.threshold])
      x.train.res <- x.train.res[, pca.pcs.sig]
      if (!is.null(x.test)){
        x.test.res <- x.test.res[, pca.pcs.sig]
      }
    }
  }
  
  return(list("x.train" = x.train.res, "x.test" = x.test.res))
}
