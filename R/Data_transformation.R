# Function to transforme datta

suppressMessages(require(stats))
suppressMessages(require(ropls))

FEATURE.TRANS.METHODS <- c("pca", "opls")

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
