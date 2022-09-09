# Signature method
# Rémy Pétremand

suppressMessages(require(Seurat))
suppressMessages(require(AUCell))
suppressMessages(require(UCell))
suppressMessages(require(singscore))
suppressMessages(require(stats))


LIMMA.METHODS <- c("limma_voom", "limma_trend")
EDGER.METHODS <- c("edgeR_LRT", "edgeR_QFL")
DESEQ.METHODS <- c("DESeq2_Wald", "DESeq2_LRT")
SEURAT.METHODS <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR")
DEA.METHODS <- c(SEURAT.METHODS, DESEQ.METHODS, LIMMA.METHODS, EDGER.METHODS)

SIGNATURE.SCORE.METHODS <- c("average", "UCell", "AUCell", "singscore")

SIGNATURE.SELECTION.METHODS <- c("logFC", "pval")

SIGNATURE.SIDE <- c("both", "up", "down")

#' Get signature from DE results
#' 
#' Function to get signature from DE results
#' 
#' @param df data.frame:  DE results, row.names = gene names
#' @param pval_col character: The name of the p adjusted column 
#' @param log2FC_col character: The name of the Log Fold Change column 
#' @param pval_limit numeric: The upper significance limit of the adjusted p-value
#' @param log2FC_limits vector of numeric: The lower and upper limits of the 
#' logFC column to define the up and down regulated genes.
#' @param gene.selection.method character; signature gene selection method.
#' possible values are "logFC" (select from logFC) or "pval" (select from pval)
#' Default = "logFC"
#' @param remove_regex character: Regular expression selecting genes to remove 
#' from the up and down gene-sets
#' 
#' @return list; list of up and down regulated gene. names of the list are 
#'   - up = n up-regulated gene
#'   - down = n down regulated genes
#' 
#' @export
GetSignature <- function(df, n_genes = 20, 
                         pval_col = "padj", log2FC_col = "logFC", 
                         pval_limit = 0.05, log2FC_limits = c(0, 0),
                         gene.selection.method = "logFC",
                         remove_regex = NULL){
  
  # Santity check
  if (!is.null(remove_regex)){
    if (remove_regex == "none"){
      remove_regex <- NULL
    }
  }
  
  # filter df
  df <- df[!is.na(df[[pval_col]]), ]
  df <- df[!is.na(df[[log2FC_col]]), ]
  
  # Keep sig genes
  df <- df[df[[pval_col]] <= pval_limit, ]
  
  # Remove genes recognized by remove_regex
  if (!is.null(remove_regex)){
    mask.2.rm <- grepl(pattern = remove_regex, x = rownames(df))
    df <- df[!mask.2.rm,]
  }
  
  # order the genes
  if (gene.selection.method == "logFC"){
    df <- df[order(df[[log2FC_col]], decreasing = T), ]
  } else if (gene.selection.method == "pval"){
    df <- df[order(df[[pval_col]], decreasing = F), ]
  } else{
    warning(paste0("GetSignature: gene.selection.method = ", gene.selection.method, " not recognized"))
  }
  
  up_genes <- rownames(df[df[[log2FC_col]] >= log2FC_limits[2], ])
  dn_genes <- rownames(df[df[[log2FC_col]] <= log2FC_limits[1], ])
  
  if (gene.selection.method == "logFC"){
    if (length(dn_genes) > 0){
      dn_genes <- dn_genes[length(dn_genes):1]
    }
  }
  
  
  if (n_genes > length(up_genes)){
    warning("n_genes is larger than the number of up regulated genes")
  }
  if (n_genes > length(dn_genes)){
    warning("n_genes is larger than the number of down regulated genes")
  }
  
  return(list(up = up_genes[1:n_genes], down = dn_genes[1:n_genes]))
}

#' Get signature lists from DE results and hyper-parameters
#' 
#' Function to get signature.list from DE results. This functions returns a list 
#' of signatures got from GetSignature() with varying lengths or side or even 
#' with and without certain genes. The name of the resulting list corresponds to 
#' the hyperparameters mentionned above. 
#' 
#' @param df data.frame:  DE results, row.names = gene names
#' @param lengths numerical vector; A list of signature lengths
#' @param sides character vector; A list of signature sides. Possible values are : 
#' up (only up-regulated genes), down (only down-regulated genes) and both (both 
#' up and down-regulated genes).
#' @param remove_regex character vector; A list of regex that defines genes names to remove.
#' if NULL or "none": No genes are removed. 
#' Default = NULL
#' @param gene.selection.methods character vector; signature gene selection method.
#' possible values are "logFC" (select from logFC) or "pval" (select from pval)
#' Default = c("logFC")
#' @param pval_col character: The name of the p adjusted column 
#' @param log2FC_col character: The name of the Log Fold Change column 
#' @param pval_limit numeric: The upper significance limit of the adjusted p-value
#' @param log2FC_limits vector of numeric: The lower and upper limits of the 
#' logFC column to define the up and down regulated genes.
#' 
#' @return list; list of up and down regulated gene. names of the list are 
#'   - up = n up-regulated gene
#'   - down = n down regulated genes
#' 
#' @export
GetSignatureList <- function(df, lengths, sides, remove_regex = NULL, 
                             gene.selection.methods = c("logFC"),
                             pval_col = "padj", log2FC_col = "logFC", 
                             pval_limit = 0.05, log2FC_limits = c(0, 0)){
  
  if (!all(sides %in% SIGNATURE.SIDE)){
    stop("GetSignatureList: sides input not recognized: Possible values are 'up', 'down', 'both'")
  }
  if(!all(gene.selection.methods %in% SIGNATURE.SELECTION.METHODS)){
    stop(paste0("GetSignatureList: methods input not recognized: Possible values are :", paste(SIGNATURE.SELECTION.METHODS, collapse = ", ")))
  }
  
  remove.regex.list <- list()
  if (!is.null(remove_regex)){
    remove_regex <- c(remove_regex)
    remove_regex[remove_regex %in% c("none", "NULL")] <- ""
    for (rm.reg in remove_regex){
      if (rm.reg == ""){
        remove.regex.list[["raw"]] <- "none"
      } else {
        remove.regex.list[[paste0("rm:", rm.reg)]] <- rm.reg
      }
    }
  } else {
    remove.regex.list[["raw"]] <- "none"
  }
  
  signature.list <- list()
  
  for (rm.reg.name in names(remove.regex.list)){ # rm.reg.name <- names(remove.regex.list)[1]
    rm.reg <- remove.regex.list[[rm.reg.name]]
    if (rm.reg == "none" | rm.reg == ""){
      rm.reg <- NULL
    }
    
    for (selection.method in gene.selection.methods){
      signature <- GetSignature(
        df = df, 
        n_genes = max(lengths), 
        pval_col = pval_col, 
        log2FC_col = log2FC_col, 
        pval_limit = pval_limit, 
        log2FC_limits = log2FC_limits,
        gene.selection.method = selection.method,
        remove_regex = rm.reg)
      signature$up <- signature$up[!is.na(signature$up)]
      signature$down <- signature$down[!is.na(signature$down)]
      
      if ("up" %in% sides){
        signature.s <- "up"
        previous.length <- -1
        for (signature.l in lengths){
          signature.key <- paste(rm.reg.name, selection.method, signature.s, signature.l, sep = "_")
          
          signature.tmp <- signature
          signature.tmp$down <- NULL
          signature.tmp$up <- signature.tmp$up[1:signature.l]
          signature.tmp$up <- signature.tmp$up[!is.na(signature.tmp$up)]
          
          if (length(signature.tmp$up) > previous.length){
            signature.list[[signature.key]] <- signature.tmp
            previous.length <- length(signature.tmp$up)
            previous.length <- ifelse(previous.length == 0, -1, previous.length)
          } else {
            break
          }
        }
      } 
      if ("down" %in% sides){
        signature.s <- "down"
        previous.length <- -1
        for (signature.l in lengths){
          signature.key <- paste(rm.reg.name, selection.method, signature.s, signature.l, sep = "_")
          
          signature.tmp <- signature
          signature.tmp$up <- NULL
          signature.tmp$down <- signature.tmp$down[1:signature.l]
          signature.tmp$down <- signature.tmp$down[!is.na(signature.tmp$down)]
          
          if (length(signature.tmp$down) > previous.length){
            signature.list[[signature.key]] <- signature.tmp
            previous.length <- length(signature.tmp$down)
            previous.length <- ifelse(previous.length == 0, -1, previous.length)
          } else {
            break
          }
        }
      } 
      if ("both" %in% sides){
        signature.s <- "both"
        previous.length <- -1
        for (signature.l in lengths){
          signature.key <- paste(rm.reg.name, selection.method, signature.s, signature.l, sep = "_")
          
          signature.tmp <- signature
          signature.tmp$down <- signature.tmp$down[1:signature.l]
          signature.tmp$down <- signature.tmp$down[!is.na(signature.tmp$down)]
          signature.tmp$up <- signature.tmp$up[1:signature.l]
          signature.tmp$up <- signature.tmp$up[!is.na(signature.tmp$up)]
          
          if ((length(signature.tmp$up) + length(signature.tmp$down)) > previous.length){
            signature.list[[signature.key]] <- signature.tmp
            previous.length <- length(signature.tmp$up) + length(signature.tmp$down)
            previous.length <- ifelse(previous.length == 0, -1, previous.length)
          } else {
            break
          }
        }
      }
    }
  }
  
  return(signature.list)
}



#' Signature Cross Validation
#' 
#' Function to train and test a signature score model with a series of parameters
#' This function first compute the differential expression analysis results to create signatures. 
#' The signature score is then computed for the signatures and the function returns 
#' the accuracy predictions. 
#' 
#' @param data.train Seurat; The Seurat training dataset
#' @param  data.test Seurat; The Seurat testing dataset
#' @param y.label character; The name of the meta.data column corresponding to the binary output
#' @param y.sample character; The name of the meta.data column corresponding to the sample information (usefull for DESeq2, EdgeR and Limma)
#' Default = NULL i.e. no sample information
#' @param y.covariates vector of he name of the meta.data column corresponding to the models covariates
#' Default = NULL i.e. no covariates
#' @param DEA.method character; Differential-Expression-Analysis method name. 
#' possible values are possible values in runDEA()
#' Default = "wilcox" from Seurat::FindMarkers()
#' @param signature.selection.method character vector; signatuire gene selection method.
#' possible values are "logFC" (select from logFC) or "pval" (select from pval)
#' Default = c("logFC")
#' @param signature.lengths numerical vector; Signature lengths
#' @param signature.sides character vector; Signature sides. Possible values ar
#' "up" (up-regulated genes), "down" (down-regulated genes) and "both" (both up 
#' and down-regulated genes)
#' Default = c("both")
#' @param signature.rm.regex character vector; A list of regex that defines 
#' genes names to remove. If genes are identified by the regex, they are removed.
#' if NULL or "none": No genes are removed. 
#' Default = NULL
#' @param signature.methods character vector; Signature score methods. 
#' Possible values are methods in GetSignatureScore()
#' @param assay character; The Seurat assay. Default = "RNA"
#' @param slot character; The Seurat slot. Default = "counts"
#' 
#' @return data.frame
#' 
#' @export
SignatureCrossValidation <- function(data.train, y.label, 
                                     data.test = NULL,
                                     y.sample = NULL, 
                                     y.covariates = NULL, 
                                     DEA.method = DEA.METHODS, 
                                     signature.selection.method = c("logFC"),
                                     signature.lengths = c("20"), 
                                     signature.sides = c("both"), 
                                     signature.rm.regex = NULL,
                                     signature.methods = c("AUCell"),
                                     assay = "RNA", slot = "counts", 
                                     save.DEA.folder = NULL){
  
  # Parameters settings and sanity checks:
  DEA.method <- match.arg(DEA.method)
  if (!all(signature.methods %in% SIGNATURE.SCORE.METHODS)){
    stop(paste0("GetSignatureList: signature.methods input not recognized: Possible values are :", paste(SIGNATURE.SCORE.METHODS, collapse = ", ")))
  }
  if (!all(signature.sides %in% SIGNATURE.SIDE)){
    stop("GetSignatureList: sides input not recognized: Possible values are 'up', 'down', 'both'")
  }
  if(!all(signature.selection.method %in% SIGNATURE.SELECTION.METHODS)){
    stop(paste0("GetSignatureList: signature.selection.method input not recognized: Possible values are :", paste(SIGNATURE.SELECTION.METHODS, collapse = ", ")))
  }
  
  # Part 1: Differential Gene expression Analysis
  message("Run Differential Gene expression Analysis")
  
  dea.res <- RunDEA(
    object = data.train, 
    col.DE_group = y.label, 
    col.sample = y.sample, 
    col.covariate = y.covariates, 
    assay = "RNA", slot = "counts", 
    method = DEA.method)
  
  if (!is.null(save.DEA.folder)){
    save.DEA.file = paste0(save.DEA.folder, "DEA_",DEA.method,".rds")
    saveRDS(dea.res, file = save.DEA.file)
  }
  
  
  # Part 2: Create signature
  message("Get Signatures")
  
  signature.list <- GetSignatureList(
    df = dea.res,
    pval_col = "padj", 
    log2FC_col = "logFC", 
    pval_limit = 0.05, 
    log2FC_limits = c(0, 0), 
    gene.selection.methods = signature.selection.method,
    lengths = signature.lengths, 
    sides = signature.sides,
    remove_regex = signature.rm.regex)
  
  # Part 3: Compute signature score
  message("Compute Signature Scores")
  
  # Create empty data-frames containers
  signature.score.train <- data.frame(row.names = colnames(data.train@assays[["RNA"]]@counts))
  if(!is.null(data.test)){
    signature.score.test <- data.frame(row.names = colnames(data.test@assays[["RNA"]]@counts))
  }
  
  # pre-compute ranks
  ranks.matrices.train <- list()
  ranks.matrices.test <- list()
  X.train <- data.train@assays[["RNA"]]@counts
  if (!is.null(data.test)){
    X.test <- data.test@assays[["RNA"]]@counts
  }
  if ("UCell" %in% c(signature.methods)) {
    ranks.matrices.train[["UCell"]] <- UCell::StoreRankings_UCell(X.train)
    if (!is.null(data.test)){
      ranks.matrices.test[["UCell"]] <- UCell::StoreRankings_UCell(X.test)
    }
  }
  if ("AUCell" %in% c(signature.methods)) {
    ranks.matrices.train[["AUCell"]] <- AUCell::AUCell_buildRankings(exprMat = X.train, plotStats = F, verbose = F)
    if (!is.null(data.test)){
      ranks.matrices.test[["AUCell"]] <- AUCell::AUCell_buildRankings(exprMat = X.test, plotStats = F, verbose = F)
    }
  }
  if ("singscore" %in% c(signature.methods)) {
    ranks.matrices.train[["singscore"]] <- singscore::rankGenes(as.matrix(X.train))
    if (!is.null(data.test)){
      ranks.matrices.test[["singscore"]] <- singscore::rankGenes(as.matrix(X.test))
    }
  }
  rm(X.train)
  if (!is.null(data.test)){
    rm(X.test)
  }
  
  for (ss.method in c(signature.methods)){ # ss.method <- c(signature.methods)[1]
    # Get the training signature score
    signature.score.tmp <- GetSignatureScore(
      object = data.train, 
      assay = "RNA", slot = "counts",
      ranks = ranks.matrices.train[[ss.method]],
      signature = signature.list, 
      method = ss.method)
    colnames(signature.score.tmp) <- paste(ss.method, colnames(signature.score.tmp), sep = "_")
    signature.score.train <- cbind(signature.score.train, signature.score.tmp)
    
    # Get the Testing signature score
    if(!is.null(data.test)){
      signature.score.tmp <- GetSignatureScore(
        object = data.test, 
        assay = "RNA", slot = "counts",
        ranks = ranks.matrices.test[[ss.method]],
        signature = signature.list,
        method = ss.method)
      colnames(signature.score.tmp) <- paste(ss.method, colnames(signature.score.tmp), sep = "_")
      signature.score.test <- cbind(signature.score.test, signature.score.tmp)
    }
  }
  
  if (!is.null(save.DEA.folder)){
    save.DEA.file = paste0(save.DEA.folder, "signature.score.train_",DEA.method,".rds")
    saveRDS(signature.score.train, file = save.DEA.file)
    if (!is.null(data.test)){
      save.DEA.file = paste0(save.DEA.folder, "signature.score.test_",DEA.method,".rds")
      saveRDS(signature.score.test, file = save.DEA.file)
    }
  }
  
  # Part 4: Compute accuracies: 
  message("Compute Accuracies")
  
  # Get the output labels 
  y.train <- as.factor(data.train@meta.data[, y.label])
  y.train <- y.train == levels(y.train)[1]
  if(!is.null(data.test)){
    y.test <- as.factor(data.test@meta.data[, y.label])
    y.test <- y.test == levels(y.test)[1]
  }
  
  res.inner <- data.frame()
  for (hyper.col in colnames(signature.score.train)){ # hyper.col <- colnames(signature.score.train)[1]
    
    if (!all(is.na(signature.score.train[,hyper.col]))){
      if (!is.null(data.test)){
        x.test = signature.score.test[,hyper.col]
      } else {
        x.test = NULL
      }
      pred.res <- GetSignaturePrediction(
        x.train = signature.score.train[,hyper.col], 
        x.test = x.test,
        method = "greedy", 
        ground.truth = y.train)
      
      # Get train & test predictions: 
      y.train.pred <- pred.res$y.train.pred
      y.test.pred <- pred.res$y.test.pred
      
      # Get training accuracy:
      res.train.metrics <- GetMetricsBinary(ground_truth = y.train, preds = y.train.pred)
      res.train.metrics.acc <- GetAccuracyMetrics(tp = res.train.metrics$TP, 
                                                  fp = res.train.metrics$FP, 
                                                  tn = res.train.metrics$TN, 
                                                  fn = res.train.metrics$FN)
      
      # Get testing accuracy
      if (!is.null(data.test)){
        res.test.metrics <- GetMetricsBinary(ground_truth = y.test, preds = y.test.pred)
        res.test.metrics.acc <- GetAccuracyMetrics(tp = res.test.metrics$TP, 
                                                   fp = res.test.metrics$FP, 
                                                   tn = res.test.metrics$TN, 
                                                   fn = res.test.metrics$FN)
      } else {
        res.test.metrics <- list(TP = NA, FP = NA, TN = NA, FN = NA)
        res.test.metrics.acc <- list("accuracy" = NA, "mcc" = NA, "F1" = NA, "kappa" = NA) 
      }
    } else {
      res.train.metrics <- list(TP = NA, FP = NA, TN = NA, FN = NA)
      res.train.metrics.acc <- list("accuracy" = NA, "mcc" = NA, "F1" = NA, "kappa" = NA) 
      res.test.metrics <- list(TP = NA, FP = NA, TN = NA, FN = NA)
      res.test.metrics.acc <- list("accuracy" = NA, "mcc" = NA, "F1" = NA, "kappa" = NA) 
    }
    
    # Save results in cv.df.outer.inner:
    # TODO split "_" replace by "::" maybe..
    str.split.res <- strsplit(x =  hyper.col, split = "_")[[1]]
    if (length(str.split.res) > 5){
      warning("What have you done... str.split.res is more than 4 elements... please check")
      cat(hyper.col)
      cat("\n")
    }
    
    # TODO add y.label, y.sample = NULL, y.covariates = NULL, DEA.method = DEA.METHODS, ASSAY AND SLOT
    res.info <- list(
      "DEA.method" = DEA.method,
      "signature.methods" = str.split.res[1],
      "signature.rm.regex" = str.split.res[2],
      "signature.selection.method" = str.split.res[3],
      "signature.side" = str.split.res[4],
      "signature.lengths" = str.split.res[5]
    )
    names(res.train.metrics) <- paste0("train.", names(res.train.metrics))
    names(res.train.metrics.acc) <- paste0("train.", names(res.train.metrics.acc))
    
    names(res.test.metrics) <- paste0("test.", names(res.test.metrics))
    names(res.test.metrics.acc) <- paste0("test.", names(res.test.metrics.acc))
    
    res.info <- c(res.info, res.train.metrics, res.train.metrics.acc, res.test.metrics, res.test.metrics.acc)
    
    res.inner <- rbind(res.inner, res.info)
  }
  
  return(res.inner)
}

#' Get Signature Accuracy
#' 
#' Function that compute accuracy based on a threshold
#' 
#' @param x numerical predictions
#' @param y logical ground truth
#' @param th numerical threhsold to apply on x to get logical predictions
#' 
#' @return accuracy for x for a threshold th
#' 
#' @export
GetSignatureAccuracy <- function(x, y, th){
  pred <-  (x >= th) == y
  return(sum(pred)/length(pred))
}


#' Get Signature Prediction
#' 
#' Function to get the prediction for some numerical score. In brief, a threshold 
#' is applied on a score vector to get binary prediction which are returned by 
#' the function
#' The threshold can be provided by the user (method = "threshold") or the 
#' threshold can be computed in a "greedy" manner by optimizing the accuracy. 
#' So there are two methods: 
#' - "threshold": unbiased method where `x.threshold` must be provided
#' - "greedy": Optimal threshold is computed to optimize maxium accuracy. Here, 
#' 'ground.truth' needs to be provided. 
#' 
#' @param x.train numerical vector; The training score vector
#' @param method character; The threhsold computation method. Possible values 
#' are "greedy" and "threshold".
#' Default = "greedy"
#' @param ground.truth logical vector; The logical ground truth needed for when method = "greedy"
#' @param x.threshold numerical; The threshold input for when method = "threshold"
#' @param x.test numerical vector; The testing score vector. If NULL, return no testing prediction
#' Default = NULL
#' 
#' @return list("y.train.pred" = ..., "y.test.pred" = ...)
#' 
#' @export
GetSignaturePrediction <- function(x.train, 
                                   method = c("greedy", "threshold"), 
                                   ground.truth = NULL,
                                   x.threshold = NULL,
                                   x.test = NULL){
  method <- match.arg(method)
  
  # sanity checks + control inputs
  if (method == "threshold" & is.null(x.threshold)){
    warning("GetSignaturePrediction: Method is 'threshold' and x.threshold is not defined. Continue with x.threshold = 0")
    x.threshold <- 0
  } else if (method == "greedy" & is.null(ground.truth)){
    if (is.null(x.threshold)){
      x.threshold <- 0
    }
    method <- "threshold"
    warning(paste0("GetSignaturePrediction: Method is 'greedy' and ground.truth is not defined. Continue with method = 'threshold' and x.threshold = ", x.threshold))
  }
  
  if (method == "greedy"){
    range_ <- range(x.train, na.rm = T)
    if ((range_[2]-range_[1]) > 0){
      opt.res <- stats::optimize(
        f = GetSignatureAccuracy, 
        interval = range(x.train, na.rm = T), 
        x = x.train, 
        y = ground.truth, 
        maximum = T)
      x.threshold <- opt.res$maximum
    } else {
      x.threshold <- 0
    }
    y.train.pred <- x.train >= x.threshold
  } else if (method == "threshold"){
    y.train.pred <- x.train >= x.threshold
  }
  
  if (!is.null(x.test)){
    y.test.pred <- x.test >= x.threshold
  } else {
    y.test.pred <- NULL
  }
  
  return(list("y.train.pred" = y.train.pred, "y.test.pred" = y.test.pred, "x.threshold" = x.threshold))
}

