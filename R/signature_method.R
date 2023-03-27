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

DEFAULT.SIGNATURE.HYPERPARAMS <- list("signature.lengths" = c(20), 
                                      "signature.side" = c("both"),
                                      "signature.rm.regex" = "none",
                                      "signature.methods" = c("AUCell"),
                                      "signature.selection.method" = c("logFC"))

EVALUATION.METRICS <- c("mcc", "accuracy", "F1", "kappa", "auc", "sensitivity", "specificity", "PPV", "NPV")

#' Get signature prediction from a given model
#' 
#' Function to get the signature model prediction when given the model. 
#' 
#' @param x Seurat; The Seurat input data
#' @param path.folder character; Folder to retrieve the model information.
#' @param DEA.file.name character; The DEA result file name
#' @param hyperparameters list; The signature hyperparmeters
#' Elements of list are 
#'   - "signature.lengths": Vector of signature lengths
#'   - "signature.side": Vector of siganture sides ("up", "down" or "both")
#'   - "signature.rm.regex": Vector of rm.regex (regex to remove genes from signature)
#'   - "signature.methods": Vector of signature score methods
#'   - "signature.selection.method": How to select genes from DEA res using the log-Fold-Change ("logFC") or p-values ("pval")
#' @param score.threshold numerical; The signature score threshold to binarize the score into positive and negative predictions
#' 
#' @return data.frame; Rows = barcode names in x. Cols = "score" and "pred"
#' 
#' @export
GetPredictionSignature <- function(x, path.folder, 
                                   DEA.file.name = "dea.res.rds", 
                                   hyperparameters = DEFAULT.SIGNATURE.HYPERPARAMS, 
                                   score.threshold = 0,
                                   score.to.scale = T, 
                                   score.scale.params = NULL){
  
  # 0. Check hyperparameters:
  for (name_ in names(DEFAULT.SIGNATURE.HYPERPARAMS)){
    if (!(name_ %in% names(hyperparameters))){
      hyperparameters[[name_]] <- DEFAULT.SIGNATURE.HYPERPARAMS[[name_]]
    }
  }
  
  # 1. Get DEA data
  DEA.path <- paste0(path.folder, DEA.file.name)
  if (file.exists(DEA.path)){
    dea.res <- readRDS(file = DEA.path)
  } else {
    stop("GetPredictionSignature: DEA file path does not exist. Please check : ", DEA.path)
  }
  
  # 2. Get signature
  signature.list <- GetSignatureList(
    df = dea.res,
    pval_col = "padj", 
    log2FC_col = "logFC", 
    pval_limit = 0.05, 
    log2FC_limits = c(0, 0), 
    gene.selection.methods = hyperparameters$signature.selection.method,
    lengths = hyperparameters$signature.lengths, 
    sides = hyperparameters$signature.side,
    remove_regex = hyperparameters$signature.rm.regex)
  
  # 3. Get signature score
  signature.score <- GetSignatureScore(
    object = x, 
    assay = "RNA", slot = "data",
    signature = signature.list, 
    method = hyperparameters$signature.methods, 
    to.scale = score.to.scale, 
    scale.params = score.scale.params)
  signature.score <- signature.score$score
  
  # 4. Get the signature score
  pred.res <- BinarizeScorePrediction(
    x.train = signature.score, 
    x.test = NULL,
    method = "threshold", 
    x.threshold = score.threshold)
  
  # 5. Get return data.frame
  signature.pred <- cbind(signature.score, pred.res$y.train.pred)
  colnames(signature.pred) <- c("score", "pred")
  
  return(signature.pred)
}

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
#' @param DEA.data data.frame or matrix; The input matrix for the DEA. 
#' Rows = samples, Columns = features. If NULL, the DEA.data is x
#' Default = NULL
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
#' @param col.aggregate character; Column in data.train@meta.data to combine cells into pseudobulk. 
#' If NULL, no pseudobulk integration
#' Default = NULL
#' @param folds.save.folder character; Folder to save the model information. 
#' If null doesn't save anything. If not null, it will save the following data into the folder: 
#' - dea.res.rds: DEA results as in RunDEA()
#' - signature.score.train.rds: signature.score.train data.frame (rows = cell name, col = score for different hyperparmaters)
#' - signature.score.test.rds: signature.score.test data.frame (rows = cell name, col = score for different hyperparmaters)
#' - y.train.pred.rds: y.train.pred data.frame (rows = cell name, col = prediction)
#' - y.test.pred.rds: y.test.pred data.frame (rows = cell name, col = prediction)
#' default = NULL
#' 
#' @return data.frame
#' 
#' @export
SignatureCrossValidation <- function(data.train, y.label, 
                                     data.test = NULL,
                                     y.sample = NULL, 
                                     y.covariates = NULL, 
                                     DEA.method = DEA.METHODS, 
                                     DEA.data = NULL,
                                     signature.selection.method = c("logFC"),
                                     signature.lengths = c("20"), 
                                     signature.sides = c("both"), 
                                     signature.rm.regex = NULL,
                                     signature.methods = c("AUCell"),
                                     assay = "RNA", slot = "counts", 
                                     col.aggregate = NULL, 
                                     sample.weights = NULL,
                                     folds.save.folder = NULL, 
                                     save.model = F, save.pred = F, save.score = F){
  
  # Sanity check on the saving of data
  if (is.null(folds.save.folder)){
    if (any(c(save.model, save.pred, save.score))){
      warning("folds.save.folder is NULL but one or more of `save.model`, `save.pred`, `save.score` is set as TRUE: Nothing will be saved")
      save.model = F;save.pred = F;save.score = F
    }
  } else {
    if (dir.exists(folds.save.folder)){
      if (any(!c(save.model, save.pred, save.score))){
        warning("folds.save.folder exists but nothing will be saved since `save.model`, `save.pred`, `save.score` are set as FALSE")
      }
    }
  }
  
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
  
  # ----------------------------------------------------------------------------
  # Part 1: Differential Gene expression Analysis
  # ----------------------------------------------------------------------------
  message("Run Differential Gene expression Analysis")
  
  if (is.null(DEA.data)){
    DEA.data <- data.train
  }
  
  dea.res <- RunDEA(
    object = DEA.data, 
    col.DE_group = y.label, 
    col.sample = y.sample, 
    col.covariate = y.covariates, 
    assay = assay, slot = slot, 
    method = DEA.method, 
    col.aggregate = col.aggregate)
  
  if (save.model){
    saveRDS(dea.res, file = paste0(folds.save.folder, "dea.res.rds"))
  }
  
  # ----------------------------------------------------------------------------
  # Part 2: Create signature + Get Scores
  # ----------------------------------------------------------------------------
  message("Get Signatures + Scores")
  
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
  
  # Create empty data-frames containers
  signature.score.train <- data.frame(row.names = colnames(data.train@assays[["RNA"]]@counts))
  if(!is.null(data.test)){
    signature.score.test <- data.frame(row.names = colnames(data.test@assays[["RNA"]]@counts))
  }
  scale.params.df <- data.frame()
  
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
      assay = assay, slot = slot,
      ranks = ranks.matrices.train[[ss.method]],
      signature = signature.list, 
      method = ss.method, 
      to.scale = T, 
      scale.params = NULL)
    
    scale.params.train <- list(
      "mean" = signature.score.tmp$scale.mean,
      "sd" = signature.score.tmp$scale.sd)
    
    signature.score.tmp <- signature.score.tmp$score
    colnames(signature.score.tmp) <- paste(ss.method, colnames(signature.score.tmp), sep = "_")
    signature.score.train <- cbind(signature.score.train, signature.score.tmp)
    
    # Get the Testing signature score
    if(!is.null(data.test)){
      signature.score.tmp <- GetSignatureScore(
        object = data.test, 
        assay = assay, slot = slot,
        ranks = ranks.matrices.test[[ss.method]],
        signature = signature.list,
        method = ss.method, 
        to.scale = T,
        scale.params = scale.params.train)
      signature.score.tmp <- signature.score.tmp$score
      colnames(signature.score.tmp) <- paste(ss.method, colnames(signature.score.tmp), sep = "_")
      signature.score.test <- cbind(signature.score.test, signature.score.tmp)
    }
    
    # Save the scale parameters:
    if (is.null(scale.params.train$mean)){
      scale.params.train$mean <- replicate(n = ncol(signature.score.tmp), expr = as.numeric(NA))
    }
    if (is.null(scale.params.train$sd)){
      scale.params.train$sd <- replicate(n = ncol(signature.score.tmp), expr = as.numeric(NA))
    }
    scale.params.df.tmp <- data.frame(
      hyperparams = colnames(signature.score.tmp),
      mean = scale.params.train$mean,
      sd = scale.params.train$sd)
    
    scale.params.df <- rbind(scale.params.df, scale.params.df.tmp)

  }
  
  rownames(scale.params.df) <- scale.params.df$hyperparams
  
  # Save Score - WARNING VERY SPACE HEAVY
  if (save.score){
    saveRDS(signature.score.train, file =  paste0(folds.save.folder, "y_train_score.rds")) # signature.score.train.rds
    saveRDS(scale.params.df, file =  paste0(folds.save.folder, "scale_params_df.rds"))
    if (!is.null(data.test)){
      saveRDS(signature.score.test, file = paste0(folds.save.folder, "y_test_score.rds")) # signature.score.test.rds
    }
  }
  
  # ----------------------------------------------------------------------------
  # Part 3: Compute performance: 
  # ----------------------------------------------------------------------------
  message("Compute Performance")
  
  # Get the output labels 
  y.train <- as.factor(data.train@meta.data[, y.label])
  y.train <- y.train == levels(y.train)[1]
  if(!is.null(data.test)){
    y.test <- as.factor(data.test@meta.data[, y.label])
    y.test <- y.test == levels(y.test)[1]
  }
  
  res.inner <- data.frame()
  if (save.pred){
    y.train.pred.df <- data.frame(row.names = colnames(data.train))
    y.test.pred.df <- data.frame(row.names = colnames(data.test))
  }
  
  for (hyper.col in colnames(signature.score.train)){ # hyper.col <- colnames(signature.score.train)[1]
    x.threshold <- NA
    if (!all(is.na(signature.score.train[,hyper.col]))){
      
      x.train <- signature.score.train[,hyper.col]
      if (!is.null(data.test)){
        x.test = signature.score.test[,hyper.col]
      } else {
        x.test = NULL
      }
      pred.res <- BinarizeScorePrediction(
        x.train = x.train, 
        x.test = x.test,
        method = "greedy", 
        weights = sample.weights[rownames(signature.score.train), 1],
        ground.truth = y.train)
      
      # Get train & test predictions: 
      y.train.pred <- pred.res$y.train.pred
      y.test.pred <- pred.res$y.test.pred
      x.threshold <- pred.res$x.threshold

      # Get training accuracy:
      res.train.metrics <- GetAccuracyMetrics(ground_truth = y.train, preds = y.train.pred, scores = x.train, metrics = EVALUATION.METRICS)
      
      # Get testing accuracy
      if (!is.null(data.test)){
        res.test.metrics <- GetAccuracyMetrics(ground_truth = y.test, preds = y.test.pred, scores = x.test, metrics = EVALUATION.METRICS)
      } else {
        res.test.metrics <- replicate(n = (4 + length(EVALUATION.METRICS)), NA)
        names(res.test.metrics) <- c("TP", "FP", "TN", "FN", EVALUATION.METRICS)
      }
      
      if (save.pred){
        y.train.pred.df.tmp <- data.frame(y.train.pred, row.names = colnames(data.train))
        colnames(y.train.pred.df.tmp) <- hyper.col
        y.train.pred.df <- cbind(y.train.pred.df, y.train.pred.df.tmp)
        
        if (!is.null(data.test)){
          y.test.pred.df.tmp <- data.frame(y.test.pred, row.names = colnames(data.test))
          colnames(y.test.pred.df.tmp) <- hyper.col
          y.test.pred.df <- cbind(y.test.pred.df, y.test.pred.df.tmp)
        }
      }
      
    } else {
      res.train.metrics <- replicate(n = (4 + length(EVALUATION.METRICS)), NA)
      names(res.train.metrics) <- c("TP", "FP", "TN", "FN", EVALUATION.METRICS)
      res.test.metrics <- replicate(n = (4 + length(EVALUATION.METRICS)), NA)
      names(res.test.metrics) <- c("TP", "FP", "TN", "FN", EVALUATION.METRICS)
    }
    
    
    
    # Save results in cv.df.outer.inner:
    # TODO split "_" replace by "::" maybe..
    str.split.res <- strsplit(x =  hyper.col, split = "_")[[1]]
    if (length(str.split.res) > 5){
      warning("What have you done... str.split.res is more than 4 elements... please check")
      cat(hyper.col)
      cat("\n")
    }
    
    res.info <- list(
      "DEA.method" = DEA.method,
      "signature.methods" = str.split.res[1],
      "signature.rm.regex" = str.split.res[2],
      "signature.selection.method" = str.split.res[3],
      "signature.side" = str.split.res[4],
      "signature.lengths" = str.split.res[5],
      "signature.x.threshold" = x.threshold,
      "scale.mean" = scale.params.df[hyper.col,]$mean,
      "scale.sd" = scale.params.df[hyper.col,]$sd
    )
    names(res.train.metrics) <- paste0("train.", names(res.train.metrics))

    names(res.test.metrics) <- paste0("test.", names(res.test.metrics))

    res.info <- c(res.info, res.train.metrics, res.test.metrics)
    
    res.inner <- rbind(res.inner, res.info)
  }
  
  if (save.pred){
    saveRDS(y.train.pred.df, file =  paste0(folds.save.folder, "y_train_pred.rds"))
    if (!is.null(data.test)){
      saveRDS(y.test.pred.df, file = paste0(folds.save.folder, "y_test_pred.rds"))
    }
  }
  
  return(res.inner)
}

