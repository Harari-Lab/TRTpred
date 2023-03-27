# Script To run DEA using Seurat as input
# Rémy Pétremand
# references:
# Those functions were inspired by the work of others: 
# - Soneson, C., Robinson, M. 
#   Bias, robustness and scalability in single-cell differential expression analysis. N
#   at Methods 15, 255–261 (2018). DOI: https://doi.org/10.1038/nmeth.4612
#   Github: https://github.com/csoneson/conquer_comparison
# - Squair, J.W., Gautier, M., Kathe, C. et al. 
#   Confronting false discoveries in single-cell differential expression. 
#   Nat Commun 12, 5692 (2021). DOI: https://doi.org/10.1038/s41467-021-25960-2
#   Github: https://github.com/neurorestore/DE-analysis

suppressMessages(require(Seurat))
suppressMessages(require(limma))
suppressMessages(require(edgeR))
suppressMessages(require(DESeq2))


LIMMA.METHODS <- c("limma_voom", "limma_trend")
EDGER.METHODS <- c("edgeR_LRT", "edgeR_QFL")
DESEQ.METHODS <- c("DESeq2_Wald", "DESeq2_LRT")
SEURAT.METHODS <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR")


.n_unique <- function(x){length(unique(x))}
.unique_1st <- function(x){unique(x)[1]}


#' get Pseudo Bulk
#' 
#' Function to compute the pseudobulk information by combining information from single cells
#' 
#' @param data matrix, assay with rows = genes and columns = observations (single-cells)
#' @param colData data.frame; columns = col data, rows = observations (single-cells)
#' @param col.aggregate character; aggregation column in colData
#' 
#' @return list of two elements
#' - data = aggregated data (matrix)
#' - colData = aggreagted colData (data.frame)
getPseudoBulk <- function(data, colData, col.aggregate){
  
  # Aggregate colData
  colData.cols <- colnames(colData)
  colData.cols <- colData.cols[colData.cols != col.aggregate]
  
  colData.n.unique <- colData %>% 
    group_by(across(all_of(col.aggregate))) %>% 
    summarise(across(.cols = colData.cols, .n_unique))
  colData.n.unique <- data.frame(colData.n.unique)
  
  if (any(colData.n.unique[, colData.cols] > 1)){
    warning("getPseudoBulk: Confusion in aggregation of colData")
  }
  
  colData.aggr <- colData %>% 
    group_by(across(all_of(col.aggregate))) %>% 
    summarise(across(.cols = colData.cols, .unique_1st))
  colData.aggr <- data.frame(colData.aggr)
  rownames(colData.aggr) <- colData.aggr[[col.aggregate]]
  colData.aggr <- colData.aggr[, colData.cols]
  
  # Aggregate expression data
  gene.names <- rownames(data)
  sample.names <- colnames(data)
  
  data.df <- data.frame(t(data), check.rows = F, check.names = F)
  data.df$aggregate_col <- colData[[col.aggregate]]
  
  data.aggr <- data.df %>% 
    group_by(aggregate_col) %>% 
    summarise(across(all_of(gene.names), mean))
  data.aggr <- data.frame(data.aggr, check.rows = F, check.names = F)
  rownames(data.aggr) <- data.aggr$aggregate_col
  data.aggr <- data.aggr[, gene.names]
  
  data.aggr <- as.matrix(t(data.aggr))
  
  return(list("data" = data.aggr, "colData" = colData.aggr))
  
}


#' Run Differential-Expression Analysis
#' 
#' Function to run the differential expression analysis. This function calls the 
#' sub-function corresponding to the method asked.
#' 
#' @param object Seurat; The Seurat input
#' @param col.DE_group character; The name of the column in meta.data 
#' corresponding to the differential expression analysis groupe
#' @param assay character; The Seurat assay. Default = NULL i.e. Seurat::DefaultAssay(object = object)
#' @param slot character; The slot of the assay. Default = "data"
#' @param method character; The Differential-Expression Analysis method. 
#' possible values are "wilcox" (default), "bimod", "roc", "t", "negbinom", 
#' "poisson", "LR", "DESeq2_Wald", "DESeq2_LRT", "edgeR_LRT", "edgeR_QFL", 
#' "limma_voom", "limma_trend"
#' @param col.sample character; The name of the column in meta.data 
#' corresponding to the sample information
#' Default = NULL
#' @param col.covariate character vector; The names of the column in meta.data 
#' corresponding to the covariates in the models
#' Default = NULL
#' @param col.aggregate character; Column in mca@meta.data to combine cells into pseudobulk. 
#' If NULL, no pseudobulk integration
#' Default = NULL
#' 
#' @return data.frame. The DE result data.frame
RunDEA <- function(object, col.DE_group, 
                   assay = NULL, slot = "data", 
                   method = c(SEURAT.METHODS, EDGER.METHODS, LIMMA.METHODS, DESEQ.METHODS), 
                   col.sample = NULL, 
                   col.covariate = NULL,
                   col.aggregate = NULL){
  
  method <- match.arg(method) # default = "wilcox"
  
  if (is.null(assay)){
    assay <- Seurat::DefaultAssay(object = object)
  }
  
  if (method %in% SEURAT.METHODS){
    de.res <- RunSeurat(mca = object, 
                         assay = assay, 
                         slot = slot,
                         method = method, 
                         col.DE_group = col.DE_group)
  } else if (method %in% EDGER.METHODS){
    de.res <- RunedgeR(mca = object, 
                       assay = assay, 
                       slot = slot,
                       method = gsub("edgeR_", "", method), # QFL or LRT
                       col.DE_group = col.DE_group, 
                       col.sample = col.sample, 
                       col.covariate = col.covariate,
                       col.aggregate = col.aggregate)
    
  } else if (method %in% DESEQ.METHODS){
    de.res <- RunDESeq2(mca = object, 
                        assay = assay, 
                        slot = slot,
                        method = gsub("DESeq2_", "", method), # Wald or LRT
                        col.DE_group = col.DE_group, 
                        col.sample = col.sample, 
                        col.covariate = col.covariate,
                        col.aggregate = col.aggregate)
    
  } else if (method %in% LIMMA.METHODS){
    de.res <- Runlimma(mca = object, 
                       assay = assay, 
                       slot = slot,
                       method = gsub("limma_", "", method), # voom or trend
                       col.DE_group = col.DE_group, 
                       col.sample = col.sample, 
                       col.covariate = col.covariate,
                       col.aggregate = col.aggregate)
  } else {
    stop("")
  } 
  
  return(de.res)
}


#' Limma DE function implementation
#' 
#' Single-cell implementation of the limma DE
#' 
#' @param mca, Seurat: The Seurat object
#' @param assay, character. The Seurat assay. Default = Seurat::DefaultAssay(object = object)
#' @param slot, character. The Seurat assay's slot. Default = "data"
#' @param method, character: The limma method. Possible value are "voom" (default) 
#' and "trend"
#' @param col.DE_group, character: The colname in mca's meta.data defining the 
#' Differential Expression group. default = "DE_group". cannot be NULL
#' @param col.sample, character: The colname in mca's meta.data defining the 
#' sample information. default = "Patient". If NULL, every sample are in the same group
#' @param col.covariate, vector of character: A vector of colnames in mca's 
#' meta.data defining the covariates to add in the DE design. 
#' If NULL, no covariates are added. Default = c("Patient"). 
#' @param col.aggregate character; Column in mca@meta.data to combine cells into pseudobulk. 
#' If NULL, no pseudobulk integration
#' Default = NULL
#' 
#' @return data.frame. The DE result data.frame
Runlimma <- function(mca, 
                     assay = NULL, slot = "data",
                     method = c("voom", "trend"), 
                     col.DE_group = "DE_group", 
                     col.sample = "Patient", 
                     col.covariate = c("Patient"),
                     col.aggregate = NULL){

  method <- match.arg(method)
  
  if (is.null(assay)){
    assay <- Seurat::DefaultAssay(object = mca)
  } 
  
  if (is.null(col.DE_group)){
    stop("col.DE_group needs to be defined")
  }
  if (is.null(col.sample)){
    mca$col.sample_ <- factor("1", levels = c("1"))
    col.sample <- "col.sample_"
  }
  
  if (!col.DE_group %in% colnames(mca@meta.data)){
    warning("DE_group_col not in metadata columns")
  } 
  if (!col.sample %in% colnames(mca@meta.data)){
    warning("DE_group_col not in metadata columns")
  }
  if (!is.null(col.covariate)){
    for (col_ in col.covariate) {
      if (!col_ %in% colnames(mca@meta.data)){
        warning(paste0("DE_group_col: ", col_, " is not found in metadata columns"))
      }
    }
  }
  
  # Remove possible NA values
  if (any(is.na(mca@meta.data[[col.DE_group]]))){
    mca@meta.data$group.tmp_ <- mca@meta.data[, col.DE_group]
    unique.values <- unique(mca@meta.data[, col.DE_group])
    unique.values <- unique.values[!is.na(unique.values)]
    mca <- subset(mca, subset = group.tmp_ %in% unique.values)
    mca@meta.data$group.tmp_ <- NULL
  }
  
  if (class(mca[[col.DE_group]][, 1]) != "factor"){
    mca[[col.DE_group]] <- factor(mca[[col.DE_group]])
  }
  
  # Define formula with colData columns
  if (!is.null(col.covariate)){
    formula_str <- paste0("~ 0 + ", col.DE_group, " + ", paste(col.covariate, collapse = " + "))
    column.colData <- c(col.DE_group, col.sample, col.covariate)
  } else {
    formula_str <- paste0("~ 0 + ", col.DE_group)
    column.colData <- c(col.DE_group, col.sample)
  }
  
  # Prepare expr.data + colData
  if (is.null(col.aggregate)){
    expr.data <- Seurat::GetAssayData(object = mca, slot = slot, assay = assay)
    
    colData <- mca@meta.data[, unique(column.colData)]
    colnames(colData) <- gsub(".*\\:", "", colnames(colData))
  } else { # col.aggregate <- "TCR_patient"
    
    expr.data <- Seurat::GetAssayData(object = mca, slot = slot, assay = assay)
    
    colData <- mca@meta.data[, unique(c(column.colData, col.aggregate))]
    colnames(colData) <- gsub(".*\\:", "", colnames(colData))
    
    pseudobulk.list <- getPseudoBulk(data = expr.data, colData = colData, col.aggregate = col.aggregate)
    expr.data <- pseudobulk.list$data
    colData <- pseudobulk.list$colData
  }
  
  # Negative values are not allowed
  if (min(expr.data, na.rm = T) < 0){
    expr.data <- expr.data + abs(min(expr.data, na.rm = T))
  }
  
  dge <- edgeR::DGEList(counts = expr.data, group = colData[[col.sample]])
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(formula(formula_str), data = colData)
  colnames(design) <- gsub(col.DE_group, "", colnames(design))
  
  if (method == "trend"){
    y <- new("EList")
    y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
    fit <- limma::lmFit(y, design = design)
  } else {
    vm <- limma::voom(dge, design = design, plot = TRUE)
    fit <- limma::lmFit(vm, design = design)
  }
  
  # TODO: CHECK IF WE NEED duplicateCorrelation
  # dupcor <- duplicateCorrelation(vm, design = design, block = colData$Patient)
  
  contr <- limma::makeContrasts(yes - no, levels = colnames(coef(fit)))
  contr.fit <- limma::contrasts.fit(fit, contr)
  
  trend_bool <- method == "trend"
  contr.fit <- limma::eBayes(contr.fit, trend = trend_bool, robust = trend_bool)
  tt <- limma::topTable(contr.fit, n = Inf, adjust.method = "fdr")
  
  res <- data.frame(
    row.names = rownames(tt),
    genes = rownames(tt),
    logFC = tt$logFC,
    AveExpr = tt$logFC,
    pval = tt$P.Value,
    padj = tt$adj.P.Val,
    t = tt$t,
    B = tt$B)
  
  return(res)
}

#' EdgeR DE function implementation
#' 
#' Single-cell implementation of the EdgeR DE
#' 
#' @param mca, Seurat: The Seurat object
#' @param assay, character. The Seurat assay. Default = Seurat::DefaultAssay(object = object)
#' @param slot, character. The Seurat assay's slot. Default = "data"
#' @param  method, character: The limma method. Possible value are "QFL" (default) 
#' and "LRT"
#' @param col.DE_group, character: The colname in mca's meta.data defining the 
#' Differential Expression group. default = "DE_group". cannot be NULL
#' @param col.sample, character: The colname in mca's meta.data defining the 
#' sample information. default = "Patient". If NULL, every sample are in the same group
#' @param col.covariate, vector of character: A vector of colnames in mca's 
#' meta.data defining the covariates to add in the DE design. 
#' If NULL, no covariates are added. Default = c("Patient"). 
#' @param col.aggregate character; Column in mca@meta.data to combine cells into pseudobulk. 
#' If NULL, no pseudobulk integration
#' Default = NULL
#' 
#' @return data.frame. The DE result data.frame
RunedgeR <- function(mca, 
                     assay = NULL, slot = "data",
                     method = c("QFL", "LRT"), 
                     col.DE_group = "DE_group", 
                     col.sample = "Patient", 
                     col.covariate = c("Patient"), 
                     col.aggregate = NULL){
  
  method <- match.arg(method)
  
  if (is.null(assay)){
    assay <- Seurat::DefaultAssay(object = mca)
  }
  
  if (is.null(col.DE_group)){
    stop("col.DE_group needs to be defined")
  }
  if (is.null(col.sample)){
    mca$col.sample_ <- factor("1", levels = c("1"))
    col.sample <- "col.sample_"
  }
  
  if (!col.DE_group %in% colnames(mca@meta.data)){
    warning("DE_group_col not in metadata columns")
  } 
  if (!col.sample %in% colnames(mca@meta.data)){
    warning("DE_group_col not in metadata columns")
  }
  if (!is.null(col.covariate)){
    for (col_ in col.covariate) {
      if (!col_ %in% colnames(mca@meta.data)){
        warning(paste0("DE_group_col: ", col_, " is not found in metadata columns"))
      }
    }
  }
  
  # Remove possible NA values
  if (any(is.na(mca@meta.data[[col.DE_group]]))){
    mca@meta.data$group.tmp_ <- mca@meta.data[, col.DE_group]
    unique.values <- unique(mca@meta.data[, col.DE_group])
    unique.values <- unique.values[!is.na(unique.values)]
    mca <- subset(mca, subset = group.tmp_ %in% unique.values)
    mca@meta.data$group.tmp_ <- NULL
  }
  
  if (class(mca[[col.DE_group]][, 1]) != "factor"){
    mca[[col.DE_group]] <- factor(mca[[col.DE_group]])
  }
  
  # Remove redundant covariates for obtaining full rank model matrix
  col.2.rm <- c()
  for (.col in unique(c(col.sample, col.covariate))){
    first.level <- levels(mca@meta.data[[col.DE_group]])[1]
    mask.pos <- mca@meta.data[[col.DE_group]] == first.level
    
    positive.val <- as.character(unique(mca@meta.data[mask.pos, ][[.col]]))
    negative.val <- as.character(unique(mca@meta.data[!mask.pos, ][[.col]]))
    if (length(intersect(positive.val, negative.val)) == 0){
      col.2.rm <- c(col.2.rm, .col)
    }
  }
  if (length(col.2.rm) > 0){
    warning("The model matrix is not full rank, so the model cannot be fit as specified.", 
            "So we remove the redundant col.covariate: ", paste(col.2.rm, collapse = ", "))
    col.covariate <- col.covariate[!col.covariate %in% col.2.rm]
    if (length(col.covariate) == 0){ col.covariate <- NULL}
  }
  
  # Define formula with colData columns
  if (!is.null(col.covariate)){
    formula_str <- paste0("~ 0 + ", col.DE_group, " + ", paste(col.covariate, collapse = " + "))
    column.colData <- c(col.DE_group, col.sample, col.covariate)
  } else {
    formula_str <- paste0("~ 0 + ", col.DE_group)
    column.colData <- c(col.DE_group, col.sample)
  }
  
  # Prepare expr.data + colData
  if (is.null(col.aggregate)){
    expr.data <- Seurat::GetAssayData(object = mca, slot = slot, assay = assay)
    
    colData <- mca@meta.data[, unique(column.colData)]
    colnames(colData) <- gsub(".*\\:", "", colnames(colData))
  } else { # col.aggregate <- "TCR_patient"
    
    expr.data <- Seurat::GetAssayData(object = mca, slot = slot, assay = assay)
    
    colData <- mca@meta.data[, unique(c(column.colData, col.aggregate))]
    colnames(colData) <- gsub(".*\\:", "", colnames(colData))
    
    pseudobulk.list <- getPseudoBulk(data = expr.data, colData = colData, col.aggregate = col.aggregate)
    expr.data <- pseudobulk.list$data
    colData <- pseudobulk.list$colData
  }
  
  # Negative values are not allowed
  if (min(expr.data, na.rm = T) < 0){
    expr.data <- expr.data + abs(min(expr.data, na.rm = T))
  }
  
  # Negative values are not allowed
  if (min(expr.data, na.rm = T) < 0){
    expr.data <- expr.data + abs(min(expr.data, na.rm = T))
  }
  
  # Drop Levels that are not used:
  for (col.2.treat in c(col.DE_group, col.covariate)){
    colData[[col.2.treat]] <- factor(colData[[col.2.treat]])
    colData[[col.2.treat]] <- droplevels(colData[[col.2.treat]])
  }

  dge <- edgeR::DGEList(counts = expr.data, group = colData[[col.sample]])
  
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(formula(formula_str), data = colData)
  colnames(design) <- gsub(col.DE_group, "", colnames(design))
  
  dge <- edgeR::estimateDisp(dge, design = design)
  
  if (method == "LRT"){
    fit <- edgeR::glmFit(dge, design = design)
    contr <- limma::makeContrasts(TR = yes - no, levels = colnames(coef(fit)))
    fit <- edgeR::glmLRT(fit, contrast = contr[, "TR"])
    tt <- edgeR::topTags(fit, n = Inf, adjust.method = "fdr")
    
    res <- data.frame(
      row.names = rownames(tt$table),
      genes = rownames(tt$table),
      logFC = tt$table$logFC,
      logCPM = tt$table$logCPM,
      padj = tt$table$PValue,
      FDR = tt$table$FDR,
      LR = tt$table$LR)
    
  } else {
    fit <- edgeR::glmQLFit(dge, design = design)
    contr <- limma::makeContrasts(TR = yes - no, levels = colnames(coef(fit)))
    fit <- edgeR::glmQLFTest(fit, contrast = contr[, "TR"])
    tt <- edgeR::topTags(fit, n = Inf, adjust.method = "fdr")
    
    res <- data.frame(
      row.names = rownames(tt$table),
      genes = rownames(tt$table),
      logFC = tt$table$logFC,
      logCPM = tt$table$logCPM,
      padj = tt$table$PValue,
      FDR = tt$table$FDR,
      F = tt$table$F)
    
  }

  return(res) 
}

#' DESeq2 DE function implementation
#' 
#' Single-cell implementation of the DESeq2 DE
#' 
#' @param mca, Seurat: The Seurat object
#' @param assay, character. The Seurat assay. Default = Seurat::DefaultAssay(object = object)
#' @param slot, character. The Seurat assay's slot. Default = "data"
#' @param  method, character: The limma method. Possible value are "Wald" (default) 
#' and "LRT"
#' @param col.DE_group, character: The colname in mca's meta.data defining the 
#' Differential Expression group. default = "DE_group". cannot be NULL
#' @param col.sample, character: The colname in mca's meta.data defining the 
#' sample information. default = "Patient". If NULL, every sample are in the same group
#' @param col.covariate, vector of character: A vector of colnames in mca's 
#' meta.data defining the covariates to add in the DE design. 
#' If NULL, no covariates are added. Default = c("Patient").
#' @param col.aggregate character; Column in mca@meta.data to combine cells into pseudobulk. 
#' If NULL, no pseudobulk integration
#' Default = NULL
#' 
#' @return data.frame. The DE result data.frame
RunDESeq2 <- function(mca, 
                      assay = NULL, slot = "data",
                      method = c("Wald", "LRT"), 
                      col.DE_group = "DE_group", 
                      col.sample = "Patient", 
                      col.covariate = c("Patient"),
                      col.aggregate = NULL){
  
  method <- match.arg(method)
  
  if (is.null(assay)){
    assay <- Seurat::DefaultAssay(object = mca)
  }
  
  if (is.null(col.DE_group)){
    stop("col.DE_group needs to be defined")
  }
  if (is.null(col.sample)){
    mca$col.sample_ <- factor("1", levels = c("1"))
    col.sample <- "col.sample_"
  }
  
  if (!col.DE_group %in% colnames(mca@meta.data)){
    warning("DE_group_col not in metadata columns")
  } 
  if (!col.sample %in% colnames(mca@meta.data)){
    warning("DE_group_col not in metadata columns")
  }
  if (!is.null(col.covariate)){
    for (col_ in col.covariate) {
      if (!col_ %in% colnames(mca@meta.data)){
        warning(paste0("DE_group_col: ", col_, " is not found in metadata columns"))
      }
    }
  }
  
  # Remove possible NA values
  if (any(is.na(mca@meta.data[[col.DE_group]]))){
    mca@meta.data$group.tmp_ <- mca@meta.data[, col.DE_group]
    unique.values <- unique(mca@meta.data[, col.DE_group])
    unique.values <- unique.values[!is.na(unique.values)]
    mca <- subset(mca, subset = group.tmp_ %in% unique.values)
    mca@meta.data$group.tmp_ <- NULL
  }
  
  
  if (class(mca[[col.DE_group]][, 1]) != "factor"){
    mca[[col.DE_group]] <- factor(mca[[col.DE_group]])
  }
  
  # Remove redundant covariates for obtaining full rank model matrix
  col.2.rm <- c()
  for (.col in unique(c(col.sample, col.covariate))){
    first.level <- levels(mca@meta.data[[col.DE_group]])[1]
    mask.pos <- mca@meta.data[[col.DE_group]] == first.level
    
    positive.val <- as.character(unique(mca@meta.data[mask.pos, ][[.col]]))
    negative.val <- as.character(unique(mca@meta.data[!mask.pos, ][[.col]]))
    if (length(intersect(positive.val, negative.val)) == 0){
      col.2.rm <- c(col.2.rm, .col)
    }
  }
  if (length(col.2.rm) > 0){
    warning("The model matrix is not full rank, so the model cannot be fit as specified.", 
            "So we remove the redundant col.covariate: ", paste(col.2.rm, collapse = ", "))
    col.covariate <- col.covariate[!col.covariate %in% col.2.rm]
    if (length(col.covariate) == 0){ col.covariate <- NULL}
  }
  
  # Define design with colData columns 
  if (!is.null(col.covariate)){
    formula_str <- paste0("~ ", col.DE_group, " + ", paste(col.covariate, collapse = " + "))
    column.colData <- c(col.DE_group, col.sample, col.covariate)
  } else {
    formula_str <- paste0("~ ", col.DE_group)
    column.colData <- c(col.DE_group, col.sample)
  }
  design <- formula(formula_str)
  
  # Prepare expr.data + colData
  if (is.null(col.aggregate)){
    expr.data <- Seurat::GetAssayData(object = mca, slot = slot, assay = assay)
    
    colData <- mca@meta.data[, unique(column.colData)]
    colnames(colData) <- gsub(".*\\:", "", colnames(colData))
  } else { # col.aggregate <- "TCR_patient"
    
    expr.data <- Seurat::GetAssayData(object = mca, slot = slot, assay = assay)
    
    colData <- mca@meta.data[, unique(c(column.colData, col.aggregate))]
    colnames(colData) <- gsub(".*\\:", "", colnames(colData))
    
    pseudobulk.list <- getPseudoBulk(data = expr.data, colData = colData, col.aggregate = col.aggregate)
    expr.data <- pseudobulk.list$data
    colData <- pseudobulk.list$colData
  }
  
  # Negative values are not allowed
  if (min(expr.data, na.rm = T) < 0){
    expr.data <- expr.data + abs(min(expr.data, na.rm = T))
  }
  
  # round the expr.data:
  expr.data <- round(expr.data, digits = 0)
  
  # m1 <- model.matrix(object = gsub("+", "*", design), data = colData)
  # all.zero <- apply(m1, 2, function(x) all(x==0))
  
  dds = DESeq2::DESeqDataSetFromMatrix(countData = expr.data,
                               colData = colData, 
                               design = design) # design = ~ design
  
  if (method == "LRT"){
    tryCatch(
      expr = {
        dds.res <- DESeq2::DESeq(object = dds, 
                                 test = "LRT", 
                                 fitType = "local", # c("parametric", "local", "mean", "glmGamPoi") 
                                 sfType = "poscounts", # c("ratio", "poscounts", "iterate")
                                 betaPrior = FALSE,  # TRUE or FALSE,
                                 reduced =  ~ 1)
      },
      error = function(cond){
        message("Error in DESeq. Original Message is: ")
        message(cond, "\n")
        message("Continue with suggested solution")
        
        dds <- DESeq2::estimateSizeFactors(object = dds)
        dds <- DESeq2::estimateDispersionsGeneEst(object = dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds.res <- DESeq2::nbinomLRT(object = dds, reduced =  ~ 1)
      }
    )
  } else {
    tryCatch(
      expr = {
        dds.res <- DESeq2::DESeq(object = dds, 
                                  test = "Wald", 
                                  fitType = "local", # c("parametric", "local", "mean", "glmGamPoi") 
                                  sfType = "poscounts", # c("ratio", "poscounts", "iterate")
                                  betaPrior = FALSE) # TRUE or FALSE
      },
      error = function(cond){
        message("Error in DESeq. Original Message is: ")
        message(cond, "\n")
        message("Continue with suggested solution")
        
        dds <- DESeq2::estimateSizeFactors(object = dds, )
        dds <- DESeq2::estimateDispersionsGeneEst(object = dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds.res <- DESeq2::nbinomWaldTest(object = dds, betaPrior = FALSE)
      }
    )
  }
  
  # note on the parmater selection above: 
  # if fitType = "parametric": I get the warning: 
  #   fitType='parametric', but the dispersion trend was not well captured by the
  #   function: y = a/x + b, and a local regression fit was automatically substituted.
  res = DESeq2::results(dds.res, contrast = c(col.DE_group, levels(mca@meta.data[[col.DE_group]])))

  res.df <- data.frame(
      row.names = rownames(res),
      genes = rownames(res),
      logFC = res$log2FoldChange,
      lfcSE = res$lfcSE,
      pval = res$pvalue,
      padj = res$padj,
      stat = res$stat,
      baseMean = res$baseMean)
  
  return(res.df)
}

#' Seurat DE function implementation
#' 
#' Implementation of the Seurat DE
#' 
#' @param mca, Seurat: The Seurat object
#' @param assay, character. The Seurat assay. Default = Seurat::DefaultAssay(object = object)
#' @param slot, character. The Seurat assay's slot. Default = "data"
#' @param  method, character: The limma method. Possible value are "wilcox" 
#' (default), "bimod", "roc", "t", "negbinom", "poisson" & "LR". The other
#' possible Seurat::FindMarkers methods are not implemented here. 
#' @param col.DE_group, character: The colname in mca's meta.data defining the 
#' Differential Expression group. default = "DE_group". cannot be NULL
#' 
#' @return data.frame. The DE result data.frame
RunSeurat <- function(mca, 
                      assay = NULL, slot = "data",
                      method = SEURAT.METHODS, 
                      col.DE_group = "DE_group"){
  
  method <- match.arg(method)
  
  if (is.null(assay)){
    assay <- Seurat::GetAssayData(object = mca, slot = slot, assay = assay)
  }
  
  if (is.null(col.DE_group)){
    stop("col.DE_group needs to be defined")
  }
  if (!col.DE_group %in% colnames(mca@meta.data)){
    warning("DE_group_col not in metadata columns")
  } 
  
  if (class(mca[[col.DE_group]][, 1]) != "factor"){
    mca[[col.DE_group]] <- factor(mca[[col.DE_group]])
  }
  
  res <- Seurat::FindMarkers(
    object = mca, 
    assay = assay, 
    slot = slot, 
    test.use = method,
    group.by = col.DE_group, 
    ident.1 = levels(mca@meta.data[[col.DE_group]])[1],
    ident.2 = levels(mca@meta.data[[col.DE_group]])[2],
    min.pct = -Inf, 
    only.pos = F, 
    logfc.threshold = -Inf, 
    min.cells.group = 0)
  
  res.df <- data.frame(
    row.names = rownames(res),
    genes = rownames(res),
    logFC = res$avg_log2FC,
    pval = res$p_val,
    padj = res$p_val_adj,
    pct.1 = res$pct.1,
    pct.2 = res$pct.2)
  
  return(res.df)
}

