# Script To compute signature score
# Rémy Pétremand

suppressMessages(require(Seurat))

suppressMessages(require(AUCell))
suppressMessages(require(UCell))
suppressMessages(require(singscore))

suppressMessages(require(GSEABase))

SIGNATURE.SCORE.METHODS <- c("average", "UCell", "AUCell", "singscore")


#' Get signature score
#' 
#' Function to get signature score from different methods
#' 
#' @param object Seurat; Seurat single-cell object
#' @param signature list: signature as defined in GetSignature()
#' @param assay character; Assay in `object` to pull the slot data from. 
#' default = "RNA"
#' @param slot character; Slot to pull data from the `assay`. 
#' default = "data"
#' @param method character: The signature score method. 
#' possible methods are "average" (default), "UCell", "AUCell", "singscore"
#' @param scale logical: Weather to scale signature score or not.
#' 
#' @return data.frame; The signature score data.frame. Row = cells and columns = names(signature)
#' 
#' @export
GetSignatureScore <- function(object, signature, ranks = NULL,
                    assay = "RNA", slot = "data",
                    method = SIGNATURE.SCORE.METHODS,
                    scale = F){
  
  
  method <- match.arg(method)
  
  for (signature.key in names(signature)){
    if (!all(names(signature[[signature.key]]) %in% c("up", "down"))){
      print("error")
    }
  }
  
  expr.data <- Seurat::GetAssayData(object = object, slot = slot, assay = assay)
  
  switch(method, 
         average = {
           res <- RunSigScoreAverage(X = expr.data, signature = signature, scale = scale)
         },
         UCell = {
           res <- RunSigScoreUCell(X = expr.data, signature = signature, ranks = ranks, scale = scale)
         },
         AUCell = {
           res <- RunSigScoreAUCell(X = expr.data, signature = signature, ranks = ranks, scale = scale)
         },
         singscore = {
           res <- RunSigScoreSingscore(X = expr.data, signature = signature, ranks = ranks, scale = scale)
         })
  
  return(res)
}


#' Run the average signature score
#' 
#' Function to run the average signature score
#' 
#' @param X matrix; The matrix. Columns = cells, rows = some features
#' @param signature list: signature as defined in GetSignature()
#' @param scale logical: Weather to scale signature score or not.
#' 
#' @return data.frame; The signature score data.frame. Row = cells and columns = names(signature)
#' 
#' @export
RunSigScoreAverage <- function(X, signature, scale = F){
  
  res <- data.frame(row.names = colnames(X))
  for (signature.key in names(signature)){
    s_up <- signature[[signature.key]]$up
    if (length(s_up) == 0){s_up <- NULL}
    s_down <- signature[[signature.key]]$down
    if (length(s_down) == 0){s_down <- NULL}
    
    if (!is.null(s_up)){
      s_up <- s_up[!is.na(s_up)]
      if (length(s_up) > 1){
        up_ss <- colSums(X[s_up, ])/length(s_up)
      } else {
        up_ss <- X[s_up, ]
        warning(paste0("RunSigScoreAverage: strange behaviour: s_up = ", paste(s_up, collapse = ", ")))
      }
    } else {
      up_ss <- replicate(n = dim(X)[2], expr = 0)
    } 
    
    if (!is.null(s_down)){
      s_down <- s_down[!is.na(s_down)]
      if (length(s_down) > 1){
        down_ss <- colSums(X[s_down, ])/length(s_down)
      } else {
        down_ss <- X[s_down, ]
        warning(paste0("RunSigScoreAverage: strange behaviour: s_down = ", paste(s_up, collapse = ", ")))
      }
    } else {
      down_ss <- replicate(n = dim(X)[2], expr = 0)
    } 
    
    res.tmp <- up_ss - down_ss
    res.tmp <- data.frame(res.tmp)
    colnames(res.tmp) <- signature.key
    
    res <- cbind(res, res.tmp)
  }
  
  # res$signature_score_scaled <- scale(res$signature_score, center = TRUE, scale = TRUE)
  
  return(res)
}

#' Run the UCell signature score
#' 
#' Function to run the UCell signature score
#' 
#' @param X matrix; The matrix. Columns = cells, rows = some features
#' @param ranks rank-matrix as defined by UCell::StoreRankings_UCell(X). 
#' If provided, the function does not compute the rankings from X
#' @param signature list: signature as defined in GetSignature()
#' @param scale logical: Weather to scale signature score or not.
#' 
#' @return data.frame; The signature score data.frame. Row = cells and columns = names(signature)
#' 
#' @export
RunSigScoreUCell <- function(X, signature, ranks = NULL, scale = F){ 
  
  if (is.null(ranks)){
    ranks <- UCell::StoreRankings_UCell(X)
  }
  
  signature.input <- list()
  for (signature.key in names(signature)){
    signature.input[[signature.key]] <- c()
    if ("up" %in% names(signature[[signature.key]])){
      if (length(signature[[signature.key]]$up) > 0){
        signature.input[[signature.key]] <- c(
          signature.input[[signature.key]], 
          paste0(signature[[signature.key]]$up, "+"))
      } else {
        signature.input[[signature.key]] <- c(
          signature.input[[signature.key]], as.character())
      }
    }
    if ("down" %in% names(signature[[signature.key]])){
      if (length(signature[[signature.key]]$down) > 0){
        signature.input[[signature.key]] <- c(
          signature.input[[signature.key]], 
          paste0(signature[[signature.key]]$down, "-"))
      } else {
        signature.input[[signature.key]] <- c(
          signature.input[[signature.key]], as.character())
      }
    }
  }
  
  missing.signatures <- c()
  for (siganture.key in names(signature.input)){
    if (length(signature.input[[siganture.key]]) == 0){
      signature.input[[siganture.key]] <- NULL
      missing.signatures <- c(missing.signatures, siganture.key)
    }
  }
  
  res <- data.frame(row.names = colnames(X))
  if (length(signature.input) > 0){
    res <- UCell::ScoreSignatures_UCell(features = signature.input, # list(res = s_up_down), 
                                        precalc.ranks = ranks, w_neg = 1)
  } 
  if (length(missing.signatures) > 0){
    for (col_ in missing.signatures){
      res[[col_]] <- 0
    }
  }
  
  colnames(res) <- gsub("_UCell$", "", colnames(res))
  
  return(res)
}

#' Run the AUCell signature score
#' 
#' Function to run the AUCell signature score
#' 
#' @param X matrix; The matrix. Columns = cells, rows = some features
#' @param ranks rank-matrix as defined by AUCell::AUCell_buildRankings(). 
#' If provided, the function does not compute the rankings from X
#' @param signature list: signature as defined in GetSignature()
#' @param scale logical: Weather to scale signature score or not.
#' 
#' @return data.frame; The signature score data.frame. Row = cells and columns = names(signature)
#' 
#' @export
RunSigScoreAUCell <- function(X, signature, ranks = NULL, scale = F){ 
  
  if (is.null(ranks)){
    ranks <- AUCell::AUCell_buildRankings(exprMat = X, plotStats = F, verbose = F)
  }
  
  signature.input <- list()
  for (signature.key in names(signature)){
    signature.key.new.up <-   paste(signature.key, "up",   sep="_")
    signature.key.new.down <- paste(signature.key, "down", sep="_")
    
    if ("up" %in% names(signature[[signature.key]])){
      signature.input[[signature.key.new.up]] <-   signature[[signature.key]]$up
    }
    if ("down" %in% names(signature[[signature.key]])){
      signature.input[[signature.key.new.down]] <- signature[[signature.key]]$down
    }
  }
  
  for (siganture.key in names(signature.input)){
    if (length(signature.input[[siganture.key]]) == 0){
      signature.input[[siganture.key]] <- NULL
    }
  }
  
  res <- data.frame(row.names = colnames(X))
  if (length(signature.input) > 0){
    res <- AUCell::AUCell_calcAUC(geneSets = signature.input, 
                                  rankings = ranks, 
                                  normAUC = T, verbose = F)
    
    
    res <- data.frame(t(AUCell::getAUC(res)))
  } 
  
  for (signature.key in names(signature)){
    signature.key.new.up <-   paste(signature.key, "up",   sep="_")
    signature.key.new.down <- paste(signature.key, "down", sep="_")
    
    if (!signature.key.new.up %in% colnames(res)){
      res[[signature.key.new.up]] <- 0
    }
    if (!signature.key.new.down %in% colnames(res)){
      res[[signature.key.new.down]] <- 0
    }
    
    res[[signature.key]] <- res[[signature.key.new.up]] - res[[signature.key.new.down]]
  }
  
  res <- res[, names(signature), drop = F]

  return(res)
}

#' Run the singscore signature score
#' 
#' Function to run the singscore signature score
#' 
#' @param X matrix; The matrix. Columns = cells, rows = some features
#' @param ranks rank-matrix as defined by singscore::rankGenes(). 
#' If provided, the function does not compute the rankings from X
#' @param signature list: signature as defined in GetSignature()
#' @param scale logical: Weather to scale signature score or not.
#' 
#' @return data.frame; The signature score data.frame. Row = cells and columns = names(signature)
#' 
#' @export
RunSigScoreSingscore <- function(X, signature, ranks = NULL, scale = F){
  
  if (is.null(ranks)){
    ranks = singscore::rankGenes(as.matrix(X))
  }
  
  missing.signature <- c()
  
  signature.one.sided.up <- list()
  signature.one.sided.down <- list()
  signature.input.up <- list()
  signature.input.down <- list()
  for (signature.key in names(signature)){
    if (length(names(signature[[signature.key]])) == 1){
      if ("up" %in% names(signature[[signature.key]])){
        if (length(signature[[signature.key]]$up) > 0){
          signature.one.sided.up[[signature.key]] <- 
            GSEABase::GeneSet(signature[[signature.key]]$up, setName = signature.key)
        } else {
          missing.signature <- c(missing.signature, signature.key)
        }
      } else {
        if (length(signature[[signature.key]]$down) > 0){
          signature.one.sided.down[[signature.key]] <- 
            GSEABase::GeneSet(signature[[signature.key]]$down, setName = signature.key)
        } else {
          missing.signature <- c(missing.signature, signature.key)
        }
      }
    } else {
      if ((length(signature[[signature.key]]$up) * length(signature[[signature.key]]$down)) > 0){
        signature.input.up[[signature.key]] <- 
          GSEABase::GeneSet(signature[[signature.key]]$up, setName = signature.key)
        signature.input.down[[signature.key]] <- 
          GSEABase::GeneSet(signature[[signature.key]]$down, setName = signature.key)
      } else {
        if (length(signature[[signature.key]]$up) > 0){
          signature.one.sided.up[[signature.key]] <- 
            GSEABase::GeneSet(signature[[signature.key]]$up, setName = signature.key)
        } else if (length(signature[[signature.key]]$down) > 0){
          signature.one.sided.down[[signature.key]] <- 
            GSEABase::GeneSet(signature[[signature.key]]$down, setName = signature.key)
        } else {
          missing.signature <- c(missing.signature, signature.key)
        }
      }
    }
  }
  signature.one.sided.up <- GSEABase::GeneSetCollection(signature.one.sided.up)
  signature.one.sided.down <- GSEABase::GeneSetCollection(signature.one.sided.down)
  signature.input.up <- GSEABase::GeneSetCollection(signature.input.up)
  signature.input.down <- GSEABase::GeneSetCollection(signature.input.down)
  
  if (F){
    test.up <- list()
    test.dn <- list()
    
    test.up[["test.1"]] <- GSEABase::GeneSet(c("CXCL13", "LAG3", "RBPJ"), setName = "test.1")
    test.dn[["test.1"]] <- GSEABase::GeneSet(c("TOX", "TNFRSF9", "CTLA4"), setName = "test.1")
    test.up[["test.2"]] <- GSEABase::GeneSet(c("CXCL13", "LAG3", "RBPJ"), setName = "test.2")
    test.dn[["test.2"]] <- GSEABase::GeneSet(c("LAG3"), setName = "test.2")
    test.up <- GSEABase::GeneSetCollection(test.up)
    test.dn <- GSEABase::GeneSetCollection(test.dn)
    res.test <- singscore::multiScore(rankData = ranks, 
                                      upSetColc = test.up, 
                                      downSetColc = test.dn, 
                                      centerScore = T, 
                                      knownDirection = T)
    res.test <- data.frame(t(res.test$Scores), check.names = F, check.rows = F)
  }
  
  
  res <- data.frame(row.names = colnames(X))
  if (length(signature.one.sided.up) > 0){
    res.tmp <- singscore::multiScore(rankData = ranks, 
                                 upSetColc = signature.one.sided.up, 
                                 centerScore = T, 
                                 knownDirection = T)
    if (!is.null(res.tmp$Scores)){
      res.tmp <- data.frame(t(res.tmp$Scores), check.names = F, check.rows = F)
    } else {
      res.tmp <- data.frame(row.names = colnames(X), check.rows = F)
    }
    missing_columns <- names(signature.one.sided.up)[!(names(signature.one.sided.up) %in% colnames(res.tmp))]
    for (col_ in missing_columns){
      res.tmp[[col_]] <- 0
    }
    res <- cbind(res, res.tmp)
  }
  if (length(signature.one.sided.down) > 0){
    res.tmp <- singscore::multiScore(rankData = ranks, 
                                 upSetColc = signature.one.sided.down, 
                                 centerScore = T, 
                                 knownDirection = T)
    
    if (!is.null(res.tmp$Scores)){
      res.tmp <- data.frame(t(res.tmp$Scores), check.names = F, check.rows = F)
      res.tmp <- -res.tmp
    } else {
      res.tmp <- data.frame(row.names = colnames(X), check.rows = F)
    }
    missing_columns <- names(signature.one.sided.down)[!(names(signature.one.sided.down) %in% colnames(res.tmp))]
    for (col_ in missing_columns){
      res.tmp[[col_]] <- 0
    }
    res <- cbind(res, res.tmp)
  }
  if (length(signature.input.up) > 0){
    res.tmp <- singscore::multiScore(rankData = ranks, 
                                 upSetColc = signature.input.up, 
                                 downSetColc = signature.input.down, 
                                 centerScore = T, 
                                 knownDirection = T)
    if (!is.null(res.tmp$Scores)){
      res.tmp <- data.frame(t(res.tmp$Scores), check.names = F, check.rows = F)
    } else {
      res.tmp <- data.frame(row.names = colnames(X), check.rows = F)
    }
    missing_columns <- names(signature.input.up)[!(names(signature.input.up) %in% colnames(res.tmp))]
    for (col_ in missing_columns){
      res.tmp[[col_]] <- 0
    }
    res <- cbind(res, res.tmp)
  }
  
  if (length(missing.signature) > 0){
    warning(paste0("In RunSigScoreSingscore: Missing genes in signature(s): ",
                   paste(missing.signature, collapse = ", ")))
    for (col_ in missing.signature){
      res[[col_]] <- 0
    }
  }
  return(res) 
}

