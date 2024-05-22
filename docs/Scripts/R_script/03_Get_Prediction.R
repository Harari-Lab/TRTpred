#' Script to the the predition of some tumour reactive model

# ------------------------------------------------------------------------------
# Prepare environement
# ------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyverse)
library(readxl)
library(getopt)
library(TRTpred)

options(expressions = 5e5) # max = 5e5

# ------------------------------------------------------------------------------
# Prepare variable
# ------------------------------------------------------------------------------
option_specification = matrix(c(
  'path.config', 'f', 1, 'character', 'The path of the config file',
  'config.active', 'g', 1, 'character', 'the config (default = "default"',
  'prediction.folder', 'j', 1, 'character', 'The path where to save the data',
  'path.mca', 'a', 1, 'character', 'The path of the seurat RDS object (default: "~/Documents/ATATIL_Remy/data/scRNA/ATATIL_CD45_normalized_v5.rds" or )',
  'mca.assay', 'b', 2, 'character', 'mca assay default = RNA',
  'mca.slot', 'c', 2, 'character', 'mca slot default = RNA',
  'help', 'h', 0, "logical", 'Help'
), byrow=TRUE, ncol=5);
in.options = getopt(option_specification)


if ( !is.null(in.options$help) ) {
  cat(getopt::getopt(option_specification, usage=TRUE))
  q(status=1)
}

prediction.folder.save <- in.options$prediction.folder
if (is.null(prediction.folder.save)){
  prediction.folder.save <- "~/Documents/Demo/MixTRTpred/scripts/Part_03_Applying_Model/results/"
}

path.config <- in.options$path.config
if (is.null(path.config)){
  path.config <- "~/Documents/Demo/MixTRTpred/scripts/Part_02_Final_Model_Training/model_configs/Signature_LOPO_CV.yml"
}

config.active <- in.options$config.active
if (is.null(config.active)){
  config.active <- "RNA_data_edgeR_QFL"
}

path.mca <- in.options$path.mca
if (is.null(path.mca)){
  path.mca <- "<Path_To_Some_External_Data>"

}
mca.assay <- in.options$mca.assay
if (is.null(mca.assay)){
  mca.assay <- "RNA"
}
mca.slot <- in.options$mca.slot
if (is.null(mca.slot)){
  mca.slot <- "scale.data"
}

cat("Running Get_Prediction.R with parameters: \n")
cat(paste0("\t- prediction.folder.save: ", prediction.folder.save, "\n"))
cat(paste0("\t- path.mca: ", path.mca, "\n"))
cat(paste0("\t- mca.assay: ", mca.assay, "\n"))
cat(paste0("\t- mca.slot: ", mca.slot, "\n"))
cat(paste0("\t- path.config: ", path.config, "\n"))
cat(paste0("\t- config.active: ", config.active, "\n")) 

# ------------------------------------------------------------------------------
# Prepare variable
# ------------------------------------------------------------------------------

if (!file.exists(path.config)){
  stop("path.config (-f) does not exist")
}

config.file <- ProcessConfigFile(file = path.config, config = config.active)
str.out <- config.file$str.out
config.file <- config.file$config.file

project.name <- config.file$name # Yes
params.data <- config.file$data # Label of TR
params.output <- config.file$output
params.feature.trans <- config.file$feature.trans
params.DA <- config.file$DA
params.model <- config.file$model
params.CV <- config.file$CV
params.run <- config.file$run

# sanity checks in configs:
if (!file.exists(params.data$path)){
  stop("path.mca does not exist")
}
if (!dir.exists(params.output$folder.save)){
  stop("folder.save does not exist")
} 
if (!is.null(params.output$folds.save.folder)){
  if (!dir.exists(params.output$folds.save.folder)){
    warning("params.output$folds.save.folder does not exist --- Continue running without this")
  }
}

cat("The config files has the following parameters: \n")
cat(str.out)

# ------------------------------------------------------------------------------
# Load data from model
# ------------------------------------------------------------------------------
mca.train <- readRDS(file = params.data$path)

mca.original.project.name <- mca.train@project.name

rm(mca.train)

folder.save <- paste0(gsub("/$", "", params.output$folder.save), "/", mca.original.project.name, "/")
if (!dir.exists(folder.save)){
  stop("folder.save does not exist: ", folder.save)
}

folds.save.folder <- paste0(gsub("/$", "", params.output$folds.save.folder), "/", 
                            mca.original.project.name, "/", 
                            project.name, "/")
if (!dir.exists(folds.save.folder)){
  stop()
}

file.save.name.res <- paste0(project.name, "_final_model_res", ".rds")

model.res <- readRDS(file = paste0(folder.save, file.save.name.res))

hyperparameter.df <- model.res
hyperparameter.df <- hyperparameter.df[, !grepl(pattern = "((train)|(test)|(n.data)|(p.data))\\.", x = colnames(hyperparameter.df))]

if (all(c("scale.mean", "scale.sd") %in% colnames(model.res))){
  scale.list <- as.list(model.res[, c("scale.mean", "scale.sd")])
  hyperparameter.df$scale.mean <- NULL; hyperparameter.df$scale.sd <- NULL
} else {
  scale.list <- NULL
}


for (hyperparam_ in colnames(hyperparameter.df)){
  params.model$hyperparameters[[hyperparam_]] <- hyperparameter.df[[hyperparam_]]
}
# ------------------------------------------------------------------------------
# Load data to get prediction form 
# ------------------------------------------------------------------------------

mca <- readRDS(file = path.mca)

if (!(mca.assay %in% names(mca@assays))){
  stop("mca.assay not in mca's assays")
}

# Get data.input
data.input <- Seurat::GetAssayData(object = mca,
                                   assay = mca.assay, 
                                   slot = mca.slot)
data.input <- data.frame(t(data.input), check.rows = F, check.names = F)

data.output <- mca@meta.data

# ------------------------------------------------------------------------------
# Define where to save the results
# ------------------------------------------------------------------------------

prediction.folder.save <- paste0(gsub("/$", "", prediction.folder.save), "/", mca@project.name, "/")
if (!dir.exists(prediction.folder.save)){
  system(command = paste("mkdir", "-p", prediction.folder.save))
}
prediction.file.path <- paste0(prediction.folder.save, project.name, "_prediction.rds")

# ------------------------------------------------------------------------------
# Get predictions
# ------------------------------------------------------------------------------
time_start <- proc.time()

prediction.res <- GetModelPrediction(data.input = data.input, 
                                     data.output = data.output, 
                                     y.label = params.data$label, 
                                     method = params.model$method, 
                                     hyperparameters = params.model$hyperparameters,
                                     path.folder = folds.save.folder,
                                     data.trans.method = params.feature.trans$method,
                                     DEA.file.name = "dea.res.rds",
                                     signature.x.threshold = params.model$hyperparameters$signature.x.threshold,
                                     signature.score.scale = T, 
                                     signature.score.scale.params = scale.list,
                                     LR.model.file.name = NULL, 
                                     LR.prob.threhsold = 0.5)

cat(paste0("Save prediction.res in :", prediction.file.path, " \n"))
saveRDS(prediction.res, file = prediction.file.path)

time_end <- proc.time()

time_elapsed <- time_end - time_start

cat("--- Information on computational time --- \n")
cat(paste(names(time_elapsed), collapse = "\t"))
cat("\n")
cat(paste(round(time_elapsed, digits = 5), collapse = "\t"))
cat("\n")

cat("FINISHED\n")
