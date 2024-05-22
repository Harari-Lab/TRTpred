# Script to run the final model

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prepare environmenta
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_start <- proc.time()

suppressMessages(library(Seurat))
suppressMessages(require(parallel))
suppressMessages(require(doParallel))
suppressMessages(library(getopt))
suppressMessages(library(TRTpred))

options(expressions = 5e5) # max = 5e5

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prepare variable
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option_specification = matrix(c(
  'path.config', 'f', 1, 'character', 'The path of the config file',
  'config.active', 'g', 2, 'character', 'the config (default = "default"',
  'read.CVres', 'a', 2, 'logical', 'Read hyperparmeters from the CV results',
  'help', 'h', 0, "logical", 'Help'
), byrow=TRUE, ncol=5);
in.options = getopt(option_specification)

if ( !is.null(in.options$help) ) {
  cat(getopt::getopt(option_specification, usage=TRUE))
  q(status=1)
}

read.CVres <- in.options$read.CVres
if (is.null(read.CVres)){
  read.CVres <- F
} else {
  read.CVres <- T
}

path.config <- in.options$path.config
if (is.null(path.config)){
  path.config <- "~/Documents/Demo/MixTRTpred/scripts/Part_02_Final_Model_Training/model_configs/LR_LOPO_CV.yml"
}

config.active <- in.options$config.active
if (is.null(config.active)){
  config.active <- "RNA_raw"
}


if (!file.exists(path.config)){
  stop("path.config (-f) does not exist")
}

config.file <- ProcessConfigFile(file = path.config, config = config.active)
str.out <- config.file$str.out
config.file <- config.file$config.file

project.name <- config.file$name
params.data <- config.file$data
params.output <- config.file$output
params.feature.trans <- config.file$feature.trans
params.DA <- config.file$DA
params.model <- config.file$model
params.CV <- config.file$CV
params.run <- config.file$run

# register the number of Cores for the computation: 
if (params.run$nCores >  parallel::detectCores()){
  message(paste0("nCores > ", parallel::detectCores(), " (= parallel::detectCores()): So nCores is set to parallel::detectCores()"))
  params.run$nCores <-  parallel::detectCores()
}
if (Sys.info()["nodename"] %in% c("HOS56290", "hos71643")){
  params.run$nCores <- 1
}
doParallel::registerDoParallel(params.run$nCores)

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

cat("Running CV_prepare_data.R with parameters: \n")
cat(paste0("\t- path.config: ", path.config, "\n"))
cat(paste0("\t- config.active: ", config.active, "\n")) 
cat("The config files has the following parameters: \n")
cat(str.out)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("Load data: \n")

# Get the Seurat single-cell data:
mca <- readRDS(file = params.data$path)

# Get data in correct format
Data.list <- PrepareTrainingDataFromSeurat(
  SeuratObject = mca, 
  label = params.data$label, 
  label.order = params.data$label.order, 
  assay = params.data$assay, 
  slot = params.data$slot, 
  DEA.assay = params.data$DEA.assay, 
  DEA.slot = params.data$DEA.slot, 
  DEA.pseudobulk.col = params.data$DEA.pseudobulk.col, 
  covariates = params.data$covariates, 
  sample.weights.col = params.model$sample.weights)

# Define the output folder
folder.save <- paste0(gsub("/$", "", params.output$folder.save), "/", mca@project.name, "/")
if (!dir.exists(folder.save)){
  system(command = paste("mkdir", "-p", folder.save))
}
if (!is.null(params.output$folds.save.folder)){
  params.output$folds.save.folder <- paste0(gsub("/$", "", params.output$folds.save.folder), "/", 
                                            mca@project.name, "/", 
                                            project.name, "/")
  if (!dir.exists(params.output$folds.save.folder)){
    system(command = paste("mkdir", "-p", params.output$folds.save.folder))
  }
}

# Remove mca (not needed)
rm(mca)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get Hyperparams from CV res
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (read.CVres){
  CV.file.name <- paste0(project.name, "_CV", ".rds")# _final_CV
  CV.file.path <- paste0(gsub("/final_model/", "/CV/", folder.save), CV.file.name)
  if (file.exists(CV.file.path)){
    CV.res <- readRDS(file = CV.file.path)
  } else {
    stop("Parameter read.CVres but no CV results in ", CV.file.path)
  }
  
  if (is.null(params.CV$selection.metric)){
    params.CV$selection.metric <- "mcc"
  }
  
  hyperparams.best <- GetBestHyperparams(
    CV.res = CV.res,
    CV.type = "CV",
    method = params.model$method,
    selection.metric = params.CV$selection.metric)
  
  hyperparams.col <- colnames(hyperparams.best)[!colnames(hyperparams.best) %in% c("acc_avg", "acc_max", "acc_min")]
  for (hyperparam_ in hyperparams.col){
    params.model$hyperparameters[[hyperparam_]] <- hyperparams.best[[hyperparam_]]
  }
}  

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Train Model
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("Create final model: \n")

file.save.name <- paste0(project.name, "_final_model_res", ".rds")
file.save.path <- paste0(folder.save, file.save.name)

res.model <- TrainModel(
  x.train = Data.list$data.input, 
  y.train = Data.list$data.output,
  x.test = NULL, 
  y.test = NULL,
  y.label = params.data$label,
  y.sample = params.data$sample,
  y.covariates = params.data$covariates,
  method = params.model$method,
  hyperparams = params.model$hyperparameters,
  data.trans.method = params.feature.trans$method,
  pca.explained.var.threshold = params.feature.trans$pca.explained.var.threshold,
  rm.corr.features = params.feature.trans$rm.mutually.corr.features,
  DA.method = params.DA$method, 
  DA.p.adjust.method = params.DA$p.adjust.method,
  DA.p.val.threshold = params.DA$p.val.threshold, 
  DA.event.per.variable = params.DA$event.per.variable,
  LR.prob.threhsold = params.model$prob.threshold,
  DEA.method = params.model$DEA.method,
  DEA.data = Data.list$DEA.data,
  sample.weights = Data.list$sample.weights,
  signature.col.aggregate = params.data$DEA.pseudobulk.col, 
  folds.save.folder = params.output$folds.save.folder,
  save.fold.model = T, 
  save.fold.pred = T, 
  save.fold.score = T)


cat(paste0("Save outer CV results in:", file.save.path, " \n"))
saveRDS(object = res.model, file = file.save.path)

time_end <- proc.time()

time_elapsed <- time_end - time_start

cat("--- Information on computational time --- \n")
cat(paste(names(time_elapsed), collapse = "\t"))
cat("\n")
cat(paste(round(time_elapsed, digits = 5), collapse = "\t"))
cat("\n")

cat("FINISHED\n")
