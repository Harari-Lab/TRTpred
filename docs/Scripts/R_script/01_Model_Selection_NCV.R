# Script to prepare different type of data

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prepare environment
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
  'help', 'h', 0, "logical", 'Help'
), byrow=TRUE, ncol=5);
in.options = getopt(option_specification)

if ( !is.null(in.options$help) ) {
  cat(getopt::getopt(option_specification, usage=TRUE))
  q(status=1)
}

path.config <- in.options$path.config
if (is.null(path.config)){
  path.config <- "~/Documents/Demo/MixTRTpred/scripts/Part_01_Model_Selection/model_configs/LR_LOPO_NCV.yml"
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
  stop("path.mca does not exist: ", params.data$path)
}
if (!dir.exists(params.output$folder.save)){
  stop("folder.save does not exist: ", params.output$folder.save)
} 
if (!is.null(params.output$folds.save.folder)){
  if (!dir.exists(params.output$folds.save.folder)){
    warning("params.output$folds.save.folder does not exist --- Continue running without this")
  }
}

if (is.null(params.run$loadComputedData)){
  params.run$loadComputedData <- FALSE
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
# Cross validation
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Define File save names
# --- Inner NCV file: 
file.save.name <- paste0(project.name, "_NCV_inner_loop", ".rds")
file.save.path.inner <- paste0(folder.save, file.save.name)
# --- Outer NCV file: 
file.save.name <- paste0(project.name, "_NCV_outer_loop", ".rds")
file.save.path.outer <- paste0(folder.save, file.save.name)

if (params.run$loadComputedData){
  if (file.exists(file.save.path.inner)){
    file.save.path.inner.arg <- file.save.path.inner
  } else {
    file.save.path.inner.arg <- NULL
  }
  if (file.exists(file.save.path.outer)){
    file.save.path.outer.arg <- file.save.path.outer
  } else {
    file.save.path.outer.arg <- NULL
  }
} else {
  file.save.path.inner.arg <- NULL
  file.save.path.outer.arg <- NULL
}

NCV.res <- NestedCrossValidation(
  x = Data.list$data.input, 
  y = Data.list$data.output,
  y.label = params.data$label,
  y.sample = params.data$sample,
  y.covariates = params.data$covariates,
  CV.k.in = params.CV$k.inner, 
  CV.k.out = params.CV$k.outer, 
  CV.y.leave.1.out.column = params.CV$leave.1.out,
  CV.replicates = params.CV$replicates,
  CV.sampling.method = params.CV$sampling.method,
  method = params.model$method,
  metric.to.optimize = params.CV$selection.metric,
  hyperparams = params.model$hyperparameters,
  data.trans.method = params.feature.trans$method,
  pca.explained.var.threshold = params.feature.trans$pca.explained.var.threshold,
  rm.corr.features = params.feature.trans$rm.mutually.corr.features,
  DA.method = params.DA$method, 
  DA.p.adjust.method = params.DA$p.adjust.method,
  DA.p.val.threshold = params.DA$p.val.threshold, 
  DA.event.per.variable = params.DA$event.per.variable,
  LR.prob.threhsold = params.model$prob.threshold,
  signature.col.aggregate = params.data$DEA.pseudobulk.col,
  DEA.method = params.model$DEA.method,
  DEA.data = Data.list$DEA.data,
  sample.weights = Data.list$sample.weights,
  folds.save.folder = params.output$folds.save.folder, 
  file.NCV.inner = file.save.path.inner.arg, 
  file.NCV.outer = file.save.path.outer.arg,
  seed = params.run$seed)

# Save results:
# --- Inner NCV file: 
cat(paste0("Save outer inner NCV results in:", file.save.path.inner, " \n"))
saveRDS(object = NCV.res$NCV.inner.res, file = file.save.path.inner)

# --- Outer NCV file: 
cat(paste0("Save outer CV results in:", file.save.path.outer, " \n"))
saveRDS(object = NCV.res$NCV.outer.res, file = file.save.path.outer)


time_end <- proc.time()

time_elapsed <- time_end - time_start

cat("--- Information on computational time --- \n")
cat(paste(names(time_elapsed), collapse = "\t"))
cat("\n")
cat(paste(round(time_elapsed, digits = 5), collapse = "\t"))
cat("\n")

cat("FINISHED\n")
