[![DOI](https://img.shields.io/badge/DOI%3A-10.1038/s41587--024--02232--0-brightgreen)](https://doi.org/10.1038/s41587-024-02232-0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11401545.svg)](https://doi.org/10.5281/zenodo.11401545)
[![License](https://img.shields.io/badge/License-LICENSE-green)](./LICENSE_TRTpred.pdf)

# TRTpred library

## Introduction

R library providing the functions to train, evaluate, and apply tumor
reactive T cell predictors referred as to TRTpred. All information about
the predictors and evaluation framework can be found at Pétremand R.,
*et al*. Nature Biotechnology (2024)
[![DOI](https://img.shields.io/badge/DOI%3A-10.1038/s41587--024--02232--0-brightgreen)](https://doi.org/10.1038/s41587-024-02232-0)

In short, TRTpred library is composed of functions to perform different
tasks:

1.  Perform models selection using NCV
2.  Train models using CV
3.  Apply a model on a new data
4.  Evaluate models

Importantly, a demo is added to illustrate TRTpred applicability to new
data.

## Quick Installation guide

For quick test of TRTpred

1.  Download TRTpred folder (e.g. in /user/username/Document/TRTpred)
2.  Download docker container with the working R environement for
    TRTpred in
    [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11401545.svg)](https://doi.org/10.5281/zenodo.11401545)
3.  Load and run the docker container in the terminal using:

``` r
docker load -i <path to TRTpred_docker.tar>
docker run -e PASSWORD=1234 -p 8787:8787 -v /user/username/Document/TRTpred:/home trtpred:devel1.0
```

1.  Open <http://localhost:8787/> in a browser and log in rstudio server
    with login:rstudio and password:1234
2.  From within the rstudio terminal, locate the TRTpred folder and
    install TRTpred as described above in the “Installation” section.

``` r
R CMD build TRTpred
R CMD INSTALL TRTpred_0.0.1.tar.gz
```

1.  Optional: Test the package in docs/TRTpred_Demo_01.rmd

## Installation

You may install the package locally or through github

### Local Installation

Otherwise, you may pull `TRTpred` and open a terminal in the directory
of the `TRTpred` folder to compile the library by using this command
line `R CMD build TRTpred`. A file called “TRTpred_0.0.1.tar.gz” should
have been created. Finally, you may install the package using this the R
command line as follows: `R CMD INSTALL TRTpred_0.0.1.tar.gz`

### Github Installation

You may install this library directly from this Git repository by
running the following code within R or RStudio:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github(repo = "https://github.com/Harari-Lab/TRTpred")
```

## System Requirements

### Hardware requirements

A computer with 8-16GB of ram and a 2.5GHz CPU can run TRTpred in a
couple of minutes. The computation and time bottlenecks are attributed
to the quantity of the data. Unfortunately, the model selection with the
nested cross-validation is highly time consuming, particularly in the
Leave-One-Patient-Out (LOPO) mode. To minimize this computation time,
you may run this part using parallel programming (see *doParallel*
library).

### Software requirements

TRTpred is built on R. Users will need a version of R higher than 4.0.3
and a compatible R-studio. The required R libraries and version are
described below.

| **Library**   | **Version** | **Source**   | **Installation line**                                                  |
|----------|---------|---------|---------------------------------------------|
| *AUCell*      | 1.16.0      | BioConductor | BiocManager::install(“AUCell”)                                         |
| *config*      | 0.3.1       | CRAN         | install.packages(“config”)                                             |
| *DESeq2*      | 1.34.0      | BioConductor | BiocManager::install(“DESeq2”)                                         |
| *doParallel*  | 1.1.17      | CRAN         | install.packages(“doParallel”, repos=“<http://R-Forge.R-project.org>”) |
| *dplyr*       | 1.1.0       | CRAN         | install.packages(“dplyr”)                                              |
| *edgeR*       | 3.36.0      | BioConductor | BiocManager::install(“edgeR”)                                          |
| *foreach*     | 1.5.2       | CRAN         | install.packages(“foreach”)                                            |
| *getopt*      | 1.20.3      | CRAN         | install.packages(“geopt”)                                              |
| *glmnet*      | 4.1-6       | CRAN         | install.packages(“glmnet”)                                             |
| *glmnetUtils* | 1.1.8       | CRAN         | install.packages(“glmnet”)                                             |
| *GSEABase*    | 1.56.0      | github       | devtools::install_github(“hongooi73/glmnetUtils”)                      |
| *limma*       | 3.50.3      | BioConductor | BiocManager::install(“limma”)                                          |
| *matrixStats* | 0.63.0      | CRAN         | install.packages(“matrixStats”)                                        |
| *readxl*      | 1.4.2       | CRAN         | install.packages(“readxl”)                                             |
| *ropls*       | 1.26.4      | BioConductor | BiocManager::install(“ropls”)                                          |
| *Seurat*      | 4.3.0       | github       | remotes::install_version(“Seurat”, version = “4.0.3”)                  |
| *singscore*   | 1.14.0      | BioConductor | BiocManager::install(“singscore”)                                      |
| *stats*       | 4.1.2       | CRAN         | install.packages(“stats”)                                              |
| *tidyr*       | 1.3.0       | github       | devtools::install_github(“hadley/tidyr”)                               |
| *tidyverse*   | 2.0.0       | CRAN         | install.packages(“tidyverse”)                                          |
| *UCell*       | 1.3.1       | BioConductor | BiocManager::install(“UCell”)                                          |

To install these libraries you will need *remotes*
(“install.packages(‘remotes’)”), *devtools*
(“install.packages(‘devtools’)”) and *BiocManager*
(“install.packages(”BiocManager”)”).

If you prefer, you may run TRTpred in a docker container. See next
section for more information.

### Docker Container

A **docker container** containing the computation environment can be
found in zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11401545.svg)](https://doi.org/10.5281/zenodo.11401545)

To install the docker image, you may first download the
TRTpred_docker.tar file and load it onto docker using this line in the
terminal.

``` r
docker load -i <path to TRTpred_docker.tar>
```

Where, “\<path to TRTpred_docker.tar\>” is the path of the
TRTpred_docker.tar file that you just downloaded. This will result in
loading a docker container called trtpred:devel1.0. you may check by
writing `docker images`.

To run this docker and have access to the TRTpred folders, we advise you
to mount a volume linking your local files to the “home” folder in the
container. To do so, you may run this line:

``` r
docker run -e PASSWORD=1234 -p 8787:8787 -v <Path to your documents>:/home trtpred:devel1.0
```

Where “<Path to your documents>” is the local path to your documents for
example (e.g. /Users/username/Document/TRTpred).

Running this docker will give you access to a rstudio server with the a
R environment necessary to run TRTpred. To access this Rstudio server,
you may open <http://localhost:8787/> in a browser and log in rstudio
server with login:rstudio and password:1234.

From within the rstudio terminal, locate the TRTpred folder and install
TRTpred as described above in the “Installation” section.

``` r
R CMD build TRTpred
R CMD INSTALL TRTpred_0.0.1.tar.gz
```

## Documentation

See a description of the functions implemented [here](docs/functions.md)

## Tutorials

A demonstration of TRTpred can be found
[here](https://html-preview.github.io/?url=https://github.com/Harari-Lab/TRTpred/blob/main/docs/TRTpred_Demo_01.html)
in docs/TRTpred_Demo_01.html

## Scripts

The “script” folder contains the scripts (shell and R scripts) that have
been used to run TRTpred.

The “script” folder contains a sub-folder called “R_script” that is
composed of three R-script:

1.  01_Model_Selection_NCV.R: A script dedicated to the model selection
    through Leave-One-Patient-Out (LOPO) Nested Cross Validation (NCV).
    It works by using configuration files of the models (see below).
2.  02_Model_Training.R: A script dedicated to the final model training
    (and optimization though Cross Validation (CV)). It works by using
    configuration files of the models (see below).
3.  03_Get_Prediction.R: A script dedicated to the final model
    application. It works by using configuration files of the models
    (see below).

The “script” folder also contains 3 sub-folders, each dedicated for one
task:

1.  Part_01_Model_Selection: Script that runs 01_Model_Selection_NCV.R
    R-script with the correct configs and data
2.  Part_02_Final_Model_Training: Script that runs 02_Model_Training.R
    R-script with the correct configs and data
3.  Part_03_Applying_Model: Script that runs 03_Get_Prediction.R
    R-script with the correct configs and data

Every sub-folder has 3 sub-folders:

1.  Model_configs: A folder with the model’s configuration files (YML
    format)
2.  Shell_script: Examples of the scripts that have been used compute
    these parts
3.  Results: Results folder

## Software

The product is provided free of charge, and, therefore, on an “as is”
basis, without warranty of any kind.

For scientific questions, please contact [Rémy
Pétremand](mailto:remy.petremand@chuv.ch). FOR-PROFIT USERS: If you plan
to use or any data provided with the script in any for-profit
application, you are required to obtain a separate license. To do so,
please contact [Nadette Bulgin](mailto:nbulgin@lcr.org) at the Ludwig
Institute for Cancer Research Ltd.

## Notes for users

TRTpred was build on CD8+ tumor-infiltrating lymphocytes (TILs). The
input of the package are the CD8+ TILs single-cell RNA-sequencing data.

For scientific questions, please contact [Rémy
Pétremand](mailto:remy.petremand@chuv.ch).

For license-related questions, please contact [Nadette
Bulgin](mailto:nbulgin@lcr.org) at the Ludwig Institute for Cancer
Research Ltd.

## Citations

-   Pétremand, R. et al. Identification of clinically relevant T cell
    receptors for personalized T cell therapy using combinatorial
    algorithms. Nat Biotechnol (2024).
    <https://doi.org/10.1038/s41587-024-02232-0>
