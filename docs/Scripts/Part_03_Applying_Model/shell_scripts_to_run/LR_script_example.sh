#!/bin/bash

#SBATCH --time=72:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mem=256G   # maximum memory per node
#SBATCH --job-name="<job-name>"
#SBATCH --account="<account-name>"
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source /dcsrsoft/spack/bin/setup_dcsrsoft
module load singularity

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define the directory of the home and sif file:
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home_dir=~/Documents/Demo/MixTRTpred/Scripts/Part_03_Applying_Model/
sif_dir=~/Documents/Demo/MixTRTpred/Scripts/Computation_Environement/MixTRTpred_singularity.sif

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define the scripts:
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
script_dir=~/Documents/Demo/MixTRTpred/Scripts/R_script/03_Get_Prediction.R

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define the parameters:
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_config=~/Documents/Demo/MixTRTpred/Scripts/Part_03_Applying_Model/Model_configs/LR_models.yml
config_active="RNA_pca"
# possible configs names are the 12 LR model names: i.e. 
# RNA_raw, RNA_raw_DA_wilcox_fdr, RNA_raw_DA_Seurat_wilcox, RNA_raw_DA_Seurat_LR, RNA_raw_DA_Seurat_negbinom, RNA_raw_DA_limma_trend, RNA_raw_DA_limma_voom, RNA_raw_DA_edgeR_QFL, RNA_raw_DA_edgeR_LRT, RNA_raw_DA_DESeq2_Wald, RNA_raw_DA_DESeq2_LRT, RNA_pca, RNA_pca_DA_wilcox_fdr

prediction_folder_save="~/Documents/Demo/MixTRTpred/~/Documents/Demo/MixTRTpred/Scripts/Part_03_Applying_Model/results/"
path_mca="<path_to_data>"
mca_assay=RNA
mca_slot=scale.data


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run scripts
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# run script with singularity:
singularity exec -B ${home_dir} ${sif_dir} \
Rscript ${script_dir} \
    --path.config ${path_config} \
    --config.active ${config_active} \
    --prediction.folder ${prediction_folder_save} \
    --path.mca ${path_mca} \
    --mca.assay ${mca_assay} \
    --mca.slot ${mca_slot}
    
exit 0

 
