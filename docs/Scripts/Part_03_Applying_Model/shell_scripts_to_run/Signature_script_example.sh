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
path_config=~/Documents/Demo/MixTRTpred/Scripts/Part_03_Applying_Model/Model_configs/Signature_models.yml
config_active="RNA_pca"
# possible configs names are the 9 Signature model names: i.e. 
# RNA_data_Seurat_wilcox, RNA_data_Seurat_LR, RNA_data_Seurat_negbinom, RNA_data_DESeq2_Wald, RNA_data_DESeq2_LRT, RNA_data_edgeR_QFL, RNA_data_edgeR_LRT, RNA_data_limma_trend, RNA_data_limma_voom

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

 