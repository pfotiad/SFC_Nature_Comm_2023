#!/bin/bash

# The purpose of this script is to register the HCP-output T1wDividedByT2w.nii.gz and T1wDividedByT2w_ribbon.nii.gz files provided with the Human Connectome Project (HCP)
# download from native to MNI152 space.

## PREREQUISITES:

# Before running the script you need to download:
# 1. The publicly available 100 unrelated HCP dataset distribution: https://db.humanconnectome.org/
# 2. The latest HCP pipelines: https://github.com/Washington-University/HCPpipelines/releases

## RUNNING THE SCRIPT:

# Adapt all paths below to the ones corresponding to your data
# In order to run the script, type into the terminal: T1w_T2w_ratio_native_to_MNI.sh

## OUTPUT:

# 1. A T1wDividedByT2w_MNI152.nii.gz volumetric file for each subject corresponding to the T1w/T2w signal intensity mapped into MNI152 standardized space.
# 2. A T1wDividedByT2w_ribbon_MNI152.nii.gz volumetric file for each subject corresponding to the T1w/T2w signal intensity only within the cortex mapped into MNI152 standardized space.

# Author: Panagiotis Fotiadis, 2021-2023

# If parts of the following code are used, please cite the following paper:

# Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD,
# Shinohara RT, and Bassett DS “Myelination and excitation-inhibition
# balance synergistically shape structure-function coupling across the
# human cortex,” (2023) Nature Communications.


## -------------------------------------------------------------------------------------------------------------------------------------------

# Define the paths of interest
subj_dir=/cbica/home/panosf/Q7/structure_function_coupling/HCP_Processing_Pipeline # path to where the HCP subjects are located
subj_file=/cbica/home/panosf/Q7/structure_function_coupling/subject_list.txt # contains the list of the 100 HCP subjects (i.e., a row per subject ID)

# Count the number of subjects
count=$(cat $subj_file | wc -l)

# Loop through the subjects
for (( i=1; i<=$count; i++ ));
do
  subj=$(cat $subj_file | head -$i | tail -1 | awk '{print $1}')

  # Input files: T1w/T2w files in native space
  t1w_t2w_ratio=$subj_dir/$subj/T1w/T1wDividedByT2w.nii.gz
  t1w_t2w_ratio_ribbon=$subj_dir/$subj/T1w/T1wDividedByT2w_ribbon.nii.gz

  # Transformation file
  transformation_file=$subj_dir/$subj/MNINonLinear/xfms/acpc_dc2standard.nii.gz

  # The MNI template
  mni_template=/cbica/home/panosf/software/HCPpipelines-4.3.0/global/templates/MNI152_T1_1mm_brain.nii.gz

  t1w_acpc_dc=$subj_dir/$subj/T1w/T1w_acpc_dc.nii.gz

  # Apply the transformation to map the input files onto MNI space
  wb_command -volume-warpfield-resample $t1w_acpc_dc $transformation_file $mni_template ENCLOSING_VOXEL $subj_dir/$subj/T1w/T1w_acpc_dc_MNI152.nii.gz -fnirt $t1w_acpc_dc
  wb_command -volume-warpfield-resample $t1w_t2w_ratio $transformation_file $mni_template ENCLOSING_VOXEL $subj_dir/$subj/T1w/T1wDividedByT2w_MNI152.nii.gz -fnirt $t1w_acpc_dc
  wb_command -volume-warpfield-resample $t1w_t2w_ratio_ribbon $transformation_file $mni_template ENCLOSING_VOXEL $subj_dir/$subj/T1w/T1wDividedByT2w_ribbon_MNI152.nii.gz -fnirt $t1w_acpc_dc
done
