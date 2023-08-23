#!/bin/bash

# The purpose of this script is to compute the average intracortical myelin content on the Human Connectome Project (HCP) subjects
# using the HCP multi-modal (aka Glasser) atlas-defined brain regions.

## PREREQUISITES:

# Before running the script you need to download:
# 1. The publicly available 100 unrelated HCP dataset distribution: https://db.humanconnectome.org/.
# 2. The Glasser parcellation file Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii from https://balsa.wustl.edu/file/3VLx

## RUNNING THE SCRIPT:

# Adapt all paths below to the ones corresponding to your data
# In order to run the script, type its name into the terminal

## OUTPUT:

# Text file containing 360 rows (one per brain region) for each subject, with each row containing that brain region's intracortical myelin estimate

# Author: Panagiotis Fotiadis, 2021-2023

# If parts of the following code are used, please cite the following paper:

# Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD,
# Shinohara RT, and Bassett DS “Myelination and excitation-inhibition
# balance synergistically shape structure-function coupling across the
# human cortex,” (2023) Nature Communications.


## -------------------------------------------------------------------------------------------------------------------------------------------

# Define the subject list (one subject ID per row)
subj_list=/cbica/home/panosf/HCP_Dataset/100_unrelated_subjects_list.txt

# Count the number of subjects
count=$(cat $subj_list | wc -l)

# HCP directory housing all subjects
hcp_dir=/cbica/projects/HCP_Data_Releases/HCP_1200

# Glasser atlas to be used
glasser_atlas=/cbica/home/panosf/software/Glasser_parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii

for (( i=1; i<=$count; i++ ));
do
  # Define the subject
  subj=$(cat $subj_list | head -$i | tail -1 | awk '{print $1}')

  # Define the output directory
  output_dir=/cbica/home/panosf/HCP_Dataset/$subj

  # Extract intracortical myelin content for each atlas-defined brain region
  wb_command -cifti-parcellate $hcp_dir/$subj/MNINonLinear/fsaverage_LR32k/${subj}.MyelinMap.32k_fs_LR.dscalar.nii $glasser_atlas COLUMN $output_dir/Glasser_myelin_content.pscalar.nii

  # Convert the pscalar data to .txt
  wb_command -cifti-convert -to-text $output_dir/Glasser_myelin_content.pscalar.nii $output_dir/Glasser_myelin_content.txt
done
