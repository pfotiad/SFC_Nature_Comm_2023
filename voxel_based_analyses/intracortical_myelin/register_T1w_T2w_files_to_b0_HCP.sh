#!/bin/bash

# The purpose of this script is to register the T1w/T2w ratio files obtained from the Human Connectome Project (HCP) Pipeline
# (located under /cbica/home/panosf/Q7/structure_function_coupling/HCP_Processing_Pipeline/$subj/T1w/*) from MNI -> b0 space.

## PREREQUISITES:

# 1. Process the T1- and T2-weighted sequences of the Penn subjects using the HCP pre-processing pipelines: (https://github.com/Washington-University/HCPpipelines/wiki/Installation-and-Usage-Instructions#running-the-hcp-pipelines-on-example-data,
# and as described in the manuscript; these scripts will generate the T1-weighted/T2-weighted ‘myelin maps’ for each subject as a proxy for their intracortical myelin content,
# 2. Run the script HCP_orig_to_b0.sh

## RUNNING THE SCRIPT:

# Adapt all paths below to the ones corresponding to your data
# In order to run the script, type the script's name into the terminal

## OUTPUT:

# 1. Each subject's T1w/T2w ratio map registered into the same subject's native diffusion (b0) space

# Author: Panagiotis Fotiadis, 2021-2023

# If parts of the following code are used, please cite the following paper:

# Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD,
# Shinohara RT, and Bassett DS “Myelination and excitation-inhibition
# balance synergistically shape structure-function coupling across the
# human cortex,” (2023) Nature Communications.


## -------------------------------------------------------------------------------------------------------------------------------------------


# Subjects list
subject_list=/cbica/home/panosf/Q7/structure_function_coupling/subject_list.txt

# Count the subjects
count=$(cat $subject_list | wc -l)

# Define the HCP directory
hcp_dir=/cbica/home/panosf/Q7/structure_function_coupling/HCP_Processing_Pipeline

# Go through each subject
for (( i=1; i<=$count; i++ ));
do
  subj=$(cat $subject_list | head -$i | tail -1 | awk '{print $1}')

  echo -e "Processing subject ${subj}: " $(date)

  # Locate the b0 file and the files of interest
  T1w_restore=$hcp_dir/$subj/T1w/T1w_acpc_dc_restore.nii.gz
  T2w_restore=$hcp_dir/$subj/T1w/T2w_acpc_dc_restore.nii.gz

  b0=$(find /cbica/home/panosf/Q7/qsiprep_results/${subj}_qsiprep/ -name "${subj}_ses-1_acq-Q7_space-T1w_dwiref.nii.gz")

# -----------------------------------------------------------------------------------------------------------------------------------------------------

  # Register each of the above files from MNI to b0 space

  # First, bring them to orig space (since the ribbon file is in orig space and not in 001.mgz space)
  mri_convert $T1w_restore $hcp_dir/$subj/T1w/T1w_acpc_dc_restore_orig.nii.gz --conform
  mri_convert $T2w_restore $hcp_dir/$subj/T1w/T2w_acpc_dc_restore_orig.nii.gz --conform

  # Then:
  [[ ! -e $hcp_dir/$subj/T1w/$subj/mri/ribbon.nii.gz ]] && mri_convert $hcp_dir/$subj/T1w/$subj/mri/ribbon.mgz $hcp_dir/$subj/T1w/$subj/mri/ribbon.nii.gz
  wb_command -volume-math "(T1w / T2w) * (((ribbon > (3 - 0.01)) * (ribbon < (3 + 0.01))) + ((ribbon > (42 - 0.01)) * (ribbon < (42 + 0.01))))" $hcp_dir/$subj/T1w/T1wDividedByT2w_ribbon_orig.nii.gz -var T1w $hcp_dir/$subj/T1w/T1w_acpc_dc_restore_orig.nii.gz -var T2w $hcp_dir/$subj/T1w/T2w_acpc_dc_restore_orig.nii.gz -var ribbon $hcp_dir/$subj/T1w/$subj/mri/ribbon.nii.gz
  wb_command -volume-palette $hcp_dir/$subj/T1w/T1wDividedByT2w_ribbon_orig.nii.gz MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false

  # Finally, co-register the T1w/T2w ratio volumes from orig space to b0 space. The lines below assume that you have ran HCP_orig_to_b0.sh
  if [[ -e $hcp_dir/$subj/T1w/$subj/mri/b0_to_orig.dat ]]; then
    export SUBJECTS_DIR=$hcp_dir/$subj/T1w

    mri_vol2vol --mov $b0 --targ $hcp_dir/$subj/T1w/T1wDividedByT2w_ribbon_orig.nii.gz --inv --interp nearest --o $hcp_dir/$subj/T1w/T1wDividedByT2w_ribbon_b0.nii.gz --reg $hcp_dir/$subj/T1w/$subj/mri/b0_to_orig.dat --no-save-reg
  else
    echo "You need to run the script HCP_orig_to_b0.sh first."
  fi

  echo -e "Subject $subj is all done: " $(date)
done
