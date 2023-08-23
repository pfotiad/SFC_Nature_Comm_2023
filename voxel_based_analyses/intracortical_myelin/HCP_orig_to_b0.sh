#!/bin/bash

## The purpose of this script is to register the HCP_Pipeline derived orig.mgz file of each subject into the subject's b0 diffusion space

## PREREQUISITES:

# 1. Process the T1- and T2-weighted sequences of the Penn subjects using the Human Connectome Project (HCP) pre-processing pipelines: (https://github.com/Washington-University/HCPpipelines/wiki/Installation-and-Usage-Instructions#running-the-hcp-pipelines-on-example-data,
# and as described in the manuscript; these scripts will generate the T1-weighted/T2-weighted ‘myelin maps’ for each subject as a proxy for their intracortical myelin content,

## RUNNING THE SCRIPT:

# Adapt all paths below to the ones corresponding to your data
# In order to run the script, type the script's name into the terminal

## OUTPUT:

# 1. Each subject's orig.mgz file registered into the same subject's native diffusion (b0) space

# Author: Panagiotis Fotiadis, 2021-2023

# If parts of the following code are used, please cite the following paper:

# Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD,
# Shinohara RT, and Bassett DS “Myelination and excitation-inhibition
# balance synergistically shape structure-function coupling across the
# human cortex,” (2023) Nature Communications.


## -------------------------------------------------------------------------------------------------------------------------------------------

# Define the subject directory
subj_dir=/cbica/home/panosf/Q7/structure_function_coupling/HCP_Processing_Pipeline

# Define the subjects file
subj_file=$subj_dir/subject_list.txt
count=$(cat $subj_file | wc -l)

# Loop through the subjects
for (( i=1; i<=$count; i++ ));
do
	subj=$(cat $subj_file | head -$i | tail -1 | awk '{print $1}')

	# Define SUBJECTS_DIR
	export SUBJECTS_DIR=$subj_dir/$subj/T1w

  # Locate the b0 file exported from QSIprep
  b0=$(find /cbica/home/panosf/Q7/qsiprep_results/${subj}_qsiprep/ -name "${subj}_ses-1_acq-Q7_space-T1w_dwiref.nii.gz")

  if [[ -e $b0 ]]; then
    # Register the b0 to orig
    bbregister --s $subj --mov $b0 --reg $SUBJECTS_DIR/$subj/mri/b0_to_orig.dat --dti --init-fsl

    # Apply the inverse of the above matrix to register the orig.mgz into the b0 space
    mri_vol2vol --mov $b0 --targ $SUBJECTS_DIR/$subj/mri/orig.mgz --inv --interp nearest --o $SUBJECTS_DIR/$subj/mri/orig_to_b0.nii.gz --reg $SUBJECTS_DIR/$subj/mri/b0_to_orig.dat --no-save-reg
  else
    echo "b0 doesn't exist for subject $subj"
  fi
done
