#!/bin/bash

# The purpose of this script is to register the BOLD files output by the CONN software from MNI to native diffusion (b0) space and prepare the necessary files to
# create voxel-based functional connectivity matrices for each Penn subject

## PREREQUISITES:

# Before running the script you need to:
# 1. Process the subjects' resting-state functional MRI data using the CONN toolbox (https://web.conn-toolbox.org/resources/installation)

## RUNNING THE SCRIPT:

# Adapt all paths below to the ones corresponding to your data
# In order to run the script, type the name of the script into the terminal

## OUTPUT:

# 1. This script will output the pre-processed and denoised average BOLD sequence in b0 space
# 2. This script will output the pre-processed and denoised BOLD sequence for each repetition time (TR) in b0 space

# Author: Panagiotis Fotiadis, 2021-2023

# If parts of the following code are used, please cite the following paper:

# Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD,
# Shinohara RT, and Bassett DS “Myelination and excitation-inhibition
# balance synergistically shape structure-function coupling across the
# human cortex,” (2023) Nature Communications.


## -------------------------------------------------------------------------------------------------------------------------------------------


# Unload the connectome_workbench module
module unload connectome_workbench/1.4.2

# Subjects list
subject_list=/cbica/home/panosf/Q7/CONN_functional_results/list.txt

# Define subjects directory
subj_dir=/cbica/home/panosf/Q7/CONN_functional_results

# Count the subjects
count=$(cat $subject_list | wc -l)

# Go through each subject
for (( i=1; i<=$count; i++ ));
do
  subj=$(cat $subject_list | head -$i | tail -1 | awk '{print $1}')

  ## Firstly, register the denoised, preprocessed BOLD file from MNI to b0 space
  echo -e "Processing subject ${subj}: " $(date)

  # Locate the denoised, preprocessed BOLD file and the b0 file
  BOLD_MNI=$subj_dir/$subj/scans/dswau${subj}_ses-1_task-rest_bold.nii
  b0=$subj_dir/$subj/${subj}_ses-1_acq-Q7_space-T1w_dwiref.nii.gz

  # Locate the transformation file required to register files from MNI -> b0
  transform_h5=$(find /cbica/home/panosf/Q7/qsiprep_results/${subj}_qsiprep/ -name "${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5")

  # Create a relevant folder
  if [[ ! -e $subj_dir/$subj/BOLD_processing ]]; then
		mkdir -p $subj_dir/$subj/BOLD_processing/BOLD_TR_volumes
	fi

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------#
  ## STEP 1: REGISTER THE BOLD FILES FROM MNI TO B0 SPACE

  # Average the 4D BOLD volume across the 4th dimension (i.e., calculate the average time series of each voxel across all TRs)
  fslmaths $BOLD_MNI -Tmean $subj_dir/$subj/BOLD_processing/${subj}_BOLD_MNI_avg.nii.gz
  mri_convert $subj_dir/$subj/BOLD_processing/${subj}_BOLD_MNI_avg.img $subj_dir/$subj/BOLD_processing/${subj}_BOLD_MNI_avg.nii.gz
  rm -f $subj_dir/$subj/BOLD_processing/${subj}_BOLD_MNI_avg.img $subj_dir/$subj/BOLD_processing/${subj}_BOLD_MNI_avg.hdr

  # Register the avg BOLD file from MNI to b0 space
  /cbica/home/panosf/software/ANTs/ANTS-build/Examples/antsApplyTransforms --default-value 0 --dimensionality 3 --float 0 --input $subj_dir/$subj/BOLD_processing/${subj}_BOLD_MNI_avg.nii.gz --output $subj_dir/$subj/BOLD_processing/${subj}_BOLD_MNI_avg_to_b0.nii.gz --reference-image $b0 --transform $transform_h5

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------#

  # Perform the analyses above for each TR (i.e., register each TR into b0 space)

  # Break the BOLD file down into each one of its individual TRs (-> 2400 files)
	fslsplit $BOLD_MNI $subj_dir/$subj/BOLD_processing/BOLD_TR_volumes/${subj}_BOLD_TR

	for j in 000{0..9} 00{10..99} 0{100..999} {1000..2399};
	do
      mri_convert $subj_dir/$subj/BOLD_processing/BOLD_TR_volumes/${subj}_BOLD_TR${j}.img $subj_dir/$subj/BOLD_processing/BOLD_TR_volumes/${subj}_BOLD_TR${j}.nii.gz
      rm -f $subj_dir/$subj/BOLD_processing/BOLD_TR_volumes/${subj}_BOLD_TR${j}.hdr $subj_dir/$subj/BOLD_processing/BOLD_TR_volumes/${subj}_BOLD_TR${j}.img

		 # Register the BOLD TR file from MNI to b0 space
	    /cbica/home/panosf/software/ANTs/ANTS-build/Examples/antsApplyTransforms --default-value 0 --dimensionality 3 --float 0 --input $subj_dir/$subj/BOLD_processing/BOLD_TR_volumes/${subj}_BOLD_TR${j}.nii.gz --output $subj_dir/$subj/BOLD_processing/BOLD_TR_volumes/${subj}_BOLD_MNI_to_b0_TR${j}.nii.gz --reference-image $b0 --transform $transform_h5
	done

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------#

  echo -e "Subject $subj is all done: " $(date)
done
