#!/bin/bash

# The purpose of this script is to generate (static and temporal) functional connectivity matrices on a pre-defined Human Connectome Project (HCP) subject, using the HCP wb_command
# i.e., https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ

## PREREQUISITES:

# Before running the script you need to download:
# 1. The publicly available 100 unrelated HCP dataset distribution: https://db.humanconnectome.org/.
# 2. The connectome_workbench module provided by the HCP: https://www.humanconnectome.org/software/connectome-workbench
# 3. The Schaefer parcellations: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations
# 4. The Glasser parcellation file Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii from https://balsa.wustl.edu/file/3VLx

## RUNNING THE SCRIPT:

# Adapt all paths below to the ones corresponding to your data
# In order to run the script, type into the terminal: functional_connectome_generation.sh <ID of the subject you would like to process>

## OUTPUT:

# 1. 1 pconn.nii file containing the static functional connectivity values of the subject using the Schaefer parcellation (i.e., the subject's average functional connectivity values across the resting-state fMRI scan).
# 2. 1 pconn.nii file containing the static functional connectivity values of the subject using the Glasser parcellation (i.e., the subject's average functional connectivity values across the resting-state fMRI scan).
# 3. 20 pconn.nii files containing the functional connectivity values of the subject for each one of the 20 time-windows using the Schaefer parcellation.
# 4. 20 pconn.nii files containing the functional connectivity values of the subject for each one of the 20 time-windows using the Glasser parcellation.

# Author: Panagiotis Fotiadis, 2021-2023

# If parts of the following code are used, please cite the following paper:

# Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD,
# Shinohara RT, and Bassett DS “Myelination and excitation-inhibition
# balance synergistically shape structure-function coupling across the
# human cortex,” (2023) Nature Communications.


## -------------------------------------------------------------------------------------------------------------------------------------------

# Load the connectome_workbench module
module load connectome_workbench/1.4.2

# Define subject ID
subj=$1

# Define pertinent directories
subj_dir=/cbica/projects/HCP_Data_Releases/HCP_1200 # HCP directory housing all subjects
mni_dir=$subj_dir/$subj/MNINonLinear/Results

[[ ! -e /cbica/home/panosf/HCP_Dataset/$subj/functional_connectivity ]] && mkdir -p /cbica/home/panosf/HCP_Dataset/$subj/functional_connectivity
output_dir=/cbica/home/panosf/HCP_Dataset/$subj/functional_connectivity

# Define the atlas directory(ies)
schaefer_atlas_dir=/cbica/home/panosf/software/Schaefer_parcellations/HCP/fslr32k/cifti
glasser_atlas_dir=/cbica/home/panosf/software/Glasser_parcellations

echo -e "Processing subject $subj on" $(date)

## -------------------------------------------------------------------------------------------------------------------------------------------

# The HCP rfMRI data were collected in 4 runs:

# 1) Demean and normalize the individual timeseries of the input subject
for dtseries in rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii rfMRI_REST2_LR/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii rfMRI_REST2_RL/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii
do
  wb_command -cifti-reduce $mni_dir/${dtseries} MEAN $output_dir/mean_${dtseries:6:8}_dtseries.dscalar.nii
  wb_command -cifti-reduce $mni_dir/${dtseries} STDEV $output_dir/stdev_${dtseries:6:8}_dtseries.dscalar.nii
  wb_command -cifti-math '(x - mean) / stdev' $output_dir/norm_${dtseries:6:8}.dtseries.nii -fixnan 0 -var x $mni_dir/${dtseries} -var mean $output_dir/mean_${dtseries:6:8}_dtseries.dscalar.nii -select 1 1 -repeat -var stdev $output_dir/stdev_${dtseries:6:8}_dtseries.dscalar.nii -select 1 1 -repeat
done

# 2) Instead of concatenating the individual timeseries into a 1-hour long scan session, I will compute an "average" run (duration: ~15 min)
wb_command -cifti-average $output_dir/avg_norm.dtseries.nii -cifti $output_dir/norm_REST1_LR.dtseries.nii -cifti $output_dir/norm_REST1_RL.dtseries.nii -cifti $output_dir/norm_REST2_LR.dtseries.nii -cifti $output_dir/norm_REST2_RL.dtseries.nii

# Split the averaged timeseries into different consecutive time-windows (let's go with 20 for right now) - This step is necessary to generate temporally-contiguous functional connectivity matrices
num_windows=20

[[ ! -e $output_dir/temporal_analyses ]] && mkdir $output_dir/temporal_analyses

for (( win=1; win<=$num_windows; win++ ));
do
  total_time_points=1200 # total number of columns in the averaged time series
  duration=$(expr $total_time_points / $num_windows)
  start_column=$(expr ${duration}*${win}-${duration}+1 | bc -l)
  end_column=$(expr ${duration}*${win} | bc -l)

  wb_command -cifti-merge $output_dir/temporal_analyses/avg_norm_time_window_${win}.dtseries.nii -cifti $output_dir/avg_norm.dtseries.nii -column $start_column -up-to $end_column
done

## PERFORM STATIC FUNCTIONAL CONNECTIVITY ANALYSES

# Generate the functional connectivity matrix for the Schaefer atlas of interest
schaefer_atlas="Schaefer2018_400Parcels_7Networks_order.dlabel.nii"

# Parcellate the voxel-wise time series data into different brain regions as per the designated atlas, by averaging all the voxel-level time series belonging to each brain region,
# in order to come up with a <number of atlas ROIs> x <time-points> matrix
wb_command -cifti-parcellate $output_dir/avg_norm.dtseries.nii $schaefer_atlas_dir/$schaefer_atlas COLUMN $output_dir/${schaefer_atlas/.dlabel.nii/_timeseries}.ptseries.nii

# Take the above matrix and generate an <ROIs> x <ROIs> connectivity matrix (Pearson's correlation)
wb_command -cifti-correlation $output_dir/${schaefer_atlas/.dlabel.nii/_timeseries}.ptseries.nii $output_dir/${schaefer_atlas/.dlabel.nii/_connectome}.pconn.nii

## -------------------------------------------------------------------------------------------------------------------------------------------

# Similarly, generate the functional connectivity matrix for the Glasser atlas

# Parcellate the voxel-wise time series data into different brain regions as per the designated atlas, by averaging all the voxel-level time series belonging to each brain region,
# in order to come up with a <number of atlas ROIs> x <time-points> matrix
wb_command -cifti-parcellate $output_dir/avg_norm.dtseries.nii $glasser_atlas_dir/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii COLUMN $output_dir/Glasser_cortical_timeseries.ptseries.nii

# Take the above matrix and generate an <ROIs> x <ROIs> connectivity matrix (Pearson's correlation)
wb_command -cifti-correlation $output_dir/Glasser_cortical_timeseries.ptseries.nii $output_dir/Glasser_cortical_connectome.pconn.nii

## -------------------------------------------------------------------------------------------------------------------------------------------

## PERFORM TEMPORAL FUNCTIONAL CONNECTIVITY ANALYSES

for (( win=1; win<=$num_windows; win++ ));
do
  # Generate the temporal functional connectivity matrices for the Schaefer atlas of interest

  # Parcellate the voxel-wise time series data into different brain regions as per the designated atlas, by averaging all the voxel-level time series belonging to each brain region,
  # in order to come up with a <number of atlas ROIs> x <time-points> matrix
  wb_command -cifti-parcellate $output_dir/temporal_analyses/avg_norm_time_window_${win}.dtseries.nii $schaefer_atlas_dir/$schaefer_atlas COLUMN $output_dir/temporal_analyses/${schaefer_atlas/.dlabel.nii/_timeseries}_${win}.ptseries.nii

  # Take the above matrix and generate an <ROIs> x <ROIs> connectivity matrix (Pearson's correlation)
  wb_command -cifti-correlation $output_dir/temporal_analyses/${schaefer_atlas/.dlabel.nii/_timeseries}_${win}.ptseries.nii $output_dir/temporal_analyses/${schaefer_atlas/.dlabel.nii/_connectome}_${win}.pconn.nii

  ## -------------------------------------------------------------------------------------------------------------------------------------------

  # Generate the temporal functional connectivity matrices for the Glasser atlas

  # Parcellate the voxel-wise time series data into different brain regions as per the designated atlas, by averaging all the voxel-level time series belonging to each brain region,
  # in order to come up with a <number of atlas ROIs> x <time-points> matrix
  wb_command -cifti-parcellate $output_dir/temporal_analyses/avg_norm_time_window_${win}.dtseries.nii $glasser_atlas_dir/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii COLUMN $output_dir/temporal_analyses/Glasser_cortical_timeseries_${win}.ptseries.nii

  # Take the above matrix and generate an <ROIs> x <ROIs> connectivity matrix (Pearson's correlation)
  wb_command -cifti-correlation $output_dir/temporal_analyses/Glasser_cortical_timeseries_${win}.ptseries.nii $output_dir/temporal_analyses/Glasser_cortical_connectome_${win}.pconn.nii
done

echo -e "Subject $subj is all done: " $(date)
