#!/bin/bash

# The purpose of this script is to generate the structural connectivity matrix of a given Human Connectome Project (HCP) subject, using two brain parcellations:
# the Schaefer 400 parcellation (400 cortical regions) and the HCP multi-modal (aka Glasser) parcellation (360 cortical regions).
# Ths script will then perform probabilistic tractography on that pre-defined subject, using the MRtrix3 toolbox.


## PREREQUISITES:

# Before running the script you need to download:
# 1. The publicly available 100 unrelated HCP dataset distribution: https://db.humanconnectome.org/.
# 2. The connectome_workbench module provided by the HCP: https://www.humanconnectome.org/software/connectome-workbench
# 3. The Mrtrix3 software: https://www.mrtrix.org

## RUNNING THE SCRIPT:

# Adapt all paths below to the ones corresponding to your data
# In order to run the script, type into the terminal: structural_connectome_generation.sh <ID of the subject you would like to process>

## OUTPUT:

# Structural connectivity matrix in .csv format (output file names: Schaefer2018_400Parcels_7Networks_order_connectome.csv and Glasser_cortical_connectome.csv)

# Author: Panagiotis Fotiadis, 2021-2023

# If parts of the following code are used, please cite the following paper:

# Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD,
# Shinohara RT, and Bassett DS “Myelination and excitation-inhibition
# balance synergistically shape structure-function coupling across the
# human cortex,” (2023) Nature Communications.


## -------------------------------------------------------------------------------------------------------------------------------------------

# Load the connectome_workbench module
module load connectome_workbench/1.4.2

# Define the subject
subj=$1

# Define the proper directories housing all subjects (change depending on the location of your subjects)
hcp_dir=/cbica/projects/HCP_Data_Releases/HCP_1200
subj_dir=/cbica/home/panosf/HCP_Dataset

# Define other pertinent directories
diff_dir=$hcp_dir/$subj/T1w/Diffusion

[[ ! -e $subj_dir/$subj/structural_connectivity ]] && mkdir -p $subj_dir/$subj/structural_connectivity
output_dir=$subj_dir/$subj/structural_connectivity

echo -e "Processing subject $subj on" $(date)

## -------------------------------------------------------------------------------------------------------------------------------------------
## DIFFUSION IMAGE PROCESSING

# Convert the HCP-derived preprocessed diffusion data to MRtrix format
mrconvert $diff_dir/data.nii.gz $output_dir/DWI.mif -fslgrad $diff_dir/bvecs $diff_dir/bvals -datatype float32 -strides 0,0,0,1 -quiet -force

# Generate a mean b=0 image for visualization purposes
dwiextract $output_dir/DWI.mif - -bzero | mrmath - mean $output_dir/meanb0.mif -axis 3

# Generate a basis response function, specific to the subject's data
dwi2response dhollander $output_dir/DWI.mif $output_dir/RF_WM_dhollander.txt $output_dir/RF_GM_dhollander.txt $output_dir/RF_CSF_dhollander.txt -voxels $output_dir/RF_voxels_dhollander.mif -force

# Create Fiber Orientation Densities (i.e., perform Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution)
dwi2fod msmt_csd $output_dir/DWI.mif $output_dir/RF_WM_dhollander.txt $output_dir/RF_WM_FOD_dhollander.mif $output_dir/RF_GM_dhollander.txt $output_dir/RF_GM_FOD_dhollander.mif $output_dir/RF_CSF_dhollander.txt $output_dir/RF_CSF_FOD_dhollander.mif -mask $diff_dir/nodif_brain_mask.nii.gz -force

# Combine the above FODs into one image, for visualization purposes
mrconvert $output_dir/RF_WM_FOD_dhollander.mif - -coord 3 0 | mrcat $output_dir/RF_CSF_FOD_dhollander.mif $output_dir/RF_GM_FOD_dhollander.mif - $output_dir/combined_FODs.mif -axis 3

# Keep track of time
echo -e "Diffusion image processing completed on" $(date)

## -------------------------------------------------------------------------------------------------------------------------------------------
## STRUCTURAL IMAGE PROCESSING

# Segment the anatomical image into tissue types, appropriate for Anatomically-Constrained Tractography
SGE_ROOT="" 5ttgen fsl $hcp_dir/$subj/T1w/T1w_acpc_dc_restore_brain.nii.gz $output_dir/5TT_in_T1w_space.mif -premasked -force

# Collapse the multi-tissue image into a 3D greyscale image for visualisation
5tt2vis $output_dir/5TT_in_T1w_space.mif $output_dir/5TT_vis.mif

# Keep track of time
echo -e "Structural image processing completed on" $(date)

## -------------------------------------------------------------------------------------------------------------------------------------------
## TRACTOGRAPHY

# Generate the initial tractogram
tckgen $output_dir/RF_WM_FOD_dhollander.mif $output_dir/tracks_10M.tck -algorithm iFOD2 -act $output_dir/5TT_in_T1w_space.mif -backtrack -crop_at_gmwmi -seed_dynamic $output_dir/RF_WM_FOD_dhollander.mif -maxlength 300 -select 10M -cutoff 0.06 -force

# Refining the Streamlines with tcksift2
tcksift2 $output_dir/tracks_10M.tck $output_dir/RF_WM_FOD_dhollander.mif $output_dir/tracks_10M_sift2.txt -act $output_dir/5TT_in_T1w_space.mif -force

# Map streamlines to the parcellated image to produce a connectome

# Perform analyses for the Schaefer atlas defined
schaefer_atlas="Schaefer2018_400Parcels_7Networks_order"

# Follow the instructions in https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/project_to_individual
# to register the Schaefer atlas in individual space (since HCP data have been processed with FreeSurfer)
export SUBJECTS_DIR=$hcp_dir/$subj/T1w

for hemi in lh rh;
do
# Considering how the HCP data release doesn't have an fsaverage subject in their $SUBJECTS_DIR I'll need to be creative
  mri_surf2surf --hemi $hemi --srcsubject ../../../../../software/external/freesurfer/centos7/6.0.0/subjects/fsaverage5 --trgsubject $subj --sval-annot /cbica/home/panosf/software/Schaefer_parcellations/FreeSurfer5.3/fsaverage5/label/${hemi}.${schaefer_atlas}.annot --tval $output_dir/${hemi}.${schaefer_atlas}.annot
done

# Generate the Schaefer2018 parcellation in volume space
# Getting creative one more time
[[ ! -e $subj_dir/$subj/mri ]] && mkdir -p $subj_dir/$subj/mri
[[ ! -e $subj_dir/$subj/surf ]] && mkdir -p $subj_dir/$subj/surf
[[ ! -e $subj_dir/$subj/label ]] && mkdir -p $subj_dir/$subj/label

mv $output_dir/*h.${schaefer_atlas}.annot $subj_dir/$subj/label/

cp -f $hcp_dir/$subj/T1w/$subj/surf/*h.pial $hcp_dir/$subj/T1w/$subj/surf/*h.white $subj_dir/$subj/surf/
cp -f $hcp_dir/$subj/T1w/$subj/mri/ribbon.mgz $hcp_dir/$subj/T1w/$subj/mri/aseg.mgz $subj_dir/$subj/mri/

export SUBJECTS_DIR=$subj_dir

mri_aparc2aseg --s $subj --o $output_dir/${schaefer_atlas}.mgz --annot $schaefer_atlas

# Convert the labels of the atlas parcellation to a format that MRtrix understands
labelconvert $output_dir/${schaefer_atlas}.mgz /cbica/home/panosf/software/Schaefer_parcellations/project_to_individual/${schaefer_atlas}_LUT.txt /cbica/home/panosf/miniconda3/pkgs/mrtrix3-3.0.2-h6bb024c_0/share/mrtrix3/labelconvert/schaefer_atlases/${schaefer_atlas}.txt $output_dir/${schaefer_atlas}_parcels.mif -force

# Perform the tractography analysis
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in $output_dir/tracks_10M_sift2.txt $output_dir/tracks_10M.tck $output_dir/${schaefer_atlas}_parcels.mif $output_dir/${schaefer_atlas}_connectome.csv -out_assignment $output_dir/${schaefer_atlas}_connectome_assignments.csv -force

## -------------------------------------------------------------------------------------------------------------------------------------------

# Perform the same analysis using the Glasser atlas
# Specifically follow the directions in: https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446?file=5528837
# This method appears to be the optimal method to register the Glasser atlas to fsaverage (and then to subject space).

# Get hemispheric gii versions of the dlabel Glasser atlas
wb_command -cifti-separate /cbica/home/panosf/software/Glasser_parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii COLUMN -label CORTEX_LEFT $subj_dir/$subj/label/L.Glasser_cortical.label.gii -label CORTEX_RIGHT $subj_dir/$subj/label/R.Glasser_cortical.label.gii

for hemi in L R;
do
  # Resample the Glasser label files from fs_32k space to fsaverage space
  wb_command -label-resample $subj_dir/$subj/label/${hemi}.Glasser_cortical.label.gii /cbica/home/panosf/software/HCPpipelines-4.3.0/global/templates/standard_mesh_atlases/${hemi}.sphere.32k_fs_LR.surf.gii /cbica/home/panosf/software/HCPpipelines-4.3.0/global/templates/standard_mesh_atlases/fs_${hemi}/fs_${hemi}-to-fs_LR_fsaverage.${hemi}_LR.spherical_std.164k_fs_${hemi}.surf.gii BARYCENTRIC $subj_dir/$subj/label/${hemi}.Glasser_cortical.fsaverage164.label.gii

  hemi2=$(echo "${hemi}"h | tr '[:upper:]' '[:lower:]')

  # Using FreeSurfer, convert the gii files to annot files
  mris_convert --annot $subj_dir/$subj/label/${hemi}.Glasser_cortical.fsaverage164.label.gii /cbica/home/panosf/software/HCPpipelines-4.3.0/global/templates/standard_mesh_atlases/fs_${hemi}/fs_${hemi}-to-fs_LR_fsaverage.${hemi}_LR.spherical_std.164k_fs_${hemi}.surf.gii $subj_dir/$subj/label/${hemi2}.Glasser_cortical_fsaverage.annot

  # Now perform the mri_surf2surf command as in the Schaefer case above, to send the annot files from fsaverage to subject space
  export SUBJECTS_DIR=$hcp_dir/$subj/T1w

  mri_surf2surf --hemi ${hemi2} --srcsubject ../../../../../software/external/freesurfer/centos7/6.0.0/subjects/fsaverage --trgsubject $subj --sval-annot $subj_dir/$subj/label/${hemi2}.Glasser_cortical_fsaverage.annot --tval $subj_dir/$subj/label/${hemi2}.Glasser_cortical_tmp.annot
done

# Before generating the volumes, I'll need to perform some MATLAB magic where I alter the ordering of the labels: The way that the above commands save label files is by also including the corresponding label values for each ROI.
# It would be great if each hemispheric label file contained only the label IDs of its own hemisphere, but for some reason each hemisphere contains label IDs in its LUT files for both hemispheres, which is messing up the aparc2aseg command.
# This can be checked with the mris_info command on an .annot file or the wb_command -file-information on any .gii label file.
# Therefore what I'm doing here is loading the *h.Glasser*_tmp.annot file and removing all labels that correspond to the other hemisphere, and resaving the .annot file (as a new file)
matlab -nodisplay -nosplash -r \
  "[v_R, L_R, ct_R] = read_annotation('$subj_dir/$subj/label/rh.Glasser_cortical_tmp.annot'); \
  ct_R.struct_names(182:end) = []; \
  ct_R.table(182:end,:) = []; \
  ct_R.numEntries = 181; \
  write_annotation('$subj_dir/$subj/label/rh.Glasser_cortical.annot',v_R,L_R,ct_R); \
  [v_L, L_L, ct_L] = read_annotation('$subj_dir/$subj/label/lh.Glasser_cortical_tmp.annot'); \
  ct_L.table(2:181,:) = []; \
  ct_L.struct_names(2:181) = []; \
  ct_L.numEntries = 181; \
  write_annotation('$subj_dir/$subj/label/lh.Glasser_cortical.annot',v_L,L_L,ct_L); exit"

# Generate the corresponding parcellation volume in subject space
export SUBJECTS_DIR=$subj_dir

mri_aparc2aseg --s $subj --o $output_dir/Glasser_atlas.mgz --annot Glasser_cortical

# Convert the labels of the atlas parcellation to a format that MRtrix understands
labelconvert $output_dir/Glasser_atlas.mgz /cbica/home/panosf/miniconda3/pkgs/mrtrix3-3.0.2-h6bb024c_0/share/mrtrix3/labelconvert/glasser_atlases/Glasser_cortical_and_subcortical_aparc2aseg_format.txt /cbica/home/panosf/miniconda3/pkgs/mrtrix3-3.0.2-h6bb024c_0/share/mrtrix3/labelconvert/glasser_atlases/Glasser_cortical.txt $output_dir/Glasser_atlas_cortical_parcels.mif -force

# Perform the actual tractography, as above
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in $output_dir/tracks_10M_sift2.txt $output_dir/tracks_10M.tck $output_dir/Glasser_atlas_cortical_parcels.mif $output_dir/Glasser_cortical_connectome.csv -out_assignment $output_dir/Glasser_cortical_connectome_assignments.csv -force

# Keep track of time
echo -e "Tractography completed on" $(date)
echo -e "Subject ${subj} is all done!"

## -------------------------------------------------------------------------------------------------------------------------------------------
