#!/usr/bin/python

# This script will create voxel-based structural connectivity matrices for a given subject using the MITTENS analytic tractography pipeline
# This script was used to analyze the subjects that were scanned here at the University of Pennsylvania

## PREREQUISITES:

# Before running the script you need to:
# 1. Diffusion imaging data of the subject of interest have been processed using QSIPrep and QSIrecon: https://qsiprep.readthedocs.io/en/latest/
# 2. Download the MITTENS Python library: https://github.com/mattcieslak/MITTENS

## RUNNING THE SCRIPT:

# Adapt all paths below to the ones corresponding to your data
# In order to run the script, type into the terminal: python create_voxel_conn_matrix.py <subject ID> <dataset ID> <file ID>

## INPUT:

# 1. subject ID -> Refers to the ID of the subject you would like to run
# 2. dataset ID -> This was mainly for my benefit since I had 2 different datasets that I was analyzing in different directories; feel free to remove this input after adjusting the
#				paths below to fit your data
# 3. file ID -> Given the sheer magnitude of the number of voxels and the intense computational power that would be required to run this script for all 60,000 - 90,000 voxels, I
#				am splitting the analysis and running it in parts. So, for instance, if you define file_ID = 1 here then the analysis will be run for voxels indexed 0 to 999, then
#				file_ID = 2 would run voxels indexed 1000 to 1999 etc. up until all voxels have been processed.

## OUTPUT:

# An .h5 (nested) file called <subject_ID>_edge_weights_n<file_ID>.h5 that contains six different types of edge weights ('raw_probs', 'null_scores', 'path_lengths', 'prob_ratio', 'mean_scores', 'mean_null_scores')
# whose detailed explanations can be found in the MITTENS documentation: https://github.com/mattcieslak/MITTENS.
# Each type of edge weight entry is of size 1000 x <number of total cortical voxels corresponding to the subject processed>
# In the follow-up scripts used to define structural connectivity, I am only focusing on 'mean_scores' as the edge weight of interest.

# Author: Panagiotis Fotiadis, 2021-2023

# If parts of the following code are used, please cite the following paper:

# Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD,
# Shinohara RT, and Bassett DS “Myelination and excitation-inhibition
# balance synergistically shape structure-function coupling across the
# human cortex,” (2023) Nature Communications.


## -------------------------------------------------------------------------------------------------------------------------------------------


import numpy as np
import math
from mittens import MITTENS
from mittens import VoxelGraph
import networkit
from networkit.algebraic import adjacencyMatrix
import os
import h5py
from tqdm import tqdm
import sys
import scipy.io as sio

# Define the subject as the first input, the dataset_ID as the second input, and the file_ID as the third input
subject = sys.argv[1]
dataset_ID = sys.argv[2]
file_ID = sys.argv[3]

print('Processing subject: '+subject+'. File ID: '+file_ID)

# Convert the start and end points into integers
file_ID = int(file_ID)

# If MITTENS hasn't been run, run it
# Check to see if the voxel null graph exists
tmp_output_save_null_file = subject+'_doubleODF_voxel_null_graph.mat'
tmp_output_save_null_path = os.path.join('/cbica/home/panosf/Q7/Q7_Dataset_'+dataset_ID, 'Nifti/', subject, 'mittens_results/', tmp_output_save_null_file)

if not os.path.exists(tmp_output_save_null_path):

	# Presumably, QSIPrep and QSIrecon have already been run
	reconstruction_file_name = subject+'_ses-1_acq-Q7_space-T1w_desc-preproc_space-T1w_gqi.fib.gz'
	reconstruction_file_path = os.path.join('/cbica/home/panosf/Q7/Q7_Dataset_'+dataset_ID, 'Nifti/', subject, 'mittens_results/', reconstruction_file_name)

	real_affine_image_file_name = subject+'_ses-1_acq-Q7_space-T1w_dwiref.nii.gz'
	real_affine_image_file_path = os.path.join('/cbica/home/panosf/Q7/Q7_Dataset_'+dataset_ID, 'Nifti/', subject, 'mittens_results/', real_affine_image_file_name)

	mitns = MITTENS(reconstruction=reconstruction_file_path, real_affine_image=real_affine_image_file_path)

	# Calculate inter-voxel fiber transition probabilities & write out NIfTI images for the transition probabilities to each neighbor using both the singleODF and doubleODF methods
	mitns.calculate_transition_probabilities(output_prefix=subject)

	# Create and save a Voxel Graph, using networkit, where edges are weighted by the tract transition expectation from one voxel to another.
	# In essence, transition probabilities are converted to edge weights
	voxel_graph = mitns.build_graph(doubleODF=True, weighting_scheme="negative_log_p")

	output_save_file = subject+'_doubleODF_voxel_graph.mat'
	output_save_path = os.path.join('/cbica/home/panosf/Q7/Q7_Dataset_'+dataset_ID, 'Nifti/', subject, 'mittens_results/', output_save_file)

	voxel_graph.save(output_save_path)

	# Create and save a null graph (i.e., a matching voxel graph using probabilities from isotropic ODFs)
	null_graph = mitns.build_null_graph(doubleODF=True, purpose="shortest paths")

	output_save_null_file = subject+'_doubleODF_voxel_null_graph.mat'
	output_save_null_path = os.path.join('/cbica/home/panosf/Q7/Q7_Dataset_'+dataset_ID, 'Nifti/', subject, 'mittens_results/', output_save_null_file)

	null_graph.save(output_save_null_path)

#-----------------------------------------------------------------------------------------------------------------------------#

# Load both regular and null voxel graphs as well as an atlas

# Voxel graph path
output_save_file = subject+'_doubleODF_voxel_graph.mat'
output_save_path = os.path.join('/cbica/home/panosf/Q7/Q7_Dataset_'+dataset_ID, 'Nifti/', subject, 'mittens_results/', output_save_file)

# Null graph path
output_save_null_file = subject+'_doubleODF_voxel_null_graph.mat'
output_save_null_path = os.path.join('/cbica/home/panosf/Q7/Q7_Dataset_'+dataset_ID, 'Nifti/', subject, 'mittens_results/', output_save_null_file)

# Load them up
voxel_graph = VoxelGraph(output_save_path)
null_graph = VoxelGraph(output_save_null_path)

voxel_graph.add_null_graph(null_graph)

# Load the altas of interest
atlas_file_name = subject+'_ses-1_acq-Q7_space-T1w_desc-preproc_space-T1w_desc-schaefer100x7_atlas_cortex_voxel_based.nii.gz'
atlas_file_path = os.path.join('/cbica/home/panosf/Q7/Q7_Dataset_'+dataset_ID, 'Nifti/', subject, 'mittens_results/', atlas_file_name)

voxel_graph.add_atlas(atlas_file_path)

#-----------------------------------------------------------------------------------------------------------------------------#

# Set up the label set
if voxel_graph.atlas_labels is None or not len(voxel_graph.atlas_labels) == voxel_graph.nvoxels:
	full_region_ids = np.arange(voxel_graph.nvoxels)
else:
	full_region_ids = np.unique(voxel_graph.atlas_labels[voxel_graph.atlas_labels > 0])

#-----------------------------------------------------------------------------------------------------------------------------#

# Nullify the output variable of interest
edge_weights = {}

# Define the start and end points for the for loop to run.
num_voxels = 1000 # Denote the number of voxels per file

start_point = (file_ID - 1) * num_voxels

if file_ID == math.floor(len(full_region_ids) / num_voxels) + 1:
	end_point = len(full_region_ids)
else:
	end_point = file_ID * num_voxels

# Run the connectivity for loop
for source_region in full_region_ids[start_point:end_point]:
    print(source_region)

    raw_scores_weight, null_scores_weight, path_lengths_weight, prob_ratio_weight, mean_scores_weight, mean_null_scores_weight = voxel_graph.voxel_connectivity_map(source_region)

    edge_weights[source_region] = {'raw_probs': raw_scores_weight, 'null_scores': null_scores_weight, 'path_lengths': path_lengths_weight, 'prob_ratio': prob_ratio_weight, 'mean_scores': mean_scores_weight, 'mean_null_scores': mean_null_scores_weight}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Save the edge_weights variable using .h5 format

edge_weights_path_h5_name = '/cbica/home/panosf/Q7/Q7_Dataset_%s/Nifti/%s/mittens_results/%s_edge_weights_n%s.h5' % (dataset_ID, subject, subject, file_ID)

hf = h5py.File(edge_weights_path_h5_name, 'w')

for i, v in edge_weights.items():
	hf.create_dataset(str(i)+'/raw_probs', data=list(edge_weights[i]['raw_probs'].values()))
	hf.create_dataset(str(i)+'/null_scores', data=list(edge_weights[i]['null_scores'].values()))
	hf.create_dataset(str(i)+'/path_lengths', data=list(edge_weights[i]['path_lengths'].values()))
	hf.create_dataset(str(i)+'/prob_ratio', data=list(edge_weights[i]['prob_ratio'].values()))
	hf.create_dataset(str(i)+'/mean_scores', data=list(edge_weights[i]['mean_scores'].values()))
	hf.create_dataset(str(i)+'/mean_null_scores', data=list(edge_weights[i]['mean_null_scores'].values()))

hf.create_dataset('Region_Labels', data=list(edge_weights.keys()))

hf.close()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

print('Done!')
