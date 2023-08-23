function [] = convert_h5_voxel_based_structural_connectivity_to_mat(subject,dataset,h5_file_ID)
% The purpose of this script is to look into the sub-*_edge_weights_n*.h5
% files generated for each subject (i.e., the voxel-based structural
% connectivity matrices) and generate the .mat files corresponding to the
% variable 'mean_scores' to be used later to calculate structure-function coupling.
%
% 
% Requirements:
% The Python script create_voxel_conn_matrix.py has been run, which 
% generates unthresholded connectivity matrices corresponding to the subject 
% of interest (in parts; .h5 format).
%
% Inputs:
% 1. subject -> The subject ID of interest
% 2. dataset -> The dataset ID of interest; this was mainly for my benefit 
%               since I had 2 different datasets that I was analyzing in different directories; 
%               feel free to remove this input after adjusting the paths 
%				below to fit your data.
% 3. h5_file_ID -> The file ID of each .h5 file being considered (check the
%                  Python script create_voxel_conn_matrix.py for more
%                  details).
%
%
% Outputs:
% Creates and saves a separate .mat file for each file_ID and for the defined
% edge weight (here, "mean_scores") containing the subject's voxel-based structural 
% connectivity data for each voxel. 
% Each mat file contains (i) a region_IDs variable which is essentially the 
% cortical voxel's ID and (ii) a data variable which is of size 
% 1000 x <number of total cortical voxels corresponding to the subject processed>
% and contains the actualy connectivity values.
%
% Author: Panagiotis Fotiadis, 2021-2023
%
% If parts of the following code are used, please cite the following paper:
%
% Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD, 
% Shinohara RT, and Bassett DS “Myelination and excitation-inhibition 
% balance synergistically shape structure-function coupling across the 
% human cortex,” (2023) Nature Communications.


% -------------------------------------------------------------------------

% Pass the values from the h5 file to the connectivity matrix

% Due to the labels not necessarily being consecutive numbers, I will
% set j below to loop through the Region_Labels component of the h5
% file
region_labels = h5read(fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset), 'Nifti', subject, 'mittens_results', string(subject)+'_edge_weights_n'+int2str(h5_file_ID)+'.h5'), '/Region_Labels');

% Initialize output variables
region_ids = [];
data = [];

for j = 1:length(region_labels)
    file_path = ['/', int2str(region_labels(j)), '.0/mean_scores'];
    
    region_ids(j,1) = region_labels(j); % corresponds to the actual region labels
    data(j,:) = h5read(fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset), 'Nifti', subject, 'mittens_results', string(subject)+'_edge_weights_n'+int2str(h5_file_ID)+'.h5'), file_path)'; % corresponds to the connectivity data
end

% Save the variable of interest for each respective .h5 file into a .mat file
save_path = ['/cbica/home/panosf/Q7/Q7_Dataset_' int2str(dataset) '/Nifti/' subject '/mittens_results/' subject '_edge_weights_mean_scores_n' int2str(h5_file_ID) '.mat'];
save(save_path, 'region_ids', 'data');