function [temporal_FC] = concatenate_voxel_based_temporal_functional_connectivity(subject,part)

% The goal of this script it to concatenate the information from the
% temporally-contiguous functional connectivity matrices into one 3d matrix 
% corresponding to each part (parts defined as in the case of the static
% voxel-based unctional connectivity matrices)

% Requirements:
% 1. The script compute_voxel_based_temporal_functional_connectivity.m has been run.
%
% Inputs:
% 1. subject -> Subject ID
% 2. part    -> part ID ranging from 1:ceil(size(BOLD,1)/1000)

% Output: 
% 1. temporal_FC -> The temporally-contiguous voxel-based functional 
% connectivity matrix of the subject, concatenated across windows and saved in parts.
% 2. the variable above saved into a mat file.
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

% Initialize variables
temporal_FC = [];
data = {};

% Concatenate
for window = 1:20
    data{window} = load(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_temporal_functional_connectivity_matrix_part_'+string(part)+'_window_'+string(window)+'.mat'));
    data{window} = data{window}.temporal_FC;
    
    temporal_FC(:,:,window) = data{window};
end

% save the variable of interest
save(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_temporal_functional_connectivity_matrix_part_'+string(part)+'.mat'),'temporal_FC','-v7.3');
