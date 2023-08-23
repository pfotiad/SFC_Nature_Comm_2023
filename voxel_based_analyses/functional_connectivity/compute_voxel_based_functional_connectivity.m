function [] = compute_voxel_based_functional_connectivity(subject,part)

% This script is meant to be the continuation of the
% extract_voxel_based_BOLD_signal.m script, and its purpose is 
% to take the BOLD output .mat file and output the subject's voxel-based 
% functional connectivity matrix as the Pearson's correlation of its signal
% time series.
% Due to memory limitations, I will do this step in parts.

% Requirements:
% 1. The script extract_voxel_based_BOLD_signal.m has been run.

% Inputs:
% subject -> subject ID
% part    -> part ID ranging from 1:ceil(size(BOLD,1)/1000) for now

% Output: 
% 1. The voxel-based functional connectivity matrix of the subject, saved
% in parts.
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

% Load the BOLD variable from the previous script
BOLD = load(char(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_BOLD_timeseries_matching_structural_connectivity_voxel_IDs.mat')));
BOLD = BOLD.BOLD;

% Set start and end points
k = str2num(part);

start_ID = (k-1)*1000 + 1;

if (k < ceil(size(BOLD,1)/1000))
    end_ID = 1000 * k;
else
    end_ID = size(BOLD,1);
end

% Calculate the functional connectivity
FC = [];
ii = 1;

for i=start_ID:end_ID
    for j=1:size(BOLD,1)
        FC_tmp = [];
        FC_tmp = corrcoef(BOLD(i,2:end)',BOLD(j,2:end)');
        FC(ii,j) = FC_tmp(2,1); % have ii go from 1:1000
    end
    ii = ii + 1;
end
    
% Save the corresponding matrix
save(char(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_functional_connectivity_matrix_part_'+string(k)+'.mat')),'FC');
