function [] = compute_voxel_based_temporal_functional_connectivity(subject,part,window)

% This script is meant to be the continuation of the
% extract_voxel_based_BOLD_signal.m and compute_voxel_based_functional_connectivity.m 
% scripts and it will generate voxel-based temporal functional connectivity 
% matrices for the designated time windows.

% Requirements:
% 1. The script extract_voxel_based_BOLD_signal.m has been run.
% 2. The script compute_voxel_based_functional_connectivity.m has been run.

% Inputs:
% subject -> subject ID
% part    -> part ID ranging from 1:ceil(size(BOLD,1)/1000) for now
% window  -> window ID ranges from 1:20 (unless I change the number of time
% windows)

% Output: 
% 1. The temporally-contiguous voxel-based functional connectivity matrix 
% of the subject, saved in parts.
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


% Load the BOLD timeseries
BOLD = load(char(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_BOLD_timeseries_matching_structural_connectivity_voxel_IDs.mat')));
BOLD = BOLD.BOLD;

% Remove the first column of the BOLD variable for simplicity
BOLD(:,1) = [];

% Define the pertinent sliding window-related parameters
window_duration = 120; % units: acquisition scans
num_windows = size(BOLD,2)/window_duration;
overlap = 0; % let this as 0 for now, units: acquisition scans

% Set start and end points
k = str2num(part);

start_ID = (k-1)*1000 + 1;

if (k < ceil(size(BOLD,1)/1000))
    end_ID = 1000 * k;
else
    end_ID = size(BOLD,1);
end

% Calculate temporal functional connectivity
temporal_FC = [];
ii = 1;

for i=start_ID:end_ID
    for j=1:size(BOLD,1)
        
        % Compute the Pearson's correlation coefficients between each pair of regions of interest (ROIs) for each time-window
        % Define the acquisition scans corresponding to each window
        window_indices = (window-1)*window_duration+1:window*window_duration; % units: acquisition scans
        
        temporal_FC_tmp = [];
        temporal_FC_tmp = corrcoef(BOLD(i,window_indices)',BOLD(j,window_indices)'); % have ii go from 1:1000
        temporal_FC(ii,j) = temporal_FC_tmp(2,1);
    end
    
    ii = ii + 1;
end

% Save the corresponding matrix
save(char(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_temporal_functional_connectivity_matrix_part_'+string(k)+'_window_'+string(window)+'.mat')),'temporal_FC');



