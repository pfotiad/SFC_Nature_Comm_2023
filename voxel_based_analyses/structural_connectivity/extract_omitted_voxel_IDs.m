function [omitted_voxel_IDs_list] = extract_omitted_voxel_IDs(subject,dataset)
% The purpose of this script is to come up with a list of voxel IDs that
% are being skipped over in the voxel-based structural connectivity
% matrices, and then use that list to remove those rows from the
% functional signal timeseries that I have, so that when I generate the voxel-based
% functional connectivity matrices afterwards, there is direct row-matching
% with the voxel-based structural connectivity matrices.
% 
% Requirements:
% The MATLAB script convert_h5_voxel_based_structural_connectivity_to_mat.m
% has been run.
%
% Inputs:
% 1. subject -> The subject ID of interest
% 2. dataset -> The dataset ID of interest; this was mainly for my benefit 
%               since I had 2 different datasets that I was analyzing in different directories; 
%               feel free to remove this input after adjusting the paths 
%				below to fit your data.
%
% Outputs:
% 1. omitted_voxel_IDs_list -> A list of all the omitted voxel IDs
% 2. The above variable saved as a .mat file
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

% Open each .mat file consecutively, identify the voxel IDs missing, and
% append them to a text file.

% total number of structural connectivity .mat files
count_path = dir(fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset), 'Nifti', subject, 'mittens_results', string(subject)+'_edge_weights_prob_ratio_n*.mat'));
count = size(count_path,1);

% Initialize output list
omitted_voxel_IDs_list = [];

for i=1:count
    % load the file
    structural_data = load(fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset), 'Nifti', subject, 'mittens_results', string(subject)+'_edge_weights_prob_ratio_n'+string(i)+'.mat'));
    consecutive_list = [structural_data.region_ids(1):structural_data.region_ids(end)]';

    % Check to see which consecutive numbers don't exist in my
    % structural_data.region_ids vector
    idx = ismember(consecutive_list, structural_data.region_ids);
    omitted_voxel_IDs = consecutive_list(~idx);

    % write to list
    omitted_voxel_IDs_list = [omitted_voxel_IDs_list; omitted_voxel_IDs];

    % also check the possibility that the first element of the file that is
    % loaded is consecutive to the last element of the previous file.
    if (i > 1)
        prior_structural_data = load(fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset), 'Nifti', subject, 'mittens_results', string(subject)+'_edge_weights_prob_ratio_n'+string(i-1)+'.mat'));

        if (structural_data.region_ids(1) ~= prior_structural_data.region_ids(end) + 1)
            difference = structural_data.region_ids(1) - prior_structural_data.region_ids(end);
            for j = 2:difference
                % write to list
                omitted_voxel_IDs_list = [omitted_voxel_IDs_list; prior_structural_data.region_ids(end) + j - 1];
            end
        end
    end
end

% NOTE: This isn't coded in this script, but you should also check the case 
% where the first element of the first file doesn't start with the index 1. 
% (I checked it manually for out data.)

% Save the voxel IDs into the subject directory
writematrix(omitted_voxel_IDs_list, fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset), 'Nifti', subject, 'mittens_results', string(subject)+'_structural_connectivity_omitted_voxel_IDs.txt'));