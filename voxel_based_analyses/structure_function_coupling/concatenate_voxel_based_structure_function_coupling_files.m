% The purpose of this script is to concatenate the thresholded 1000x1
% structure-function coupling .mat files into one file, for each subject

% Requirements:
% 1. The script extract_voxel_based_structure_function_coupling.m has been run.
%
% Output: 
% One .mat file per subject containing the subject's voxel-based structure-function 
% coupling indices for each of their cortical voxels.
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

clear; clc;

% Define the subject IDs (this should be changed depending on your dataset)
subject_IDs = ["sub-EHJ", "sub-HJK", "sub-JFG", "sub-LT", "sub-NWT", "sub-VW", "sub-YH", "sub-s10", "sub-s11"];

% Define the threshold used
threshold = 70;

% Load the files
for i = 1:numel(subject_IDs)
    
    % total number of parts
    num_parts = importdata(fullfile('/cbica/home/panosf/Q7/CONN_functional_results', subject_IDs(i), 'BOLD_processing', subject_IDs+'_number_of_parts.txt'));
    
    % Define the output variable
    coupling = {};
    
    for j = 1:num_parts
        % Load the parsed up structure-function coupling data files
        coupling_file{i}{j} = load(fullfile('/cbica/home/panosf/Q7/structure_function_coupling',subject_IDs(i),subject_IDs(i)+'_voxel_based_threshold_'+string(threshold)+'_percent_structure_function_coupling_part_'+string(j)+'.mat'));
        
        % Define indices
        start_id = (j-1) * 1000 + 1;
        
        if (j < num_parts)
            end_id = 1000 * j;
        else
            end_id = start_id + size(coupling_file{i}{j}.coupling.raw_probs.static_correlational_approach_all_data_pearson_rho,1) - 1;
        end
        
        % STATIC CASE

        % mean_scores        
        coupling.mean_scores.static_correlational_approach_non_zero_data_spearman_rho(start_id:end_id,:) = coupling_file{i}{j}.coupling.mean_scores.static_correlational_approach_non_zero_data_spearman_rho;
        coupling.mean_scores.static_correlational_approach_non_zero_data_spearman_pvalue(start_id:end_id,:) = coupling_file{i}{j}.coupling.mean_scores.static_correlational_approach_non_zero_data_spearman_pvalue;
        
        % TEMPORAL CASE

        % mean_scores      
        coupling.mean_scores.temporal_correlational_approach_non_zero_data_spearman_rho(start_id:end_id,:) = coupling_file{i}{j}.coupling.mean_scores.temporal_correlational_approach_non_zero_data_spearman_rho;
        coupling.mean_scores.temporal_correlational_approach_non_zero_data_spearman_pvalue(start_id:end_id,:) = coupling_file{i}{j}.coupling.mean_scores.temporal_correlational_approach_non_zero_data_spearman_pvalue;
    end
    
    % Save the concatenated file
    save(fullfile('/cbica/home/panosf/Q7/structure_function_coupling',subject_IDs(i),subject_IDs(i)+'_voxel_based_threshold_'+string(threshold)+'_percent_structure_function_coupling_concatenated.mat'),'coupling');
end
