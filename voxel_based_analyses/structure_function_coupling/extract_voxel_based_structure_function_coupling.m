function [coupling] = extract_voxel_based_structure_function_coupling(subject, dataset, part)
% The goal of this script is to calculate the structure-function coupling
% across each voxel. It assumes that the voxel-based (thresholded) structural 
% and functional connectivity matrices for that subject have been generated. 
% Currently, I should have an n-number of .mat files for each type of 
% connectivity matrix, each containing information for 1000 voxels, with a 
% 1-1 correspondence between them (i.e., structural to functional matrix 
% row correpsondence).
% 
% Requirements:
% 1. All scripts under the voxel-based structural and functional connectivity 
%    sections of the README.md document have been run.

% Inputs: 
% 1. subject -> subject ID 
% 2. dataset -> The dataset ID of interest; this was mainly for my benefit 
%               since I had 2 different datasets that I was analyzing in different directories; 
%               feel free to remove this input after adjusting the paths 
%				below to fit your data.
% 3. part    -> part ID ranging from 1:ceil(size(BOLD,1)/1000)

% Output: 
% 1. coupling -> Nested variable containing structure-function coupling 
% indices corresponding to the subject. Each coupling variable will be of 
% size: number of cortical voxels x 1.
% 2. The above variable saved as a .mat file; structure-function coupling 
% is saved in parts once again due to the sheer amount of data and
% computational power required to run these analyses.
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


% Define the sliding window-related parameters
num_windows = 20;

% Initialize variables
clear coupling

% Load the subject's thresholded structural connectivity matrix, with the edge weights 
% defined using the variable 'mean_scores' 
structural_data_mean_scores = load(fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset), 'Nifti', subject, 'mittens_results', string(subject)+'_edge_weights_mean_scores_threshold_70_percent_n'+string(part)+'.mat'));

structural_data_mean_scores = structural_data_mean_scores.data;

% Load the subject's functional connectivity matrix
% Static data
functional_data = load(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_functional_connectivity_matrix_part_'+string(part)+'.mat'));
functional_data = functional_data.FC;

% Temporal data
temporal_functional_data = load(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_temporal_functional_connectivity_matrix_part_'+string(part)+'.mat'));
temporal_functional_data = temporal_functional_data.temporal_FC;

% Compute the structure-function coupling
for row = 1:size(functional_data,1)
    
    % STATIC CASE

    % Calculate structure-function coupling given the static functional
    % connectivity data - select elements where both structural
    % and functional connectivity are non-zero
    
    % Start with a clean slate
    clear X_mean_scores Y IdxBothNonZero_mean_scores;
    
    X_mean_scores = structural_data_mean_scores(row,:)';
    
    Y = functional_data(row,:)';
    
    IdxBothNonZero_mean_scores = X_mean_scores ~= 0 & Y ~= 0;
    IdxBothNonZero_mean_scores(row) = 0; % also set the (row,row) index into 0, for it not to skew the data
    
    % Calculate the coupling based on the mean_scores values
    if (isempty(X_mean_scores(IdxBothNonZero_mean_scores))) || (isempty(Y(IdxBothNonZero_mean_scores)))
        coupling.mean_scores.static_correlational_approach_non_zero_data_spearman_rho(row,1) = NaN;
        coupling.mean_scores.static_correlational_approach_non_zero_data_spearman_pvalue(row,1) = NaN;
    else
        % Non-parametrically
        [coupling.mean_scores.static_correlational_approach_non_zero_data_spearman_rho(row,1), coupling.mean_scores.static_correlational_approach_non_zero_data_spearman_pvalue(row,1)] = corr(X_mean_scores(IdxBothNonZero_mean_scores),Y(IdxBothNonZero_mean_scores),'Type','Spearman','rows','pairwise');
    end
    
    % TEMPORAL CASE: Calculate coupling given the temporal functional
    % connectivity data

    % Calculate structure-function coupling given the temporal functional
    % connectivity data - select elements where both structural
    % and functional connectivity are non-zero
    for win = 1:num_windows
        
        % Start with a clean slate
        clear X_mean_scores Y IdxBothNonZero_mean_scores;
        
        % Re-define the parameters
        X_mean_scores = structural_data_mean_scores(row,:)';
        
        Y = temporal_functional_data(row,:,win)';
        
        IdxBothNonZero_mean_scores = X_mean_scores ~= 0 & Y ~= 0;
        IdxBothNonZero_mean_scores(row) = 0; % also set the (row,row) index into 0, for it not to skew the data
        
        % Calculate the coupling based on the mean_scores values
        if (isempty(X_mean_scores(IdxBothNonZero_mean_scores))) || (isempty(Y(IdxBothNonZero_mean_scores)))
            coupling.mean_scores.temporal_correlational_approach_non_zero_data_spearman_rho(row,win) = NaN;
            coupling.mean_scores.temporal_correlational_approach_non_zero_data_spearman_pvalue(row,win) = NaN;
        else
            % Non-parametrically
            [coupling.mean_scores.temporal_correlational_approach_non_zero_data_spearman_rho(row,win), coupling.mean_scores.temporal_correlational_approach_non_zero_data_spearman_pvalue(row,win)] = corr(X_mean_scores(IdxBothNonZero_mean_scores),Y(IdxBothNonZero_mean_scores),'Type','Spearman','rows','pairwise');
        end
    end
end

% Save the variable
save(fullfile('/cbica/home/panosf/Q7/structure_function_coupling',subject,string(subject)+'_voxel_based_threshold_70_percent_structure_function_coupling_part_'+string(part)+'.mat'),'coupling');
