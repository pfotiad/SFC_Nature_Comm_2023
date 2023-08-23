% The purpose of this script is to estimate the intracortical myelin content for each
% voxel of the Penn subjects, using the T1w/T2w method

% Requirements:
% 1. The scripts HCP_orig_to_b0.sh and register_T1w_T2w_files_to_b0_HCP.sh have been run.
%
% Output: 
% The voxel-based intracortical myelin estimates (within 1 standard deviation) for each subject's
% cortical voxel, saved as a .mat file.
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

% Clear things up
clear; clc;

% Define the subject IDs and the dataset they correspond to (this should be
% adjusted to your dataset)
subject_IDs = ["sub-EHJ", "sub-HJK", "sub-JFG", "sub-LT", "sub-NWT", "sub-VW", "sub-YH", "sub-s10", "sub-s11"];
dataset = [2 2 2 2 2 2 2 1 1];

% Loop through the subjects
for i=1:numel(subject_IDs)
    
    % Define the subject's cortical voxel-based atlas
    atlas_file{i} = MRIread(char(fullfile('/cbica/home/panosf/Q7/CONN_functional_results', subject_IDs(i), subject_IDs(i)+'_ses-1_acq-Q7_space-T1w_desc-preproc_space-T1w_desc-schaefer100x7_atlas_cortex_voxel_based.nii.gz')));

    % Convert the 3D matrix into a vector
    atlas_tmp{i} = permute(atlas_file{i}.vol,[2 1 3]);
    atlas{i} = atlas_tmp{i}(:);
    
    % Load the myelin maps and convert them into columns
    
    % HCP-derived maps
    hcp_myelin_map{i} = MRIread(char(fullfile('/cbica/home/panosf/Q7/structure_function_coupling/HCP_Processing_Pipeline',subject_IDs(i),'T1w/T1wDividedByT2w_ribbon_b0.nii.gz')));
    
    hcp_myelin_tmp{i} = permute(hcp_myelin_map{i}.vol,[2 1 3]);
    hcp_myelin{i} = hcp_myelin_tmp{i}(:);

    % Append the atlas and myelin maps column-wise
    myelin{i} = [atlas{i} hcp_myelin{i}];
    
    % Exclude all rows corresponding to an atlas region of interest (ROI) = 0 
    zero_indices{i} = [];
    zero_indices{i} = find(myelin{i}(:,1) == 0);
    myelin{i}(zero_indices{i},:) = [];
    
    % Sort the myelin{i} matrix based on the voxels ID on column 1
    myelin{i} = sortrows(myelin{i});
    
    % Load the voxel IDs to omit so that we can have a 1-1 matching with
    % the structural and functional connectivity voxel IDs.
    omitted_voxels{i} = importdata(fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset(i)), 'Nifti', subject_IDs(i), 'mittens_results', string(subject_IDs(i))+'_structural_connectivity_omitted_voxel_IDs.txt'));
    
    % Use that list to create a new cleaned up version of the myelin data for
    % each subject that does not include these omitted voxel ID rows.
    myelin{i}(omitted_voxels{i},:) = [];
      
    % Keep only the myelin estimates within a pre-defined standard deviation,
    % as suggested by Glasser et al 2011.
    num_stdev = 1;

    % find indices with myelin content out of range
    hcp_idx{i} = find((myelin{i}(:,2) == 0) | (myelin{i}(:,2) > (mean(nonzeros(myelin{i}(:,2)),'omitnan') + num_stdev * std(nonzeros(myelin{i}(:,2)),'omitnan'))) | (myelin{i}(:,2) < (mean(nonzeros(myelin{i}(:,2)),'omitnan') - num_stdev * std(nonzeros(myelin{i}(:,2)),'omitnan'))));

    myelin{i}(hcp_idx{i},2) = NaN;
end

% Save the results
for i=1:numel(subject_IDs)
    
    % Temporarilly save the output variable
    clear myelin_output;
    myelin_output = myelin{i};
    
    % Save the resulting myelin values for each subject
    save(fullfile('/cbica/home/panosf/Q7/structure_function_coupling',subject_IDs(i),'myelin_results/voxel_based_myelin_content.mat'),'myelin_output');
    
    % Clear the output variable
    clear myelin_output;
end
    