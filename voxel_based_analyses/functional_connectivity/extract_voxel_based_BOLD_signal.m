function [BOLD] = extract_voxel_based_BOLD_signal(subject,dataset,atlas_path)
% The purpose of this script is to extract the BOLD timeseries from each TR
% as designated by the prepare_functional_data.sh script, for each of the 
% subject's cortical voxels.

% Requirements:
% 1. Download the gifti and cifti-matlab toolboxes. Instructions can be found in:
% https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ
% 2. The script prepare_functional_data.sh has been run.

% Input:
% 1. subject -> Subject ID
% 2. dataset -> The dataset ID of interest; this was mainly for my benefit 
%               since I had 2 different datasets that I was analyzing in different directories; 
%               feel free to remove this input after adjusting the paths 
%				below to fit your data.
% 2. atlas_path -> Full path to the subject's cortical voxel atlas.

% Output: 
% 1. BOLD -> voxel-based BOLD_signal with size <number of cortical voxels> x 2400 (number of TRs) 
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

% Add the appropriate paths
addpath('/cbica/home/panosf/software/gifti');
addpath('/cbica/home/panosf/software/cifti-matlab');

% Load the cortical atlas 
atlas_file = MRIread(atlas_path);

% Convert the 3D atlas map into a vector
atlas_tmp = permute(atlas_file.vol,[2 1 3]);
atlas = atlas_tmp(:);

% Define a new variable that will contain the atlas and BOLD TR info
BOLD(:,1) = atlas;

% Load all TRs and convert them into columns
for i=1:2400
    vol = i-1;
    BOLD_TR_file{i} = MRIread(char(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing/BOLD_TR_volumes',string(subject)+'_BOLD_MNI_to_b0_TR'+num2str(vol, '%04.f')+'.nii.gz')));
    BOLD_TR_tmp{i} = permute(BOLD_TR_file{i}.vol,[2 1 3]);
    BOLD_TR{i} = BOLD_TR_tmp{i}(:);
    BOLD(:,i+1) = BOLD_TR{i};
end

% Exclude all rows corresponding to a voxel with a zero value 
zero_indices = find(BOLD(:,1) == 0);
BOLD(zero_indices,:) = [];

% Sort the BOLD array based on the voxels ID on column 1
BOLD = sortrows(BOLD);

% Save the full BOLD signal time-series matrix to have a record of it
save(char(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_BOLD_timeseries.mat')),'BOLD');

% Import the subject's omitted voxels as defined by the script
% extract_omitted_voxel_IDs.m defined in the voxel-based structural connectivity
% analyses
omitted_voxels = importdata(fullfile('/cbica/home/panosf/Q7/Q7_Dataset_'+string(dataset), 'Nifti', subject, 'mittens_results', string(subject)+'_structural_connectivity_omitted_voxel_IDs.txt'));

% Use that list to create a new cleaned up version of the subject's BOLD_timeseries 
% that does not include these voxel ID rows.
BOLD(omitted_voxels,:) = [];

% Save the cleaned BOLD timeseries
save(char(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_BOLD_timeseries_matching_structural_connectivity_voxel_IDs.mat')),'BOLD');