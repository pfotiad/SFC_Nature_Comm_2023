function [myelin_content] = average_myelin_data_schaefer_atlas(subject)
% The purpose of this script is to compute the average intracortical myelin 
% content on a provided Human Connectome Project (HCP) subject using the
% Schaefer atlas-defined brain regions.
% 
% Requirements:
% 1. The publicly available 100 unrelated HCP dataset distribution: 
% https://db.humanconnectome.org/
% 2. The downloaded Schaefer parcellations: 
% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations
% 3. To have run the Bash shell script T1w_T2w_ratio_native_to_MNI.sh
% script that registers the subject's myelin map from native to MNI152 
% standardized space

% Input:
% subject -> Subject ID of interest

% Output: 
% 1. myelin_content -> matrix of size <400x3> containing the brain region 
% ID in the first column, the mean estimate of intracortical myelin of each
% brain region in the second column, and the standard deviation of
% intracortical myelin of each brain region in the third column
% 2. Saved .mat file containing that variable, saved at each subject's directory
%
% Author: Panagiotis Fotiadis, 2021-2023
%
% If parts of the following code are used, please cite the following paper:
%
% Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD, 
% Shinohara RT, and Bassett DS “Myelination and excitation-inhibition 
% balance synergistically shape structure-function coupling across the 
% human cortex,” (2023) Nature Communications.


% Load the atlas
atlas_file = MRIread('/cbica/home/panosf/software/Schaefer_atlases_MNI/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz');
    
% Convert the 3D atlas map into a vector
atlas_tmp = permute(atlas_file.vol,[2 1 3]);
atlas = atlas_tmp(:);

% Load the subject's myelin map in MNI152 standardized space
hcp_myelin_map = MRIread(fullfile('/cbica/home/panosf/HCP_Dataset/',subject,'T1wDividedByT2w_ribbon_MNI152.nii.gz'));

% Convert the 3D myelin map into a vector
hcp_myelin_tmp = permute(hcp_myelin_map.vol,[2 1 3]);
hcp_myelin = hcp_myelin_tmp(:);

% Append the atlas and myelin maps column-wise
myelin = [atlas hcp_myelin];
    
% Exclude all rows corresponding to an atlas region of interest [ROI] = 0
zero_indices = [];
zero_indices = find(myelin(:,1) == 0);
myelin(zero_indices,:) = [];
    
% number of brain regions
num_ROIs = 400;
    
% number of standard deviations (see below)
num_stdev = 1;
    
% for each ROI, find all voxels that correspond to it and average their 
% T1w/T2w ratio signal intensity while only keeping values within a
% certain range
for k=1:num_ROIs
        
    % initialize indices
    hcp_idx = [];
        
    % indices within range that contain positive values for both scans
    hcp_idx = find((myelin(:,1) == k) & (myelin(:,2) > 0) & (myelin(:,2) <= (mean(nonzeros(myelin(:,2)),'omitnan') + num_stdev * std(nonzeros(myelin(:,2)),'omitnan'))) & (myelin(:,2) >= (mean(nonzeros(myelin(:,2)),'omitnan') - num_stdev * std(nonzeros(myelin(:,2)),'omitnan'))));
        
    % calculate the average T1w/T2w signal intensity values
    myelin_content(k,1) = k; % brain region ID
    myelin_content(k,2) = mean(myelin(hcp_idx,2)); % Mean myelin for ROI ID
    myelin_content(k,3) = std(myelin(hcp_idx,2)); % STDev myelin for ROI ID
end
    
% Save the resulting myelin values in the subject's directory
save(fullfile('/cbica/home/panosf/HCP_Dataset/',subject,'Schaefer2018_400Parcels_7Networks_myelin_content.mat'),'myelin_content');

