function [coupling] = structure_function_coupling(subject)
% The purpose of this script is to compute the structure-function coupling
% based on the correlational approach, on each HCP subject's brain region.
% 
% Requirements:
% 1. Download the gifti and cifti-matlab toolboxes. Instructions can be found in:
% https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ
% 2. Structural and functional connectivity matrices have been generated
% for each subject, using the scripts structural_connectome_generation.sh
% and functional_connectome_generation.sh, respectively.
%
% Input:
% subject -> The subject ID of interest
%
% Outputs:
% coupling -> Nested variable containing structure-function coupling 
% indices corresponding to the subject and the atlas of interest. Each 
% coupling variable will be of size: number of brain regions x 1. A .mat 
% file is also saved on the defined output directory. 
%
% Author: Panagiotis Fotiadis, 2021-2023
%
% If parts of the following code are used, please cite the following paper:
%
% Fotiadis P, Cieslak M, He X, Caciagli L, Ouellet M, Satterthwaite TD, 
% Shinohara RT, and Bassett DS “Myelination and excitation-inhibition 
% balance synergistically shape structure-function coupling across the 
% human cortex,” (2023) Nature Communications.


% Add the appropriate paths pointing to the downloaded software
addpath('/cbica/home/panosf/software/gifti');
addpath('/cbica/home/panosf/software/cifti-matlab');

% Define the output directory
output_dir = '/cbica/home/panosf/HCP_Dataset';

%% Compute the structure-function coupling for each brain region

% Define the atlases of interest
atlases = ["Schaefer2018_400Parcels_7Networks_order", "Glasser_cortical"];

% The total number of time windows used in the analysis
num_windows = 20;

% Perform the structure-function coupling calculations
for j=1:numel(atlases) % Go through each atlas
    
    % Load the subject's structural connectivity matrix
    str_conn_matrix{j} = readmatrix(fullfile(output_dir,subject,'structural_connectivity',atlases(j)+'_connectome.csv'));
    
    % Load the subject's static functional connectivity matrices
    static_func_conn_matrix{j}.raw = ciftiopen(fullfile(output_dir,subject,'functional_connectivity',atlases(j)+'_connectome.pconn.nii'),'/cbica/software/external/connectome_workbench/1.4.2/bin/wb_command'); % Pearson's correlation values
    static_func_conn_matrix{j}.raw = static_func_conn_matrix{j}.raw.cdata;

    % Load the subject's temporally-contiguous functional connectivity matrices for
    % all time windows and concatenate them into one 3D variable
    for win = 1:num_windows
        temporal_func_conn_matrix{j}.raw_tmp{win} = ciftiopen(fullfile(output_dir,subject,'functional_connectivity/temporal_analyses',atlases(j)+'_connectome_'+win+'.pconn.nii'),'/cbica/software/external/connectome_workbench/1.4.2/bin/wb_command');
        temporal_func_conn_matrix{j}.raw(:,:,win) = temporal_func_conn_matrix{j}.raw_tmp{win}.cdata;
    end
    
    % Define some housekeeping variable names
    coupling.atlas{j}.atlas = atlases(j);
    
    % Calculate structure-function coupling for each brain region
    for m=1:size(str_conn_matrix{j},1)

        % Calculate structure-function coupling given the static functional
        % connectivity data; select elements where both structural and
        % functional connectivities are non-zero
        
        % Define some tmp variables
        X = str_conn_matrix{j}(:,m);
        Y = static_func_conn_matrix{j}.raw(:,m);
        
        IdxBothNonZero = X ~= 0 & Y ~= 0;
        IdxBothNonZero(m) = 0; % also set the (m,m) index into 0, for it not to skew the data
        
        % Calculate the coupling
        if (isempty(X(IdxBothNonZero))) || (isempty(Y(IdxBothNonZero)))
            coupling.atlas{j}.raw.static_correlational_approach_non_zero_data_pearson_rho(m,1) = NaN;
            coupling.atlas{j}.raw.static_correlational_approach_non_zero_data_pearson_pvalue(m,1) = NaN;
        else
            % Parametrically
            [coupling.atlas{j}.raw.static_correlational_approach_non_zero_data_pearson_rho(m,1), coupling.atlas{j}.raw.static_correlational_approach_non_zero_data_pearson_pvalue(m,1)] = corr(X(IdxBothNonZero),Y(IdxBothNonZero),'rows','pairwise');
        end

        % Clear tmp variables
        clear X Y IdxBothNonZero;
        
        % Calculate coupling given the temporal functional connectivity data
        
        % Calculate structure-function coupling given the temporal functional 
        % connectivity data - select elements where both structural and 
        % functional connectivities are non-zero
        for n=1:num_windows
            
            % Define some tmp variables
            X = str_conn_matrix{j}(:,m);
            Y = temporal_func_conn_matrix{j}.raw(:,m,n);

            IdxBothNonZero = X ~= 0 & Y ~= 0;
            IdxBothNonZero(m) = 0; % also set the (m,m) index into 0, for it not to skew the data
            
            % Calculate the coupling
            if (isempty(X(IdxBothNonZero))) || (isempty(Y(IdxBothNonZero)))
                coupling.atlas{j}.raw.temporal_correlational_approach_non_zero_data_pearson_rho(m,n) = NaN;
                coupling.atlas{j}.raw.temporal_correlational_approach_non_zero_data_pearson_pvalue(m,n) = NaN;
            else
                % Parametrically
                [coupling.atlas{j}.raw.temporal_correlational_approach_non_zero_data_pearson_rho(m,n), coupling.atlas{j}.raw.temporal_correlational_approach_non_zero_data_pearson_pvalue(m,n)] = corr(X(IdxBothNonZero),Y(IdxBothNonZero),'rows','pairwise');
            end

            % Clear the tmp variables
            clear X Y IdxBothNonZero;
        end
    end
end

% Also save the output variable into a .mat file
save(fullfile(output_dir,subject,'structure_function_coupling.mat'),'coupling');

