% The purpose of this script is to compute the Hurst exponent of all of the
% 100 unrelated HCP subjects' BOLD timeseries—-as a proxy of their 
% excitation-inhibition ratio--using the nonfractal toolbox:  
% https://github.com/wonsang/nonfractal 
% described in the paper:
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7402681/pdf/elife-55684.pdf
% 
% Requirements:
% 1. The publicly available 100 unrelated HCP dataset distribution: 
% https://db.humanconnectome.org/
% 2. Download the gifti and cifti-matlab toolboxes. Instructions can be found
% in: https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ
% 3. Download the wmtsa-matlab-0.2.6 and nonfractal toolboxes. Instructions
% can be found in: https://github.com/wonsang/nonfractal
% 4. Functional connectivity matrices have been generated for each subject
% using the script functional_connectome_generation.sh.
%
% Output: 
% The resulting csv file will contain the Hurst exponent
% corresponding to each brain region (depending on the atlas you define 
% below). The size of the resulting file will be <number of ROIs> x 1. One
% csv file per subject will be generated.
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
addpath(genpath('/cbica/home/panosf/software/wmtsa-matlab-0.2.6'));
addpath(genpath('/cbica/home/panosf/software/nonfractal'));

% Arguments for the bfn_mfin_ml.m command used later on, as defined in the
% paper cited above
lb = [-0.5,0];
ub = [1.5,10];

% Initialize variables
hcp_BOLD_timeseries = {};
hcp_H = {};

disp('Processing the HCP subjects:')

% Load the HCP subject IDs: The text file below contains a subject ID per
% row. This text file is provided by the HCP download. 
hcp_subject_IDs = importdata('/cbica/home/panosf/HCP_Dataset/100_unrelated_subjects_list.txt');

% Define the atlases of interest
hcp_atlases = ["Schaefer2018_400Parcels_7Networks_order", "Glasser_cortical"];

for i=1:numel(hcp_subject_IDs) % loop through the subjects
    
    for a=1:numel(hcp_atlases) % loop through the atlases
        
        % Load the BOLD timeseries: After the following two lines are ran
        % the size of the BOLD timeseries will be <number of TRs> x <number
        % of ROIs>
        % TR is repetition time and ROIs is regions of interest
        hcp_BOLD_timeseries{i}{a} = ciftiopen(fullfile('/cbica/home/panosf/HCP_Dataset',num2str(hcp_subject_IDs(i)),'functional_connectivity',hcp_atlases(a)+'_timeseries.ptseries.nii'),'/cbica/software/external/connectome_workbench/1.4.2/bin/wb_command');
        
        hcp_BOLD_timeseries{i}{a} = double(hcp_BOLD_timeseries{i}{a}.cdata');
        
        % Run the nonfractal toolbox command "bfn_mfin_ml" to estimate the Hurst exponent, nonfractal connectivity, and fractal connectivity
        % in multivariate time series with long memory using the maximum likelihood method
        [hcp_H{i}{a}, ~, ~] = bfn_mfin_ml(hcp_BOLD_timeseries{i}{a}, 'filter', 'Haar', 'lb', lb, 'ub', ub);
        
        % Concatenate the HCP results -> size(hcp_results) = <number of brain regions> x 1
        hcp_results = [];
        hcp_results = hcp_H{i}{a}';
        
        % Save the results as a csv file per subject
        writematrix(hcp_results,fullfile('/cbica/home/panosf/HCP_Dataset',num2str(hcp_subject_IDs(i)),string(hcp_subject_IDs(i))+'_'+hcp_atlases(a)+'_ei_hurst_results.csv'));        
    end
end

disp('Done!')
