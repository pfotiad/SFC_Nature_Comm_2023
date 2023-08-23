function [] = compute_voxel_based_hurst_exponent(subject,part)

% The purpose of this script is to compute the Hurst exponent of the Penn
% sample's diffusion voxel-based functional timeseries, using the nonfractal toolbox:
% https://github.com/wonsang/nonfractal as descrived in the paper:
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7402681/pdf/elife-55684.pdf

% Requirements:
% 1. The subjects' resting-state functional MRI data have been preprocessed 
% using the CONN toolbox (https://web.conn-toolbox.org/resources/installation), 
% as described in the manuscript
% 2. Download the gifti and cifti-matlab toolboxes. Instructions can be found
% in: https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ
% 3. Download the wmtsa-matlab-0.2.6 and nonfractal toolboxes. Instructions
% can be found in: https://github.com/wonsang/nonfractal
% 4. The script extract_voxel_based_BOLD_signal.m has been run.
%
% Input:
% 1. subject -> Subject ID
% 2. part    -> integer from 1 to 20 (run it in 20 parts due to memory limitations)
%
% Output:
% 20 csv files containing the Hurst exponent corresponding to the designated 
% subject's cortical voxels.
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

% Arguments for the bfn_mfin_ml.m command later on, as defined in the paper
% cited above
lb = [-0.5,0];
ub = [1.5,10];

% Initialize variables
BOLD_timeseries = [];
H = [];

% Load the voxel-based BOLD signal time series
BOLD_timeseries = load(fullfile('/cbica/home/panosf/Q7/CONN_functional_results',subject,'BOLD_processing',string(subject)+'_BOLD_timeseries_matching_structural_connectivity_voxel_IDs.mat'));
BOLD_timeseries = BOLD_timeseries.BOLD;
BOLD_timeseries(:,1) = []; % Remove the first column as it displays the voxel ID 

% Trim the BOLD timeseries in parts (because otherwise the algorithm won't
% work)
start_id = (part - 1) * (floor(size(BOLD_timeseries,1)/20)) + 1;

if (part < 20)
    end_id = part * (floor(size(BOLD_timeseries,1)/20));
else
    end_id = size(BOLD_timeseries,1);
end

% define the trimmed BOLD signal time series
trimmed_BOLD_timeseries = BOLD_timeseries(start_id:end_id,:);

% Transpose the matrix
trimmed_BOLD_timeseries = trimmed_BOLD_timeseries';

% Run the nonfractal toolbox command "bfn_mfin_ml" to estimate the Hurst exponent, nonfractal connectivity, and fractal connectivity
% in multivariate time series with long memory using the maximum likelihood method
[H, ~, ~] = bfn_mfin_ml(trimmed_BOLD_timeseries, 'filter', 'Haar', 'lb', lb, 'ub', ub);

% Output results
results = [];
results = H';

% Save the results as a csv file
writematrix(results,fullfile('/cbica/home/panosf/Q7/structure_function_coupling',subject,string(subject)+'_voxel_based_ei_hurst_results_part_'+string(part)+'.csv'));
