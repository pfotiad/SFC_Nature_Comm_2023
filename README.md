# SFC_Nature_Comm_2023
**Author:** Panagiotis Fotiadis, 2021-2023

This repository contains code that was created to conduct the main analyses described in our paper:

**"Myelination and excitation-inhibition balance synergistically shape structure-function coupling across the human cortex," 
currently accepted for publication in *Nature Communications* (2023).** Link to be soon provided.

Bash Shell Scripting, *MATLAB* (version R2021a: The MathWorks, Inc.), and *Python* (version 3.7) were used to write the code discussed below. 

We separate the code into two sections (folders): 
one corresponding to our atlas-based analyses and one corresponding to our voxel-based analyses.

## **A. Atlas-based analyses:**

  1. <ins>Structural connectivity</ins>

     To generate the structural connectivity matrix corresponding to the subject of interest, run the Bash shell script: `structural_connectome_generation.sh`.
     Prerequisites and further instructions can be found at the top of the script.

  2. <ins>Functional connectivity</ins>

     To generate the static and temporally-contiguous functional connectivity matrices corresponding to the subject of interest, run the Bash shell script: `functional_connectome_generation.sh`.
     Prerequisites and further instructions can be found at the top of the script.

  3. <ins>Structure-Function Coupling</ins>

     To compute the overall structure-function coupling and temporal structure-function coupling variance corresponding to the subject of interest, run the MATLAB script: `structure_function_coupling.m`.
     Prerequisites and further instructions can be found at the top of the script.

  4. <ins>Intracortical Myelination</ins>

     a. To estimate the intracortical myelin content of the subjects of interest using the HCP multi-modal (aka Glasser) atlas, run the Bash shell script: average_myelin_data_glasser_atlas.sh.
     
     b. To estimate the intracortical myelin content of the subjects of interest using the Schaefer atlas, first run the Bash shell script: T1w_T2w_ratio_native_to_MNI.sh, and then run the MATLAB script: average_myelin_data_schaefer_atlas.m.
     Prerequisites and further instructions can be found at the top of each script.

  5. <ins>Excitation-Inhibition Ratio (i.e., Hurst Exponent)</ins>

     To compute the Hurst exponent (and by proxy assess the excitation-inhibition ratio) of the subject of interest, run the MATLAB script: `compute_hurst_exponent.m`.
     Prerequisites and further instructions can be found at the top of the script.


## **B.  Voxel-based analyses:**

  1. <ins>Structural connectivity</ins>

     To generate the thresholded voxel-based structural connectivity matrix corresponding to the subject of interest:
     
          a. Run the Python script: create_voxel_conn_matrix.py to generate un-thresholded connectivity matrices (using multiple definitions of edge weights),
             for the subject of interest (in parts; .h5 format),
     
          b. Run the MATLAB script: convert_h5_voxel_based_structural_connectivity_to_mat.m to convert the .h5 files output from (a) into .mat files,
     
          c. Run the MATLAB script: extract_omitted_voxel_IDs.m in order to identify any cortical voxels that might have been skipped over in the calculation of the structural connectivity matrices,so that we can make sure that all subsequent analyses contain and refer to the exact same number of cortical voxels (the output file from this script will be used later on), and
     
          d. Run the MATLAB script: threshold_voxel_based_structural_connectivity_matrices.m which will generate thresholded voxel-based structural connectivity matrices for each subject.
     Prerequisites and further instructions can be found at the top of each script.

  2. <ins>Functional connectivity</ins>

     To generate the voxel-based static and temporally-contiguous resting-state functional connectivity matrices corresponding to the subject of interest:
     
          a. Process the subjects' resting-state functional MRI data using the CONN toolbox (https://web.conn-toolbox.org/resources/installation), as described in the manuscript,
     
          b. Run the MATLAB script: prepare_functional_data.sh to generate some files that will be necessary later on when we generate the voxel-based functional connectivity matrices,
     
          c. Run the MATLAB script: extract_voxel_based_BOLD_signal.m to generate the voxel-based BOLD signal time series of the subject of interest,
     
          d. Run the MATLAB script: compute_voxel_based_functional_connectivity.m to generate the voxel-based functional connectivity matrix of the subject of interest, and
     
          e. Run the MATLAB scripts: compute_voxel_based_temporal_functional_connectivity.m  and concatenate_voxel_based_temporal_functional_connectivity.m to generate the temporally-contiguous voxel-based functional connectivity matrices for the subject of interest.
     Prerequisites and further instructions can be found at the top of each script.

  3. <ins>Structure-Function Coupling</ins>

     To compute the voxel-based structure-function coupling and temporal structure-function coupling variance corresponding to the subject of interest, run the MATLAB scripts: `extract_voxel_based_structure_function_coupling.m` and `concatenate_voxel_based_structure_function_coupling_files.m`.
     Prerequisites and further instructions can be found at the top of each script.

  4. <ins>Intracortical Myelination</ins>

     To estimate the intracortical myelin content of the subject(s) of interest:
     
          a. Process the T1- and T2-weighted sequences of the Penn subjects using the HCP pre-processing pipelines: (https://github.com/Washington-University/HCPpipelines/wiki/Installation-and-Usage-Instructions#running-the-hcp-pipelines-on-example-data),as described in the manuscript; these scripts will generate the T1-weighted/T2-weighted 'myelin maps' for each subject as a proxy of their intracortical myelin content,
       
          b. Run the Bash shell scripts: HCP_orig_to_b0.sh and register_T1w_T2w_files_to_b0_HCP.sh, to register the T1-weighted/T2-weighted map from (a) into the subject's native diffusion (b0) space, and
     
          c. Run the MATLAB script: extract_voxel_based_myelin_content.m to generate the intracortical myelin estimates for each subject, at the voxel-level.

  5. <ins>Excitation-Inhibition Ratio (i.e., Hurst Exponent)</ins>

     To compute the voxel-based Hurst exponent (and, by proxy, assess the excitation-inhibition ratio) of the subject of interest, run the MATLAB scripts: `compute_voxel_based_hurst_exponent.m` and `concatenate_voxel_based_hurst_exponent.m`.
     Prerequisites and further instructions can be found at the top of each script.
