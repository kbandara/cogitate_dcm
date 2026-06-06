%% ALL SOURCE ANALYSIS UP TO DCM ANALYSIS 
% Kav Bandara, University of Melbourne, 2025

%{

This script runs source reconstruction and relevant contrasts; including
model inversion, image extraction, GLMs, extracting peak mni coords

%}

clear all; close all

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12')

spm('defaults', 'eeg');

%%%%%%%%%%           SECTION 0 - LOAD FILES & DEFINE VARIABLES          %%%%%%%%%%%%

%% Settings and filenames for pre-processing

base_filepath  = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/'; 
cd(base_filepath)

%Load participant names
P_Names_struct = load('p_names.mat');  
P_Names_table =  P_Names_struct.uniqueToFirstList; %{'CA107'}; 
P_Names_str = P_Names_table;
p_names     = cellstr(P_Names_str);  % 1x52 cell array of char

% for reference: epochTimeWindow = [-100 2500];
epochSourceTimeWindow = [0 2000]; %time window for source reconstruction

filepath = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessed_files/';
cd(filepath)

image_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'onset_contrasts');
cd(image_dir) 

%% Run source inversion process
% Steps include: head template, coregistrations, define and invert forward
% model, specify time/frequency window and whether to use taper, create
% images of the results

for si = 1:length(p_names) 

    filename = fullfile(filepath, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' p_names{si} '_ses-1_task-dur_run-01_meg.mat']); 
    load(filename);    
    
    %load anatomical scans and fiducials 
    anat_mri = fullfile(base_filepath,['sub-' p_names{si}], 'ses-1', 'anat', ['sub-' p_names{si} '_ses-1_T1w.nii']);     
    anat_fid = fullfile(base_filepath,['sub-' p_names{si}], 'ses-1', 'anat', ['sub-' p_names{si} '_ses-1_T1w.json']);
    
    % Check if subject has individual MRI
    has_mri = isfile(anat_mri);

    spm_jobman('initcfg');

    %Source space modelling (using templates), Coregister, Forward Model
    matlabbatch{1}.spm.meeg.source.headmodel.D = {filename};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    if has_mri
        % Create mesh using subject's sMRI
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {[anat_mri ',1']};
    else
        % Use template mesh (canonical head)
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
    end
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'LPA';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'lpa'; % enter 1 × 3 vector of 3d coordinates for each fiducial
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'Nasion';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'nas'; % enter 1 × 3 vector of 3d coordinates for each fiducial
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa'; % enter 1 × 3 vector of 3d coordinates for each fiducial
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0; % 0 = 'no, no is advised for MEG with sMRI; 1 is advised for EEG
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';

    %Model inversion and specify time (and frequency) window, using Multiple Sparse Priors algorithm 
    matlabbatch{2}.spm.meeg.source.invert.D = {filename};
    matlabbatch{2}.spm.meeg.source.invert.val = 1; %save inversion to new index
    matlabbatch{2}.spm.meeg.source.invert.whatconditions.all = 1; %use all conditions
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS'; %Multiple sparse priors (greedy search)
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.woi = epochSourceTimeWindow; %time window of interest 
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.foi = [0 256]; %frequency window of interest 
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.hanning = 0; %1 = Hanning taper at start and end of trial (0 = no taper, 1 = yes)
    matlabbatch{2}.spm.meeg.source.invert.modality = {'All'}; %note that all other channels have been cropped out already

    spm_jobman('run',matlabbatch);
        
    clear D 
    clear spm_jobman
    clear matlabbatch
 
end

%% Create images

onset_timewindow = [0 600];
offset_timewindow = [];

for onset = 1:length(p_names) 

    filename = fullfile(filepath, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' p_names{onset} '_ses-1_task-dur_run-01_meg.mat']); 
        
    load(filename);
    spm_jobman('initcfg');
    %Display results within time (and frequency) window and create images
    matlabbatch{1}.spm.meeg.source.results.D = {filename};
    matlabbatch{1}.spm.meeg.source.results.val = 1; %display inversion results from index
    matlabbatch{1}.spm.meeg.source.results.woi = onset_timewindow; % time of interest
    matlabbatch{1}.spm.meeg.source.results.foi = [0 256]; % frequency window specify
    matlabbatch{1}.spm.meeg.source.results.ctype = 'evoked'; % ORIGINALLY 'trials' %'evoked' 'induced' or single 'trials'
    matlabbatch{1}.spm.meeg.source.results.space = 1; % 1 = MNI or Native
    matlabbatch{1}.spm.meeg.source.results.format = 'image';
    matlabbatch{1}.spm.meeg.source.results.smoothing = 12; %Smoothing mm^3
        
    spm_jobman('run',matlabbatch);
        
    clear D 
    clear spm_jobman
    clear matlabbatch
 
end

%% MOVE FILES TO ONSET CONTRASTS FOLDER 

% setup filepaths
base_filepath  = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/'; 

input_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files');
image_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'onset_contrasts');

% Create the output directory if it doesn't exist
if ~exist(image_dir, 'dir')
    mkdir(image_dir);
end

cd(image_dir)

for move = 1:length(p_names)

    for i = 1:16 %for each condition/image generated above 
    
    og_file = fullfile(input_dir, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' p_names{move} '_ses-1_task-dur_run-01_meg_1_t0_600_f0_256_' num2str(i) '.nii']);
    target_file = fullfile(image_dir, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' p_names{move} '_ses-1_task-dur_run-01_meg_1_t0_600_f0_256_' num2str(i) '.nii']); 
       
    %move the files
    movefile(og_file, target_file, 'f'); % 'f' to force overwrite if target exists

    clear og_file target_file; 

    end
end

%% RENAME SOURCE RECONSTRUCTED IMAGES WITH CONDITION LABELS from D.Trials SCRIPT


%%%%%%%%%%           SECTION 0 - DEFINE VARIABLES          %%%%%%%%%%%%

% names to construct each participants mat file and nii image file paths
meg_prefix = 'ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-';
meg_suffix = '_ses-1_task-dur_run-01_meg.mat';
nii_file_prefix = 'ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-';
nii_file_suffix = '_ses-1_task-dur_run-01_meg_1_t0_600_f0_256_';

cd(image_dir)

meg_path = fullfile(base_filepath, 'derivatives', 'preprocessed_files');

for i = 1:numel(p_names)
    
    current_subject = p_names{i};

    meg_filename = [meg_prefix, current_subject, meg_suffix];
    meg_filepath = fullfile(meg_path, meg_filename);
 
    % --- Load the MEG mat file to get condition labels ---
    load(meg_filepath);
    
    % Get the list of condition labels from the object
    condition_labels = {D.trials.label};

    % --- Loop through conditions and rename files ---
    for j = 1:numel(condition_labels)
        
        label = condition_labels{j};
        
        % Construct the old and new filenames
        old_filename = sprintf('%s%s%s%d.nii', nii_file_prefix, current_subject, nii_file_suffix, j);
        old_filepath = fullfile(image_dir, old_filename);

        new_filename = [current_subject, '_ONSET_', label, '.nii'];
        new_filepath = fullfile(image_dir, new_filename);

        copyfile(old_filepath, new_filepath);
        
    end
end

%% RENAME SOURCE RECONSTRUCTED IMAGES WITH CONDITION LABELS from D.Trials SCRIPT

%%%%%%%%%%           SECTION 0 - DEFINE VARIABLES          %%%%%%%%%%%%

% names to construct each participants mat file and nii image file paths
meg_prefix = 'ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-';
meg_suffix = '_ses-1_task-dur_run-01_meg.mat';
nii_file_prefix = 'ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-';
nii_file_suffix = '_ses-1_task-dur_run-01_meg_1_t0_600_f0_256_';

cd(image_dir)

meg_path = fullfile(base_filepath, 'derivatives', 'preprocessed_files');

for i = 1:numel(p_names)
    
    current_subject = p_names{i};
  
    meg_filename = [meg_prefix, current_subject, meg_suffix];
    meg_filepath = fullfile(meg_path, meg_filename);
 
    % --- Load the MEG mat file to get condition labels ---
    load(meg_filepath);
    
    % Get the list of condition labels from the object
    condition_labels = {D.trials.label};

    % --- Loop through conditions and rename files ---
    for j = 1:numel(condition_labels)
        
        label = condition_labels{j};
        
        % Construct the old and new filenames
        old_filename = sprintf('%s%s%s%d.nii', nii_file_prefix, current_subject, nii_file_suffix, j);
        old_filepath = fullfile(meg_path, old_filename);

        new_filename = [current_subject, '_ONSET_', label, '.nii'];
        new_filepath = fullfile(image_dir, new_filename);

        copyfile(old_filepath, new_filepath);
        
    end
end

%% RUN GLM CONTRASTS OF SOURCE RECONSTRUCTED IMAGES - 1 GLM per participant

relevant = '_relevant';
irrelevant = '_irrelevant';

clear matlabbatch;

for c = 1:length(p_names)

    new_folder = [p_names{c} '_GLM'];
    GLM_dir = fullfile(image_dir, new_folder);  
    
    if exist(GLM_dir, 'dir')
        rmdir(GLM_dir, 's'); 
    end
    mkdir(GLM_dir);

                % --- STEP 1: define factorial design --- 
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {GLM_dir}; %dir to save SPM.mat file
    
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'face_other';
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; %1 = dependent;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0; %1 = equal 
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
    
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'relevance';
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1; %1 = dependent;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0; %1 = equal 
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
    
    % Get all .nii files for the current participant
    participant_pattern = [p_names{c}, '*.nii'];
    files_struct = dir(fullfile(image_dir, participant_pattern));
    nii_files = {files_struct.name}'; % Get filenames as a cell array

    % 2. Create logical masks for each condition
    is_object = contains(nii_files, '_object');
    is_face   = contains(nii_files, '_face');
    is_letter = contains(nii_files, '_letter');
    is_false  = contains(nii_files, '_false');
    is_other = is_object | is_letter | is_false;
    is_relevant   = contains(nii_files, relevant);
    is_irrelevant = contains(nii_files, irrelevant);

    % Cell 1: Face - Relevant
    clean_files_fr = nii_files(is_face & is_relevant);
    
    % Cell 2: Face - Irrelevant
    clean_files_fi = nii_files(is_face & is_irrelevant);
    
    % Cell 3: Other - Relevant
    clean_files_or = nii_files(is_other & is_relevant);

    % Cell 4: Other - Irrelevant
    clean_files_oi = nii_files(is_other & is_irrelevant);

    % --- Populate the cells in the matlabbatch ---
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1 1]; % Type 1, Relevance 1
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = fullfile(image_dir, clean_files_fr);
    fprintf(' -> Found %d files for Face-Relevant\n', numel(clean_files_fr));
    
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1 2]; % Type 1, Relevance 2
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = fullfile(image_dir, clean_files_fi);
    fprintf(' -> Found %d files for Face-Irrelevant\n', numel(clean_files_fi));
    
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [2 1]; % Type 2, Relevance 1
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = fullfile(image_dir, clean_files_or);
    fprintf(' -> Found %d files for Other-Relevant\n', numel(clean_files_or));
    
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2 2]; % Type 2, Relevance 2
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = fullfile(image_dir, clean_files_oi);
    fprintf(' -> Found %d files for Other-Irrelevant\n', numel(clean_files_oi));

    %rest of the settings
    
    matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % --- STEP 2: Run Model Estimation ---
    spmmat_file = fullfile(GLM_dir, 'SPM.mat'); %change this line if running on all participants in one go
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {spmmat_file};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    spm_jobman('run', matlabbatch);
    
    clear matlabbatch

end

%% NOW EXTRACT PEAK MNI COORDINATES FROM GLM CONTRASTS

%                    --- Configuration ---

dcm_folder = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_onset');

if ~exist(dcm_folder, 'dir')
    mkdir(dcm_folder);
end

outfile = fullfile(dcm_folder, 'indiv_mni_coords.m');

% ROI labels in order 
roi_labels = {'L_V1', 'R_V1', 'L_FFA', 'R_FFA','L_PFC', 'R_PFC'};

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');
spm_jobman('initcfg');
clear matlabbatch;

%                    --- Run Contrasts ---

for c = 1:length(p_names)

    %folder for this participant
    GLM_dir = fullfile(image_dir, [p_names{c} '_GLM']); 
    filename = fullfile(GLM_dir, 'SPM.mat'); 
    
    for i = 1:length(roi_labels)
    
        %path to mask for this roi
        mask_dir = '/data/gpfs/projects/punim2118/Cogitate_DCM/roi_masks/jubrain_masks';
        mask_file = fullfile(mask_dir, [roi_labels{i} '.nii']);
    
        %contrast manager
        matlabbatch{1}.spm.stats.con.spmmat = {filename};
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'Main effect of face_other';
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 1 -1 -1];
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 0;
        
        %create results 
        matlabbatch{2}.spm.stats.results.spmmat = {filename};
        matlabbatch{2}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{2}.spm.stats.results.conspec.contrasts = 2; % specify the contrast index -- im pretty sure this is 2? see SPM.xCon 
        matlabbatch{2}.spm.stats.results.conspec.threshdesc = 'none'; %none is ok because multiple comparisons isn't a problem here - we just want to find where there is the most robust activation
        matlabbatch{2}.spm.stats.results.conspec.thresh = 1.0; %note 0.05 is 95%. setting to 1 means spm will get any voxel that is different/positive for a participant 
        matlabbatch{2}.spm.stats.results.conspec.extent = 0;
        matlabbatch{2}.spm.stats.results.conspec.conjunction = 1;
        %matlabbatch{2}.spm.stats.results.conspec.mask.none = 1; %if no mask 
        matlabbatch{2}.spm.stats.results.conspec.mask.image.name = {mask_file}; %filepath to mask .nii
        matlabbatch{2}.spm.stats.results.conspec.mask.image.mtype = 0; %0 is inclusive mask 
        matlabbatch{2}.spm.stats.results.units = 1;
        matlabbatch{2}.spm.stats.results.export{1}.ps = true;
        matlabbatch{2}.spm.stats.results.export{2}.jpg = true;
        matlabbatch{2}.spm.stats.results.export{3}.csv = true;
        
        startTime = now; % Get a timestamp -- used for renaming files

        spm_jobman('run', matlabbatch);

        % --- RENAME THE OUTPUT CSV FILE ---
        % Find all CSV files matching the SPM default format in the GLM directory
        spm_csv_files = dir(fullfile(GLM_dir, 'spm_*.csv'));
        [~, latest_idx] = max([spm_csv_files.datenum]); %get latest file        
        old_filename = fullfile(GLM_dir, spm_csv_files(latest_idx).name);
        new_filename = fullfile(GLM_dir, [roi_labels{i} '.csv']);
        movefile(old_filename, new_filename); %now rename        

        clear matlabbatch

    end
end

%                    --- Get Peak MNI Coords ---


indiv_mni_coords = cell(length(p_names), 2);

for p = 1:length(p_names)
    
    this_p = p_names{p};
    
    GLM_dir = fullfile(image_dir, [this_p '_GLM']);
    
    % Create a 3xN matrix to hold the coordinates for this participant
    % N is the number of ROIs. Each column will be one [x; y; z] coordinate.
    participant_coords_matrix = nan(3, length(roi_labels));

    for r = 1:length(roi_labels)
        
        this_roi = roi_labels{r};
        csv_file = fullfile(GLM_dir, [this_roi '.csv']);
    
        try
            data = readmatrix(csv_file, 'HeaderLines', 2);
        catch
            data = [];
        end

        if ~isempty(data)
            % The x,y,z coordinates are in columns 12, 13, and 14 and we
            % only want the first row - the peak
            peak_coord = data(1, 12:14)'; 
            
        else
            peak_coord = [NaN; NaN; NaN];
        end

        participant_coords_matrix(:, r) = peak_coord;

    end 
    
    % Store this participant's 3x6 coordinate matrix in the main cell array
    indiv_mni_coords{p, 1} = this_p;
    indiv_mni_coords{p, 2} = participant_coords_matrix;
    
end
    
% Save the final cell array to a .mat file

output_mat_file = fullfile(dcm_folder, 'indiv_peak_coords.mat');

save(output_mat_file, 'indiv_mni_coords');    
    
   

%% RUN GLM CONTRASTS OF SOURCE RECONSTRUCTED IMAGES at group level for missing peaks 

% CA101 had no value for L_PFC so we will find the group mean value for
% this participant to use as their PFC coordinate

clear matlabbatch
face = '_face';
other = '_other';

relevant = '_relevant';
irrelevant = '_irrelevant';

output_dir    = fullfile(image_dir, 'GROUP_AVERAGE_L_PFC');
mask_file    = '/data/gpfs/projects/punim2118/Cogitate_DCM/roi_masks/jubrain_masks/L_PFC.nii'; % change mask file here!!!

if exist(output_dir, 'dir')
    delete(fullfile(output_dir, '*.*')); % Deletes the old SPM.mat and images safely
else
    mkdir(output_dir);
end
FR = {}; 
FI = {}; 
OR = {}; 
OI = {};

p_names(strcmp(p_names, 'CA101')) = [];

% get images
for i = 1:length(p_names)
    
    sub = p_names{i};
    
    pattern = fullfile(image_dir, [sub, '*.nii']);
    d = dir(pattern);
    all_files = fullfile({d.folder}, {d.name})'; % Get full paths
       
    % logical indexing to find specific conditions
    is_face = contains(all_files, face);
    is_other   = contains(all_files, other);
    is_rel = contains(all_files, relevant);
    is_irrel  = contains(all_files, irrelevant);
    
    % Sort into lists
    FR = [FR; all_files(is_face & is_rel)];
    FI = [FI; all_files(is_face & is_irrel)];
    OR = [OR; all_files(is_other   & is_rel)];
    OI = [OI; all_files(is_other   & is_irrel)];
end

            % --- STEP 1: define factorial design --- 

matlabbatch{1}.spm.stats.factorial_design.dir = {output_dir}; %dir to save SPM.mat file

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'face_other';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; %1 = dependent;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0; %1 = equal 
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'relevance';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1; %1 = dependent;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0; %1 = equal 
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;


% --- Populate the cells in the matlabbatch ---
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1 1]; % Type 1, Relevance 1
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = FR;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1 2]; % Type 1, Relevance 2
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = FI;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [2 1]; % Type 2, Relevance 1
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = OR;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2 2]; % Type 2, Relevance 2
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = OI;


matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spmmat_file = fullfile(output_dir, 'SPM.mat'); %change this line if running on all participants in one go

matlabbatch{2}.spm.stats.fmri_est.spmmat = {spmmat_file};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat = {spmmat_file};
matlabbatch{3}.spm.stats.con.delete = 1;
matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Main effect of face_other';
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1 1 -1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';

        %create results 
matlabbatch{4}.spm.stats.results.spmmat = {spmmat_file};
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1; % specify the contrast index -- im pretty sure this is 2? see SPM.xCon 
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none'; %none is ok because multiple comparisons isn't a problem here - we just want to find where there is the most robust activation
matlabbatch{4}.spm.stats.results.conspec.thresh = 1.0; %note 0.05 is 95%. setting to 1 means spm will get any voxel that is different/positive for a participant 
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1; 
%matlabbatch{4}.spm.stats.results.conspec.mask.none = 1; %if no mask 
matlabbatch{4}.spm.stats.results.conspec.mask.image.name = {mask_file}; %filepath to mask .nii
matlabbatch{4}.spm.stats.results.conspec.mask.image.mtype = 0; %0 is inclusive mask 
matlabbatch{4}.spm.stats.results.units = 1;

matlabbatch{4}.spm.stats.results.export{1}.jpg = true;
matlabbatch{4}.spm.stats.results.export{2}.csv = true;

spm_jobman('run', matlabbatch);

clear matlabbatch
