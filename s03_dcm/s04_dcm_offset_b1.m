function s02_dcm_offset_b1(time_window, p_names_idx)

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');
spm_jobman('initcfg');


    %{
    kav 2025
    unimelb
        
    % this takes a time window from the SLURM scheduler to run in multiple
    % job arrays - i.e. parallel

    This script performs DCM across an expanding time window for OFFSET only!
    
    It runs DCMs/GCMs iteratively in expanding windows (0-100ms, 0-200ms, etc.) 
    from 0ms up to 600ms post offset

    %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        
%                            TESTING BLOCK
%                        uncomment to run on slurm

%   time_window = 1; p_names_idx = 1; 

%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT: Variables and paths setup

base_filepath  = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH1/'; 

%dir containing spmeeg files 
filepath = fullfile(base_filepath, 'derivatives', 'preprocessed_files');
cd (filepath);

%Load participant names
P_Names_struct = load(fullfile(base_filepath, 'participants.mat'));
P_Names_table =  P_Names_struct.participants;
p_names = P_Names_table{:, 1}';%{'CA107'}; 


%load each participants peak mni coordinates 
indiv_mni_coords = load(fullfile(filepath, 'dcm_onset', 'indiv_peak_coords.mat'));
indiv_mni_coords = indiv_mni_coords.indiv_mni_coords;

% ROI labels in order 
roi_labels = {'L_V1', 'R_V1', 'L_FFA', 'R_FFA','L_PFC', 'R_PFC'};

% output directory for DCM results
output_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw', 'offset'); 

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

PP = p_names_idx;
% Define the end points for each sliding window in milliseconds
time_window_increments = 100:100:600;
current_time_window_end = time_window_increments(time_window);

fprintf('\n\n--- TW: 0–%dms | Participant: %s (%d/%d) ---\n\n', ...
    current_time_window_end, p_names{PP}, PP, length(p_names));

dcm_filename = ['DCM_offset_' p_names{PP} '_' num2str(current_time_window_end) 'ms.mat'];
%}

%{ 
for move = 1:length(p_names)
    
    source_file = fullfile(filepath, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' p_names{move} '_ses-1_task-dur_run-01_meg.mat']);
    
    spm_jobman('initcfg');

    % 900ms is -100ms relative to Short Offset
    matlabbatch{1}.spm.meeg.preproc.crop.D = {source_file};
    matlabbatch{1}.spm.meeg.preproc.crop.timewin = [900 1500]; 
    matlabbatch{1}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
    matlabbatch{1}.spm.meeg.preproc.crop.channels{1}.all = 'all';
    matlabbatch{1}.spm.meeg.preproc.crop.prefix = 'short_offset_';
    spm_jobman('run', matlabbatch); clear matlabbatch;

    % 1400ms is -100ms relative to Long Offset
    matlabbatch{1}.spm.meeg.preproc.crop.D = {source_file};
    matlabbatch{1}.spm.meeg.preproc.crop.timewin = [1400 2000]; 
    matlabbatch{1}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
    matlabbatch{1}.spm.meeg.preproc.crop.channels{1}.all = 'all';
    matlabbatch{1}.spm.meeg.preproc.crop.prefix = 'long_offset_';
    spm_jobman('run', matlabbatch); clear matlabbatch;

    % shift both to start from -100ms 
    D_short = spm_eeg_load(fullfile(filepath, ['short_offset_ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' p_names{move} '_ses-1_task-dur_run-01_meg.mat']));
    D_short = timeonset(D_short, -0.1);
    D_short.save();
    
    D_long = spm_eeg_load(fullfile(filepath, ['long_offset_ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' p_names{move} '_ses-1_task-dur_run-01_meg.mat']));
    D_long = timeonset(D_long, -0.1);  
    D_long.save();

    % merge
    matlabbatch{1}.spm.meeg.preproc.merge.D = {D_short.fullfile; D_long.fullfile};
    matlabbatch{1}.spm.meeg.preproc.merge.recode(1).file = D_short.fname;
    matlabbatch{1}.spm.meeg.preproc.merge.recode(1).labelorg = '.*';
    matlabbatch{1}.spm.meeg.preproc.merge.recode(1).labelnew = 'Short_#labelorg#';
    matlabbatch{1}.spm.meeg.preproc.merge.recode(2).file = D_long.fname;
    matlabbatch{1}.spm.meeg.preproc.merge.recode(2).labelorg = '.*';
    matlabbatch{1}.spm.meeg.preproc.merge.recode(2).labelnew = 'Long_#labelorg#';
    matlabbatch{1}.spm.meeg.preproc.merge.prefix = 'OFFSET_COMBINED_';
    spm_jobman('run', matlabbatch); clear matlabbatch;

    % delete temps
    delete(D_short.fname); delete(D_short.fnamedat);
    delete(D_long.fname); delete(D_long.fnamedat);

    fprintf('Created OFFSET_COMBINED for %s\n', p_names{move});
end
%}
%% DCM 


        % Data filename
        %--------------------------------------------------------------------------
        filename = ['OFFSET_COMBINED_short_offset_ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' p_names{PP} '_ses-1_task-dur_run-01_meg.mat'];
        spmeeg_file = fullfile(output_dir, filename);
    
        DCM.xY.Dfile = spmeeg_file;
    
        % Parameters and options used for setting up model
        %--------------------------------------------------------------------------
        DCM.options.analysis = 'ERP'; % analyze evoked responses
        DCM.options.model    = 'ERP'; % ERP model
        DCM.options.spatial  = 'ECD'; % spatial model
        DCM.options.trials   = [2, 12];     % index of condition - because the merging processed doubled up we must select the correctly baseline corrected trials from the first 8 conds (baseline corrected short trials) and then from the last 8 (the long trials)
        DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
        DCM.options.Tdcm(2)  = current_time_window_end;   % end of peri-stimulus time to be modelled
        DCM.options.Nmodes   = 8;     % nr of modes for data selection 
        DCM.options.D        = 1;     % downsampling factor
        DCM.options.h        = 1;     % nr of DCT components %'1' for 'detrend' i.e. model the mean       
        DCM.options.onset    = 64;    % selection of onset in ms (prior mean); 64 is SPM default PLUS 1000 because it's relative to epoch start          
        DCM.options.dur      = 16;    % and dispersion (sd); 16ms is SPM default
        
        DCM.xY.modality = 'MEGPLANAR'; % Specify modality -- MEGPLANAR IS COMBINED
    
        %--------------------------------------------------------------------------
        % Data and spatial model
        %--------------------------------------------------------------------------
        DCM  = spm_dcm_erp_data(DCM);
        
        %--------------------------------------------------------------------------
        % Location priors for dipoles
        %--------------------------------------------------------------------------
        DCM.Lpos  = indiv_mni_coords{PP, 2}; %only select the second col which contains the actual coordinates
        DCM.Sname = roi_labels;
        Nareas    = size(DCM.Lpos,2); %number of dipoles, used to create A/B/C matrices 
       
        %--------------------------------------------------------------------------
        % Spatial model
        %--------------------------------------------------------------------------
        DCM = spm_dcm_erp_dipfit(DCM);
        
        %--------------------------------------------------------------------------
        % Specify connectivity model
        %--------------------------------------------------------------------------
        
        %Setup forward and backward connection matrices
        DCM.A{1} = zeros(Nareas,Nareas); %forward matrix
        DCM.A{2} = zeros(Nareas,Nareas); %backward matrix
        DCM.A{3} = zeros(Nareas,Nareas); %lateral matrix
        
        DCM.B = {}; %let B matrix be an empty cell array 
        %DCM.B{1} = zeros(Nareas,Nareas); %Setup the empty B matrix 
        
        DCM.C = zeros(Nareas); %setup empty C matrix
    
        %Forwards connections
        DCM.A{1}(3,1) = 1; % Locc - Lffa
        DCM.A{1}(5,1) = 1; % Locc - Lpfc
        DCM.A{1}(4,2) = 1; % Rocc - Rffa
        DCM.A{1}(6,2) = 1; % Rocc - Rpfc
        DCM.A{1}(5,3) = 1; % Lffa - Lpfc
        DCM.A{1}(6,4) = 1; % Rffa - Rpfc
        
        %Backward connections
        DCM.A{2}(1,3) = 1; % Lffa - Locc
        DCM.A{2}(1,5) = 1; % Lpfc - Locc
        DCM.A{2}(2,4) = 1; % Rffa - Rocc
        DCM.A{2}(2,6) = 1; % Rpfc - Rocc
        DCM.A{2}(3,5) = 1; % Lpfc - Lffa
        DCM.A{2}(4,6) = 1; % Rpfc - Rffa
        
        %Lateral connections
        DCM.A{3}(2,1) = 1; %Locc - Rocc
        DCM.A{3}(1,2) = 1; %Rocc - Locc
        DCM.A{3}(4,3) = 1; %Lffa - Rffa
        DCM.A{3}(3,4) = 1; %Rffa - Lffa
        DCM.A{3}(6,5) = 1; %Lpfc - Rpfc
        DCM.A{3}(5,6) = 1; %Rpfc - Lpfc
            
      
        %DCM C matrix simply specified where the original input source is in the brain (V1 for this exp)
        DCM.C = [1; 1; 0; 0; 0; 0]; %Input sources (LV1 and RV1)
       
        
        %--------------------------------------------------------------------------
        % Between trial effects 
        %--------------------------------------------------------------------------

        DCM.xU.X = [1; 1]; %model mean of both trials
        DCM.xU.name = {};
    
        %--------------------------------------------------------------------------
        % Invert and save
        %--------------------------------------------------------------------------
        DCM.name = dcm_filename; 
        
        DCM = spm_dcm_erp(DCM);
    
        disp(dcm_filename);

        save(fullfile(output_dir, dcm_filename), 'DCM');

end

