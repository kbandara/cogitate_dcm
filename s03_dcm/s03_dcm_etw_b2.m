function s03_dcm_etw_b2(time_window, p_names_idx)

    addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
    spm('defaults', 'eeg');
    spm_jobman('initcfg');
    clear matlabbatch;

    %{

    NOTE: this is a function that takes a time window/participant from the SLURM scheduler to run in multiple
    (parallel) job arrays 

    THIS SCRIPT IS IDENTICAL TO b1 BUT FOR BATCH 2 
    
    kav 2025
    unimelb
    
    this script performs an expanding time window analysis of DCMs. 
    It runs DCMs iteratively in expanding windows (0-100ms, 0-200ms, etc.) 
    from 0ms up to 1000ms. For each time window, it 
    generates DCMs for all participants, then collates them into a GCM, 
    fits the GCM, and saves it.
    
    %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        
%                            TESTING BLOCK
%                      uncomment to run on slurm
%
%   time_window = 1; p_names_idx = 1; 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    time_window_increments = 100:100:1000;

    current_time_window_end = time_window_increments(time_window);

    fprintf('\n\n--- processing time window: 0 to %dms ---\n\n', current_time_window_end);

    addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
    spm('defaults', 'eeg');
    spm_jobman('initcfg');
    clear matlabbatch;
    
    base_filepath  = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/'; 
    
    %dir containing spmeeg files 
    spm_path = fullfile(base_filepath, 'derivatives', 'preprocessed_files');
   
    % Load participant names
    P_Names_struct = load(fullfile(base_filepath, 'p_names.mat'));  
    P_Names_table =  P_Names_struct.uniqueToFirstList; %{'CA107'}; 
    P_Names_str = P_Names_table;
    p_names     = cellstr(P_Names_str);  % 1x52 cell array of char

    PP = p_names_idx;
    fprintf('\n\n--- TW: 0–%dms | Participant: %s (%d/%d) ---\n\n', ...
        current_time_window_end, p_names{PP}, PP, length(p_names));
    
    % output directory for DCM results
    output_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw'); 
    
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    dcm_filename = ['DCM_' p_names{PP} '_' num2str(current_time_window_end) 'ms.mat'];
    dcm_filepath = fullfile(output_dir, dcm_filename);

    %load each participants peak mni coordinates 
    indiv_mni_coords = load(fullfile(spm_path, 'dcm_onset', 'indiv_peak_coords.mat'));
    indiv_mni_coords = indiv_mni_coords.indiv_mni_coords;
    
    % ROI labels in order 
    roi_labels = {'L_V1', 'R_V1', 'L_FG', 'R_FG','L_PFC', 'R_PFC'};
   
%% DCM 


    % Data filename
    %--------------------------------------------------------------------------
    filename = ['ronly_meg_bc_lf_ra_r_ms_re_faces_a_d_spmeeg_sub-' p_names{PP} '_ses-1_task-dur_run-01_meg.mat'];
    spmeeg_file = fullfile(spm_path, filename);

    DCM.xY.Dfile = spmeeg_file;

    % Parameters and options used for setting up model
    %--------------------------------------------------------------------------
    DCM.options.analysis = 'ERP'; % analyze evoked responses
    DCM.options.model    = 'ERP'; % ERP model
    DCM.options.spatial  = 'ECD'; % spatial model
    DCM.options.trials   = [2, 4];     % index of ERPs within ERP/ERF file - 2 is the index for irrelevant face 1000ms trials
    DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
    DCM.options.Tdcm(2)  = current_time_window_end;   % end of peri-stimulus time to be modelled
    DCM.options.Nmodes   = 8;     % nr of modes for data selection 
    DCM.options.D        = 1;     % downsampling factor
    DCM.options.h        = 1;     % nr of DCT components %'1' for 'detrend' i.e. model the mean       
    DCM.options.onset    = 64;    % selection of onset in ms (prior mean); 64 is SPM default          
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
    
    % Add in the contrasts here - length of this vector should match amount of conditions specified above 

    %OR to model just the A matrix leave between trial effects as empty as
    %below
    DCM.xU.X = [1; 1]; 
    DCM.xU.name = {};

    %--------------------------------------------------------------------------
    % Invert and save
    %--------------------------------------------------------------------------
    DCM.name = dcm_filename; 
    
    DCM = spm_dcm_erp(DCM);

    disp(dcm_filename);
    save(dcm_filepath, 'DCM');

end