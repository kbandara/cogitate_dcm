% GNWT PEB tests - batch 1 
% kav 2025

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');
spm_jobman('initcfg');

clear all;

%% Setup paths
base_filepath = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH1/';
spm_path = fullfile(base_filepath, 'derivatives', 'preprocessed_files');

P_Names_struct = load(fullfile(base_filepath, 'participants.mat'));
P_Names_table = P_Names_struct.participants;
p_names = P_Names_table{:, 1}';
p_names = p_names(~strcmp(p_names, 'CA136')); %exclude CA136 for too few hits

indiv_mni_coords = load(fullfile(spm_path, 'dcm_onset', 'indiv_peak_coords.mat'));
indiv_mni_coords = indiv_mni_coords.indiv_mni_coords;

roi_labels = {'L_V1', 'R_V1', 'L_FFA', 'R_FFA', 'L_PFC', 'R_PFC'};

output_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw');

offset_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw', 'offset');

post_offset_increments = 300:100:500;  % 300-500ms post-offset

cd(output_dir)

time_window_increments = 100:100:1000;

%% ONSET ANALYSES 

for t = 1:length(time_window_increments)
    
    current_time_window_end = time_window_increments(t);
    
    fprintf('\n=== Processing %dms ===\n', current_time_window_end);
    
    gcm_file = fullfile(output_dir, ['GCM_' num2str(current_time_window_end) '.mat']);
    
    load(gcm_file, 'GCM');
    n_total = length(GCM);

    % Specify PEB model settings
    M   = struct();
    M.Q = 'all'; 
    M.X = ones(n_total, 1);  % model mean
    
    field = {'A'};
    
    % Estimate model
    PEB = spm_dcm_peb(GCM, M, field);
    
    peb_filename = ['PEB_combined_' num2str(current_time_window_end) 'ms.mat'];
    save(fullfile(output_dir, peb_filename), 'PEB');
    
    clear GCM PEB
end

% Run BMC and save

for t = 1:length(time_window_increments)
    current_time_window_end = time_window_increments(t);
    
    peb_filename = fullfile(output_dir, ['PEB_combined_' num2str(current_time_window_end) 'ms.mat']);
    bma_filename = fullfile(output_dir, ['BMA_combined_' num2str(current_time_window_end) 'ms.mat']);
    
    load(peb_filename, 'PEB');

    BMA = spm_dcm_peb_bmc(PEB);
        
    % Save the BMA results
    save(bma_filename, 'BMA');
    fprintf('Saved %s\n', bma_filename);
    
    clear PEB BMA
end

% Store posterior probabilities over time
GNWT_prob = nan(1, length(time_window_increments));
GNWT_df   = nan(1, length(time_window_increments));

GNWT_onset = [300, 500];
stimulus_end = 1000; 

for t = 1:length(time_window_increments)
    
    current_time = time_window_increments(t);
    
 
    % Load the computed PEB for this time window
    peb_filename = fullfile(output_dir, ['PEB_combined_' num2str(current_time) 'ms.mat']);
    gcm_filename = fullfile(output_dir, ['GCM_combined_' num2str(current_time) 'ms.mat']);
       
    load(peb_filename, 'PEB');
    load(gcm_filename, 'GCM');
    
    % Use first subject's DCM as a template
    DCM_template_raw = GCM{1};

    if ischar(DCM_template_raw) || isstring(DCM_template_raw)
        loaded_file = load(DCM_template_raw);
        DCM_template = loaded_file.DCM;
    else
        DCM_template = DCM_template_raw;
    end

    if isfield(DCM_template,'M'), DCM_template = rmfield(DCM_template,'M'); end
    
    % Model 1: FULL
    DCM_Full = DCM_template; 
    
    % Model 2: No-GNWT (Remove PFC Feedback)

    is_onset = (current_time >= GNWT_onset(1) && current_time <= GNWT_onset(2));
  
    if is_onset 
        DCM_GNWT = DCM_template;
        
        % Node order: 1=LV1, 2=RV1, 3=LFFA, 4=RFFA, 5=LPFC, 6=RPFC
        DCM_GNWT.A{2}(3, 5) = 0; % LPFC -> LFFA
        DCM_GNWT.A{2}(4, 6) = 0; % RPFC -> RFFA

        % Run BMC
        [~, BMR_GNWT] = spm_dcm_peb_bmc(PEB, {DCM_Full, DCM_GNWT});
        
        % Extract Free Energy Difference
        F_full     = BMR_GNWT{1}.F;
        F_reduced  = BMR_GNWT{2}.F;
        
        dF = F_full - F_reduced;
        
        % Store F and pps
        GNWT_df(t) = dF;
        GNWT_prob(t) = 1 / (1 + exp(-dF));
        
        fprintf('  GNWT: dF = %.2f, Pp = %.3f\n', dF, GNWT_prob(t));
    end

    clear PEB GCM DCM_template DCM_Full DCM_GNWT BMR_GNWT
end

% Save results
save(fullfile(output_dir, 'hypothesis_test_results.mat'), ...
     'time_window_increments', 'GNWT_prob', 'GNWT_df');

%% OFFSET PEB: Combine durations and run PEB at each offset time 

cd(offset_dir)

GNWT_prob_offset = nan(1, length(post_offset_increments));
GNWT_df_offset   = nan(1, length(post_offset_increments));

for t = 1:length(post_offset_increments)
    
    post_offset_time = post_offset_increments(t);

    gcm_offset_file = fullfile(offset_dir, ['GCM_OFFSET_' num2str(post_offset_time) 'ms.mat']);
    load(gcm_offset_file, 'GCM');
    n_offset= length(GCM);
    
    M = struct();
    M.Q = 'all';
    M.X = ones(n_offset, 1);  % Just model the mean 
    
    field = {'A'};
    
    PEB = spm_dcm_peb(GCM, M, field);
    
    peb_filename = ['PEB_OFFSET_' num2str(post_offset_time) 'ms_postoffset.mat'];
    save(fullfile(offset_dir, peb_filename), 'PEB');
    

    % Get template DCM
    DCM_template = GCM{1};
    if isfield(DCM_template, 'M')
        DCM_template = rmfield(DCM_template, 'M');
    end
    
    % Model 1: FULL (Baseline)
    DCM_Full = DCM_template;
    
    % Model 2: No-GNWT (Remove PFC Feedback)
    DCM_GNWT = DCM_template;
    
    % Turn OFF Backward connections from PFC to posterior areas
    DCM_GNWT.A{2}(3, 5) = 0;  % LPFC -> LFFA
    DCM_GNWT.A{2}(4, 6) = 0;  % RPFC -> RFFA
    
    % Run BMC
    [~, BMR_GNWT] = spm_dcm_peb_bmc(PEB, {DCM_Full, DCM_GNWT});
    
    % Extract Free Energy Difference
    F_full = BMR_GNWT{1}.F;
    F_reduced = BMR_GNWT{2}.F;
    
    dF = F_full - F_reduced;
    
    % Store results
    GNWT_df_offset(t) = dF;
    GNWT_prob_offset(t) = 1 / (1 + exp(-dF));
    
    fprintf('  GNWT: dF = %.2f, P(PFC necessary) = %.3f\n', dF, GNWT_prob_offset(t));
    
    clear GCM GCM_short GCM_long GCM_combined PEB
end

% Save results
save(fullfile(offset_dir, 'GNWT_offset_combined_results.mat'), ...
     'post_offset_increments', 'GNWT_prob_offset', 'GNWT_df_offset');
