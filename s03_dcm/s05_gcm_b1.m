%s02_gcm_b1

%kav bandara unimelb 2025

clear; clc;
addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');

% 1. Setup paths
base_filepath  = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH1/'; 
output_dir     = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw'); 

% 2. Load participant names
P_Names_struct = load(fullfile(base_filepath, 'participants.mat'));
P_Names_table  = P_Names_struct.participants;
p_names        = P_Names_table{:, 1}'; 

% Define the end points for each sliding window in milliseconds
time_window_increments = 100:100:1000;

% for onset
for t = 1:length(time_window_increments)
    
    current_time_window_end = time_window_increments(t);
    fprintf('\n\n--- Processing Time Window: 0 to %dms ---\n\n', current_time_window_end);

    GCM = {};
    for PP = 1:length(p_names)    
        dcm_filename = ['DCM_' p_names{PP} '_' num2str(current_time_window_end) 'ms.mat']; 
        dcm_path = fullfile(output_dir, dcm_filename);
        GCM = [GCM, dcm_path];
    end
    GCM = GCM';
    
    % Write results for this time window
    gcm_filename = ['GCM_' num2str(current_time_window_end) 'ms.mat'];
    fprintf('Saving fitted GCM: %s\n', gcm_filename);
    save(fullfile(output_dir, gcm_filename), 'GCM');

end

clear GCM gcm_filename

%%
% FOR OFFSET

time_window_increments = 100:100:600;
output_dir     = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw', 'offset'); 

for t = 1:length(time_window_increments)
    
    current_time_window_end = time_window_increments(t);
    fprintf('\n\n--- Processing Time Window: 0 to %dms ---\n\n', current_time_window_end);

    GCM = {};
    for PP = 1:length(p_names)    
        dcm_filename = ['DCM_offset_' p_names{PP} '_' num2str(current_time_window_end) 'ms.mat']; 
        dcm_path = fullfile(output_dir, dcm_filename);
        GCM = [GCM, dcm_path];
    end
    GCM = GCM';
    
    % Write results for this time window
    gcm_filename = ['GCM_OFFSET_' num2str(current_time_window_end) 'ms.mat'];
    fprintf('Saving fitted GCM: %s\n', gcm_filename);
    save(fullfile(output_dir, gcm_filename), 'GCM');

end



