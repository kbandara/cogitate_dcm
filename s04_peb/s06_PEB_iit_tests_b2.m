% PEB analysis for IIT - batch 2 
% kav bandara unimelb 2025

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');
spm_jobman('initcfg');
clear matlabbatch;

clear all; clc;

% INPUT: Variables and paths setup

base_filepath  = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/'; 

% Load participant names
P_Names_struct = load(fullfile(base_filepath, 'p_names.mat'));  
P_Names_table =  P_Names_struct.uniqueToFirstList; %{'CA107'}; 
P_Names_str = P_Names_table;
p_names     = cellstr(P_Names_str);  
p_names = p_names(~strcmp(p_names, 'CB082')); %exclude CB082 for too few hits

% output directory for DCM results
output_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw'); 
cd(output_dir)

% Define the end points for each sliding window in milliseconds
time_window_increments = 300:100:1000;

%% extrinsic tests 

% Store results
IIT_prob = nan(1, length(time_window_increments));
IIT_df = nan(1, length(time_window_increments));

% Run Combined PEB for each time window

for t = 1:length(time_window_increments)
    
    current_time = time_window_increments(t);
    
    fprintf('\n=== Processing %dms ===\n', current_time);

    gcm_file = fullfile(output_dir, ['GCM_' num2str(current_time) 'ms.mat']);

    load(gcm_file, 'GCM');
    n_total = length(GCM);

    save(fullfile(output_dir, ['GCM_combined_' num2str(current_time) 'ms.mat']), 'GCM');
    
    % Run PEB
    M = struct();
    M.Q = 'all';
    M.X = ones(n_total, 1);
    field = {'A'};
    
    PEB = spm_dcm_peb(GCM, M, field);
    save(fullfile(output_dir, ['PEB_combined_' num2str(current_time) 'ms.mat']), 'PEB');
    
    if current_time >= 100
        
        DCM_template_raw = GCM{1};
        loaded_file = load(DCM_template_raw);
        DCM_template = loaded_file.DCM;

        if isfield(DCM_template, 'M')
            DCM_template = rmfield(DCM_template, 'M');
        end
        
        % Full model
        DCM_Full = DCM_template;
        
        % No-IIT model (remove FG to V1 and lateral posterior connections
        DCM_IIT = DCM_template;
        DCM_IIT.A{2}(1, 3) = 0;  % L_FFA -> L_V1
        DCM_IIT.A{2}(2, 4) = 0;  % R_FFA -> R_V1

        % Run BMC
        [~, BMR_IIT] = spm_dcm_peb_bmc(PEB, {DCM_Full, DCM_IIT});
        
        % Extract results
        F_full = BMR_IIT{1}.F;
        F_reduced = BMR_IIT{2}.F;
        dF = F_full - F_reduced;
        
        IIT_df(t) = dF;
        IIT_prob(t) = 1 / (1 + exp(-dF));
        
        fprintf('  IIT: dF = %.2f, Pp = %.3f\n', dF, IIT_prob(t));
    end
    
    clear GCM PEB
end

% Save results
save(fullfile(output_dir, 'IIT_combined_results.mat'), ...
     'time_window_increments', 'IIT_prob', 'IIT_df');

%% intrinsic connectivity tests 

short_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw'); 
cd(short_dir)

% Run PEB on G-matrix for each time window 
for t = 1:length(time_window_increments)

    current_time = time_window_increments(t);

    gcm_filename = ['GCM_combined_' num2str(current_time) 'ms.mat'];
    load(fullfile(output_dir, gcm_filename))
    n_total = length(GCM);
    
    % Run PEB
    M = struct();
    M.Q = 'all';
    M.X = ones(n_total, 1);
    field = {'G'};

    PEB_G = spm_dcm_peb(GCM, M, field);
    
    peb_g_filename = ['PEB_G_' num2str(current_time) 'ms.mat'];

    save(fullfile(output_dir, peb_g_filename), 'PEB_G');
    
    clear GCM PEB_G
    
end

short_time_window = 100:100:1000;

for t = 1:length(time_window_increments)

    current_time = time_window_increments(t);

    peb_g_filename = fullfile(output_dir, ['PEB_G_' num2str(current_time) 'ms.mat']);

    load(peb_g_filename, 'PEB_G');

    BMA_G = spm_dcm_peb_bmc(PEB_G);

    bma_g_filename = ['BMA_G_' num2str(current_time) 'ms.mat'];
    save(fullfile(output_dir, bma_g_filename), 'BMA_G');

    FG_Pp(1, t) = BMA_G.Pp(3); %left and right FG
    FG_Pp(2, t) = BMA_G.Pp(4);
    clear PEB_G BMA_G
end

% Calculate combined FG measure (average of L_FFA and R_FFA)
FG_intrinsic_prob = mean(FG_Pp([1, 2], :), 1, 'omitnan');

save(fullfile(output_dir, 'G_matrix_results.mat'), ...
     'short_time_window', 'FG_Pp', 'FG_intrinsic_prob');
