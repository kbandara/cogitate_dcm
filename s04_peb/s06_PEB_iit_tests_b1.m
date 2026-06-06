% PEB analysis for IIT - batch 1 
% kav bandara unimelb 2025

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');
spm_jobman('initcfg');
clear matlabbatch;

clear all; clc;

%{

this analysis runs PEBs on the GCMs created in the previous step to test
for IIT related predictions 

%}

% INPUT: Variables and paths setup

base_filepath  = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH1/'; 

% output directory for DCM results
output_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw'); 
cd(output_dir)

base_filepath = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH1/';
spm_path = fullfile(base_filepath, 'derivatives', 'preprocessed_files');

P_Names_struct = load(fullfile(base_filepath, 'participants.mat'));
P_Names_table = P_Names_struct.participants;
p_names = P_Names_table{:, 1}';
p_names = p_names(~strcmp(p_names, 'CA136')); %exclude CA136 for too few hits

indiv_mni_coords = load(fullfile(spm_path, 'dcm_onset', 'indiv_peak_coords.mat'));
indiv_mni_coords = indiv_mni_coords.indiv_mni_coords;

roi_labels = {'L_V1', 'R_V1', 'L_FFA', 'R_FFA', 'L_PFC', 'R_PFC'};

% Define the end points for each sliding window in milliseconds
time_window_increments = 300:100:1000;

%% extrinsic connectivity test

% store results
IIT_prob = nan(1, length(time_window_increments));
IIT_df = nan(1, length(time_window_increments));

% run PEB for each time window

for t = 1:length(time_window_increments)
    
    current_time = time_window_increments(t);
    
    fprintf('\n=== Processing %dms ===\n', current_time);

    gcm_file = fullfile(output_dir, ['GCM_' num2str(current_time) '.mat']);

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

%% REVIEW: PEB results in GUI (interactive step)

clear PEB GCM BMA 
%load GCM & PEB you want to view results for 
load(fullfile(output_dir,'GCM_combined_100ms.mat')); %change these depending on time window
load(fullfile(output_dir,'PEB_combined_100ms.mat')); %change these depending on time window

% Search over nested PEB models.
BMA = spm_dcm_peb_bmc(PEB);

spm_dcm_peb_review(BMA,GCM)

%% intrinsic connectivity tests 

short_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw'); 
cd(short_dir)

short_time_window = time_window_increments;

% Run PEB on G-matrix for each time window 
for t = 1:length(short_time_window)

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

for t = 1:length(short_time_window)

    current_time = short_time_window(t);

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
V1_intrinsic_prob = mean(V1_Pp([1, 2], :), 1, 'omitnan');

save(fullfile(output_dir, 'G_matrix_results.mat'), ...
     'short_time_window', 'FG_Pp', 'FG_intrinsic_prob');
%% SAME FOR LONG TRIALS 

long_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files', 'dcm_etw', 'long_trials'); 
cd(long_dir)

long_time_window = 100:100:1500;

% Extract Pp values
long_G_Pp = nan(n_nodes, length(long_time_window));
long_G_Ep = nan(n_nodes, length(long_time_window));

for t = 1:length(long_time_window)
    
    current_time = long_time_window(t);
    
    peb_g_filename = fullfile(long_dir, ['long_trials_PEB_G_' num2str(current_time) 'ms.mat']);
    
    if ~exist(peb_g_filename, 'file')
        continue;
    end
    
    load(peb_g_filename, 'PEB_G');
    
    BMA_G = spm_dcm_peb_bmc(PEB_G);
    
    bma_g_filename = ['long_trials_BMA_G_' num2str(current_time) 'ms.mat'];
    save(fullfile(long_dir, bma_g_filename), 'BMA_G');
    
    for p = 1:length(BMA_G.Pnames)
        pname = BMA_G.Pnames{p};
        
        for node = 1:n_nodes
            pattern1 = sprintf('G(%d)', node);
            pattern2 = sprintf('G(%d,', node);
            
            if contains(pname, pattern1) || contains(pname, pattern2)
                long_G_Pp(node, t) = BMA_G.Pp(p);
                long_G_Ep(node, t) = BMA_G.Ep(p);
                break;
            end
        end
    end
    
    fprintf('  L_FFA Pp = %.3f, R_FFA Pp = %.3f\n', long_G_Pp(3, t), long_G_Pp(4, t));
    
    clear PEB_G BMA_G
end

long_FG_intrinsic_prob = mean(long_G_Pp([3, 4], :), 1, 'omitnan');

save(fullfile(long_dir, 'G_matrix_results.mat'), ...
     'long_time_window', 'long_G_Pp', 'long_G_Ep', 'long_FG_intrinsic_prob', 'node_names');
%}
