function s01_prepro_b2(p_name) 

%% Pre-Processing of MEG data 

% pre processing batch 2 

% Kav Bandara, University of Melbourne, 2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        
%                            TESTING BLOCK
%                        uncomment to run on slurm
%p_name = 1; 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12')

spm('defaults', 'eeg');

%%%%%%%%%%           SECTION 0 - LOAD FILES & DEFINE VARIABLES          %%%%%%%%%%%%

%% Settings and filenames for pre-processing
epochTimeWindow = [-100 2000]; %epoch around this time window (ms)

base_filepath  = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/'; 
cd(base_filepath)

%Load participant names
P_Names_struct = load('p_names.mat');  
P_Names_table =  P_Names_struct.uniqueToFirstList; %{'CA101'}; 
P_Names_str = P_Names_table;
P_Names = cellstr(P_Names_str);  % 1x52 cell array of char
P_Names = {P_Names{p_name}}; %single p-name for this job
P_Names = {'CB082'}; 

runs = {'01','02','03','04','05'};


%%
%%%%%%%%%%             SECTION I - PREPROCESSING            %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   STEP 1: Convert the datafile to SPM readable format   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cc = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{cc} '/ses-1/meg/'];
    cd(filepath)
    new_folder = 'spm_preprocessed_files';
    if ~exist(fullfile(filepath, new_folder), 'dir')
        mkdir(filepath, new_folder)
    end

    for i = 1:length(runs)

        filename = fullfile(filepath, ['sub-' P_Names{cc} '_ses-1_task-dur_run-' runs{i} '_filt.fif']);

        spm_jobman('initcfg');
        
        % %%Convert (from EDF to dat/mat file)
        matlabbatch{1}.spm.meeg.convert.dataset = {filename};
        matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
        matlabbatch{1}.spm.meeg.convert.channels{1}.all = 'all';
        matlabbatch{1}.spm.meeg.convert.outfile = fullfile(filepath, new_folder, ['spmeeg_sub-' P_Names{cc} '_ses-1_task-dur_run-' runs{i} '_meg']);
        matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
        matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
        matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
        matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
        matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';

        spm_jobman('run',matlabbatch);
               
        clear matlabbatch
        
    end
end


%% Downsample step  

for ds = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{ds} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)

    for i = 1:length(runs)

        filename = fullfile(filepath, ['spmeeg_sub-' P_Names{ds} '_ses-1_task-dur_run-' runs{i} '_meg.mat']);
    
        fprintf('Downsampling Participant %d of %d (Sub: %s), Run %s...\n', ds, length(P_Names), P_Names{ds}, runs{i});

        spm_jobman('initcfg');

        matlabbatch{1}.spm.meeg.preproc.downsample.D = {filename};
        matlabbatch{1}.spm.meeg.preproc.downsample.fsample_new = 200;
        matlabbatch{1}.spm.meeg.preproc.downsample.method = 'downsample';
        matlabbatch{1}.spm.meeg.preproc.downsample.prefix = 'd_';   
            
        spm_jobman('run',matlabbatch);
        
        clear D
       
        clear matlabbatch
    end
end

%% RENAME TRIGGERS TO EASE PREPROCESSING AND SUBSEQENT SOURCE RECON/DCM

%intermediate step where we create a new trigger value based on what is contained in the trial 

for xx = 1:length(P_Names)

    %navigate to this subject folder
    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{xx} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)

    for l = 1:length(runs)

        %navigate to this subjects run
        filename = fullfile(filepath, ['d_spmeeg_sub-' P_Names{xx} '_ses-1_task-dur_run-' runs{l} '_meg.mat']);
        load(filename);
        
        % Get the event structure
        events = D.trials.events;
	
        fprintf('Loaded %s. Processing %d events...\n', filename, numel(events));
               
        for i = 1:numel(events)
            % Get current event
            current_event = events(i);
            
            % Find Trial Onsets using trigger value and type
            is_sti101_up = false; %initialise
            if ischar(current_event.type)
                is_sti101_up = strcmp(current_event.type, 'STI101_up');
            end
        
            is_valid_stimulus = false; %initialise
            if ~isempty(current_event.value)
                val = current_event.value(1);
                is_valid_stimulus = (val >= 1 && val <= 80);
            end

            if is_sti101_up && is_valid_stimulus
                % This is a valid trial onset, let's process it :D
                onset_idx = i;
                onset_time = current_event.time;
                original_onset_value = current_event.value(1);
        
                 % change stimulus onset to triggers 
                transformed_onset_value = 0;
                if original_onset_value >= 1 && original_onset_value <= 20
                    transformed_onset_value = 90;  % faces
                elseif original_onset_value >= 21 && original_onset_value <= 40
                    transformed_onset_value = 91;  % objects
                elseif original_onset_value >= 41 && original_onset_value <= 60
                    transformed_onset_value = 92;  % letters
                elseif original_onset_value >= 61 && original_onset_value <= 80
                    transformed_onset_value = 93;  % falses
                end

                % This array will hold all triggers for concatenation for this trial
                final_triggers = [transformed_onset_value];
                
                % Find Subsequent Triggers when type is STI101_up
                codes_to_find = [152, 153, 201, 202, 203]; %151 excluded so it cant be marked as bad later
                
                for j = (i + 1):numel(events)
                    if events(j).time > (onset_time + 0.4), break; end %break loop if we search outside 100ms from stim onset
                    
                    % Check if the subsequent event meets BOTH strict conditions
                    is_valid_type = strcmp(events(j).type, 'STI101_up');
                    is_valid_value = ismember(events(j).value(1), codes_to_find);
        
                    if is_valid_type && is_valid_value
                        final_triggers = [final_triggers, events(j).value(1)];
                        events(j).value = 0; % Nullify trigger to prevent reuse
                    end
                end
                
                % --- 3. Concatenate and Update the Onset Trigger ---
                if numel(final_triggers) > 1
                    % Convert all numbers to a single concatenated string
                    new_trigger_str = sprintf('%d', final_triggers);
                    % Convert string back to a number
                    new_trigger_val = str2double(new_trigger_str);
                    
                    % Update the original onset event with the new value
                    events(onset_idx).value = new_trigger_val;              
                    fprintf('Updated onset at %.3fs. Old value: %d, New value: %ld\n', onset_time, original_onset_value, new_trigger_val);
                else
                    % If no other triggers were found, just update with 90 or 95
                    events(onset_idx).value = transformed_onset_value;
                end
            end
        end
        
        % Replace the old events structure with our modified one
        D.trials.events = events;
        
        % Save the modified D object
        save(D.fname, 'D');

        fprintf('\nProcessing complete.\n');
        fprintf('Modified data saved to: %s\n', filename);
    end
end

%% Artefact detection -- rejection happens later

for aa = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{aa} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)
%{  
    %Initialize list of files to delete for this participant after this stage
    files_to_delete_after_artefact_detection = {};
%}
    for i = 1:length(runs)

        filename = fullfile(filepath, ['d_spmeeg_sub-' P_Names{aa} '_ses-1_task-dur_run-' runs{i} '_meg.mat']);
        filename_dat = fullfile(filepath, ['d_spmeeg_sub-' P_Names{aa} '_ses-1_task-dur_run-' runs{i} '_meg.dat']);

        spm_jobman('initcfg');
%{
        files_to_delete_after_artefact_detection{end+1} = filename;
        files_to_delete_after_artefact_detection{end+1} = filename_dat;
%}  
        %Setup artefact detection
        matlabbatch{1}.spm.meeg.preproc.artefact.D = {filename};
        matlabbatch{1}.spm.meeg.preproc.artefact.mode = 'mark';
        matlabbatch{1}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
        matlabbatch{1}.spm.meeg.preproc.artefact.append = true;

        matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).channels{1}.type = 'MEGMAG';
        matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).fun.threshchan.threshold = 5000; %5000 fT
        matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).fun.threshchan.excwin = 1000;
        
        matlabbatch{1}.spm.meeg.preproc.artefact.methods(2).channels{1}.type = 'MEGPLANAR';
        matlabbatch{1}.spm.meeg.preproc.artefact.methods(2).fun.threshchan.threshold = 500; %500 ft/mm or 5000ft/cm
        matlabbatch{1}.spm.meeg.preproc.artefact.methods(2).fun.threshchan.excwin = 1000;
        
        matlabbatch{1}.spm.meeg.preproc.artefact.prefix = 'a_';
    
        spm_jobman('run',matlabbatch);
        
        clear D
       
        clear matlabbatch
    end
    
%{
 for f_idx = 1:length(files_to_delete_after_artefact_detection)
        file_to_del = files_to_delete_after_artefact_detection{f_idx};
        if exist(file_to_del, 'file')
            try
                delete(file_to_del);
                fprintf('    Deleted: %s\n', file_to_del);
            catch ME_del
                fprintf('    WARNING: Could not delete %s. Reason: %s\n', file_to_del, ME_del.message);
            end
        else
        end
    end
    clear files_to_delete_after_artefact_detection;
  %}
end

%% Now Epoching!!

for pp = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{pp} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)

    for i = 1:length(runs)
        
        filename = fullfile(filepath, ['a_d_spmeeg_sub-' P_Names{pp} '_ses-1_task-dur_run-' runs{i} '_meg.mat']);
        
        spm_jobman('initcfg');

        matlabbatch{1}.spm.meeg.preproc.epoch.D = {filename};
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.timewin = epochTimeWindow;
        %  FACES = 90
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'face_target_1000ms'; %90 = face; 152 = 1000ms; 201 = target
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 90152201;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).conditionlabel = 'face_relevant_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventvalue = 90152202;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).conditionlabel = 'face_irrelevant_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).eventvalue = 90152203;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(4).conditionlabel = 'face_target_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(4).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(4).eventvalue = 90153201;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(4).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(5).conditionlabel = 'face_relevant_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(5).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(5).eventvalue = 90153202;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(5).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(6).conditionlabel = 'face_irrelevant_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(6).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(6).eventvalue = 90153203;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(6).trlshift = 0;
        %  OBJECTS = 91
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(7).conditionlabel = 'object_target_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(7).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(7).eventvalue = 91152201;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(7).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(8).conditionlabel = 'object_relevant_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(8).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(8).eventvalue = 91152202;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(8).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(9).conditionlabel = 'object_irrelevant_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(9).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(9).eventvalue = 91152203;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(9).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(10).conditionlabel = 'object_target_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(10).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(10).eventvalue = 91153201;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(10).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(11).conditionlabel = 'object_relevant_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(11).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(11).eventvalue = 91153202;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(11).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(12).conditionlabel = 'object_irrelevant_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(12).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(12).eventvalue = 91153203;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(12).trlshift = 0;
        %  LETTERS = 92
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(13).conditionlabel = 'letter_target_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(13).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(13).eventvalue = 92152201;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(13).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(14).conditionlabel = 'letter_relevant_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(14).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(14).eventvalue = 92152202;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(14).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(15).conditionlabel = 'letter_irrelevant_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(15).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(15).eventvalue = 92152203;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(15).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(16).conditionlabel = 'letter_target_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(16).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(16).eventvalue = 92153201;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(16).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(17).conditionlabel = 'letter_relevant_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(17).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(17).eventvalue = 92153202;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(17).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(18).conditionlabel = 'letter_irrelevant_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(18).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(18).eventvalue = 92153203;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(18).trlshift = 0;
        %  FALSES = 93 
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(19).conditionlabel = 'false_target_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(19).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(19).eventvalue = 93152201;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(19).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(20).conditionlabel = 'false_relevant_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(20).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(20).eventvalue = 93152202;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(20).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(21).conditionlabel = 'false_irrelevant_1000ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(21).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(21).eventvalue = 93152203;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(21).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(22).conditionlabel = 'false_target_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(22).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(22).eventvalue = 93153201;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(22).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(23).conditionlabel = 'false_relevant_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(23).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(23).eventvalue = 93153202;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(23).trlshift = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(24).conditionlabel = 'false_irrelevant_1500ms';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(24).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(24).eventvalue = 93153203;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(24).trlshift = 0;
        % --- 500ms (catch, marked bad later) ---
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(25).conditionlabel = '500ms_trials';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(25).eventtype = 'STI101_up';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(25).eventvalue = 151;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(25).trlshift = -0.033; %stimulus onset is 33ms before stimulus presented (or 4 frames)
        matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.prefix = 'e_ALL_';

        spm_jobman('run',matlabbatch);

        load(['e_ALL_a_d_spmeeg_sub-' P_Names{pp} '_ses-1_task-dur_run-' runs{i} '_meg.mat']);
    
        %mark trials with the trigger 151 as bad - these are 500ms trials
        for k = 1:length(D.trials)  
            % Initialize bad as 0 (not bad trial)
            D.trials(k).bad = 0;             
            for l = 1:length(D.trials(k).events)
                % Check if the value is '151'
                if D.trials(k).events(l).value == 151
                    % Mark the trial as bad
                    D.trials(k).bad = 1;
                    break;  % No need to check further if '151' is found
                end
            end
        end
        
        save(D.fname, 'D'); % save D to disk with 151 marked as bad

        %remove these bad trials
        S = [];
        S.D = D.fname; % load 'e_faces_a_d_spmeeg_sub-' P_Names{pp} '_ses-1_task-dur_run-' runs{i} '_meg.mat'
        D = spm_eeg_remove_bad_trials(S);          

        clear D 
       
        clear matlabbatch     
    end
end

% Merge FACE sessions into one file 

for ms = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{ms} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)
%{  
    %Initialize list of files to delete for this participant after this stage <<<
    files_to_delete_after_merge_faces = {};
    for r_idx = 1:length(runs)
        run_file_mat = fullfile(filepath, ['re_ALL_faces_a_d_spmeeg_sub-' P_Names{ms} '_ses-1_task-dur_run-' runs{r_idx} '_meg.mat']);
        run_file_dat = fullfile(filepath, ['re_ALL_faces_a_d_spmeeg_sub-' P_Names{ms} '_ses-1_task-dur_run-' runs{r_idx} '_meg.dat']);
        files_to_delete_after_merge_faces{end+1} = run_file_mat;
        files_to_delete_after_merge_faces{end+1} = run_file_dat;
    end
%}
    spm_jobman('initcfg');
    matlabbatch{1}.spm.meeg.preproc.merge.D = {['re_ALL_a_d_spmeeg_sub-' P_Names{ms} '_ses-1_task-dur_run-01_meg.mat']; ...
                                               ['re_ALL_a_d_spmeeg_sub-' P_Names{ms} '_ses-1_task-dur_run-02_meg.mat']; ...
                                               ['re_ALL_a_d_spmeeg_sub-' P_Names{ms} '_ses-1_task-dur_run-03_meg.mat']; ...
                                               ['re_ALL_a_d_spmeeg_sub-' P_Names{ms} '_ses-1_task-dur_run-04_meg.mat']; ...
                                               ['re_ALL_a_d_spmeeg_sub-' P_Names{ms} '_ses-1_task-dur_run-05_meg.mat']};
    matlabbatch{1}.spm.meeg.preproc.merge.recode.file = '.*';
    matlabbatch{1}.spm.meeg.preproc.merge.recode.labelorg = '.*';
    matlabbatch{1}.spm.meeg.preproc.merge.recode.labelnew = '#labelorg#';
    matlabbatch{1}.spm.meeg.preproc.merge.prefix = 'ms_';

    spm_jobman('run',matlabbatch);

    clear D i      
   
    clear matlabbatch

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Mark missing triggers as bad and sort   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Now reject bad trials

for rr = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{rr} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)

    filename = fullfile(filepath, ['ms_re_ALL_a_d_spmeeg_sub-' P_Names{rr} '_ses-1_task-dur_run-01_meg.mat']);

    spm_jobman('initcfg');

    % Remove bad trials
    matlabbatch{1}.spm.meeg.preproc.artefact.D = {filename};  
    matlabbatch{1}.spm.meeg.preproc.artefact.mode = 'reject';  

    matlabbatch{1}.spm.meeg.preproc.artefact.badchanthresh = 0.2; 
    matlabbatch{1}.spm.meeg.preproc.artefact.append = true;
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).channels{1}.type = 'MEGMAG';
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).fun.events.whatevents.artefacts = 1; % 1 = "all"
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(2).channels{1}.type = 'MEGPLANAR';
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(2).fun.events.whatevents.artefacts = 1; % 1 = "all"

    matlabbatch{1}.spm.meeg.preproc.artefact.prefix = 'r_';

    spm_jobman('run',matlabbatch);
        
    clear D
   
    clear matlabbatch

end 
 

%% robust averaging, low-pass filter and baseline correction 
for ra = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{ra} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)
  
    filename = fullfile(filepath, ['r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{ra} '_ses-1_task-dur_run-01_meg.mat']);

    spm_jobman('initcfg');
 
    %Robustly average 
    matlabbatch{1}.spm.meeg.averaging.average.D = {filename};
    matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.ks = 3;
    matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.bycondition = true;
    matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.savew = false;
    matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.removebad = true; 
    matlabbatch{1}.spm.meeg.averaging.average.plv = false;
    matlabbatch{1}.spm.meeg.averaging.average.prefix = 'ra_';
    
    spm_jobman('run',matlabbatch);
    
   
    clear matlabbatch

    clear D; % clear for next loop
%{  
    % Delete the r_ms_re_ALL_faces_a_... files for this participant
    files_to_delete_after_robuavg_faces = {filename, filename_dat};
    for f_idx = 1:length(files_to_delete_after_robuavg_faces)
        file_to_del = files_to_delete_after_robuavg_faces{f_idx};
        if exist(file_to_del, 'file')
            try
                delete(file_to_del);
                fprintf('    Deleted: %s\n', file_to_del);
            catch ME_del
                fprintf('    WARNING: Could not delete %s. Reason: %s\n', file_to_del, ME_del.message);
            end
        else
            % fprintf('    INFO: File not found for deletion: %s\n', file_to_del);
        end
    end
    clear files_to_delete_after_robuavg_faces;
%}
end

for lf = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{lf} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)

    filename = fullfile(filepath, ['ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{lf} '_ses-1_task-dur_run-01_meg.mat']);
    %{  
    filename_dat = fullfile(filepath, ['ra_r_ms_re_ALL_faces_a_d_spmeeg_sub-' P_Names{lf} '_ses-1_task-dur_run-01_meg.dat']);
%}
    spm_jobman('initcfg');

    %Lowpass filter a second time (robust averaging introduces high frequency noise)
    matlabbatch{1}.spm.meeg.preproc.filter.D = {filename};
    matlabbatch{1}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{1}.spm.meeg.preproc.filter.band = 'low';
    matlabbatch{1}.spm.meeg.preproc.filter.freq = 40; % 40 Hz
    matlabbatch{1}.spm.meeg.preproc.filter.dir = 'twopass'; %default
    matlabbatch{1}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{1}.spm.meeg.preproc.filter.prefix = 'lf_';
    
    spm_jobman('run',matlabbatch);
    
   
    clear matlabbatch
    clear D 
%{  
    % Delete the ra_r_ms_re_ALL_faces_a_... files for this participant
    fprintf('>>> Deleting intermediate files for participant %s after Low-pass Filter (faces) (lf loop):\n', P_Names{lf});
    files_to_delete_after_lpfilter_faces = {filename, filename_dat};
    for f_idx = 1:length(files_to_delete_after_lpfilter_faces)
        file_to_del = files_to_delete_after_lpfilter_faces{f_idx};
        if exist(file_to_del, 'file')
            try
                delete(file_to_del);
                fprintf('    Deleted: %s\n', file_to_del);
            catch ME_del
                fprintf('    WARNING: Could not delete %s. Reason: %s\n', file_to_del, ME_del.message);
            end
        else
            % fprintf('    INFO: File not found for deletion: %s\n', file_to_del);
        end
    end
    clear files_to_delete_after_lpfilter_faces;
%}
end

for bc = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{bc} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)

    filename = fullfile(filepath, ['lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{bc} '_ses-1_task-dur_run-01_meg.mat']);
    %{  
    filename_dat = fullfile(filepath, ['lf_ra_r_ms_re_ALL_faces_a_d_spmeeg_sub-' P_Names{bc} '_ses-1_task-dur_run-01_meg.dat']);
    %}
    spm_jobman('initcfg');

    matlabbatch{1}.spm.meeg.preproc.bc.D =  {filename};
    matlabbatch{1}.spm.meeg.preproc.bc.timewin = [-100 0];
    matlabbatch{1}.spm.meeg.preproc.bc.prefix = 'bc_';
    
    spm_jobman('run',matlabbatch);

    clear D       
   
    clear matlabbatch    
    %{  
    % Delete the lf_ra_r_ms_re_ALL_faces_a_... files for this participant
    fprintf('>>> Deleting intermediate files for participant %s after baseline correction (faces):\n', P_Names{bc});
    files_to_delete_after_bc_faces = {filename, filename_dat};
    for f_idx = 1:length(files_to_delete_after_bc_faces)
        file_to_del = files_to_delete_after_bc_faces{f_idx};
        if exist(file_to_del, 'file')
            try
                delete(file_to_del);
                fprintf('    Deleted: %s\n', file_to_del);
            catch ME_del
                fprintf('    WARNING: Could not delete %s. Reason: %s\n', file_to_del, ME_del.message);
            end
        else
        end
    end
    clear files_to_delete_after_bc_faces;
    %}
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   select MEG channels from final files  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now select MEG channels from preprocessed files 

for select = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{select} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)
    filename = fullfile(filepath, ['bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{select} '_ses-1_task-dur_run-01_meg.mat']);

    spm_jobman('initcfg');
 
    % Define the channel selection task to run on the new file
    matlabbatch{1}.spm.meeg.preproc.crop.D = {filename};
    matlabbatch{1}.spm.meeg.preproc.crop.timewin = [-Inf Inf]; %select all
    matlabbatch{1}.spm.meeg.preproc.crop.freqwin = [-Inf Inf]; %select all 
    matlabbatch{1}.spm.meeg.preproc.crop.channels{1}.type = 'MEGPLANAR';
    matlabbatch{1}.spm.meeg.preproc.crop.channels{2}.type = 'MEGMAG';
    matlabbatch{1}.spm.meeg.preproc.crop.prefix = 'only_meg_';
   
    spm_jobman('run',matlabbatch);
 
    clear D matlabbatch; 

end

valid_conditions = {
    'face_relevant_1000ms',   'face_irrelevant_1000ms', ...
    'face_relevant_1500ms',   'face_irrelevant_1500ms', ...
    'object_relevant_1000ms', 'object_irrelevant_1000ms', ...
    'object_relevant_1500ms', 'object_irrelevant_1500ms', ...
    'letter_relevant_1000ms', 'letter_irrelevant_1000ms', ...
    'letter_relevant_1500ms', 'letter_irrelevant_1500ms', ...
    'false_relevant_1000ms',  'false_irrelevant_1000ms', ...
    'false_relevant_1500ms',  'false_irrelevant_1500ms'
};


%remove target trials 
for notarget = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{notarget} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)
    filename = fullfile(filepath, ['only_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{notarget} '_ses-1_task-dur_run-01_meg.mat']);

    load(filename);
      
    marked_count = 0;
    % Loop through each trial in the dataset
    for i = 1:length(D.trials)
        trial_label = D.trials(i).label;
        
        % ismember is an efficient way to check if the trial's label exists in your list of valid conditions.
        is_valid = ismember(trial_label, valid_conditions);
        
        % If the trial's label is NOT in the valid list, mark the trial as bad.
        if ~is_valid
            D.trials(i).bad = 1; % Mark the trial as bad
            marked_count = marked_count + 1;
        end
    end
    
    save(D.fname, 'D'); % save D to disk with trials marked as bad

    % and remove these bad trials
    S = [];
    S.D = D.fname; 
    D = spm_eeg_remove_bad_trials(S); 

    clear D S;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%     Move files    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update filepath - this is where the preprocessed files are now
filepath = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessed_files/';
cd(filepath)
% Define the new centralized output directory
central_output_dir = fullfile(base_filepath, 'derivatives', 'preprocessed_files');

% Create the centralized output directory if it doesn't exist
if ~exist(central_output_dir, 'dir')
    fprintf('Creating directory: %s\n', central_output_dir);
    mkdir(central_output_dir);
else
    fprintf('Directory already exists: %s\n', central_output_dir);
end

for move = 1:length(P_Names)

    filepath = ['/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessing/sub-' P_Names{move} '/ses-1/meg/spm_preprocessed_files/'];
    cd(filepath)

    %set up variables
    merged_mat = fullfile(filepath, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{move} '_ses-1_task-dur_run-01_meg.mat']);
    merged_dat = fullfile(filepath, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{move} '_ses-1_task-dur_run-01_meg.dat']); 
    target_mat = fullfile(central_output_dir, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{move} '_ses-1_task-dur_run-01_meg.mat']);
    target_dat = fullfile(central_output_dir, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{move} '_ses-1_task-dur_run-01_meg.dat']); 
    
    %move the files
    movefile(merged_mat, target_mat, 'f'); % 'f' to force overwrite if target exists
    movefile(merged_dat, target_dat, 'f');
    
    Dmerged = spm_eeg_load(target_mat); % Load from new location
    Dmerged = Dmerged.fname(target_mat); % Update internal filename and path
    Dmerged.save(); % Save the changes
    
    clear D Dmerged merged_mat merged_dat target_mat target_dat; 
    
end

%% update D.condlist to match D.trials 


for conds = 1:length(P_Names)

    filepath = '/data/gpfs/projects/punim2118/Cogitate_DCM/Cogitate_Data/COG_MEEG_EXP1_BIDS_BATCH2/derivatives/preprocessed_files/';
    cd(filepath)
    filename = fullfile(filepath, ['ronly_meg_bc_lf_ra_r_ms_re_ALL_a_d_spmeeg_sub-' P_Names{conds} '_ses-1_task-dur_run-01_meg.mat']);

    load(filename);
    good_condlist = valid_conditions;
    D.condlist = good_condlist;

    save(D.fname, 'D'); % save D to disk

    clear D;

end
end