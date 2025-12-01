%% ===============================================================
%  PHASE 1 : ERP COMPUTATION PER SUBJECT & TASK GROUP
%  (Updated mapping for T0, T1, T2 meaning)
%  ===============================================================
clear; clc;

% ---------- CONFIGURATION ----------
processed_root = './processed/';
results_folder = './results/erp/';
if ~exist(results_folder, 'dir'), mkdir(results_folder); end

% -------------------- RUN GROUP DEFINITIONS --------------------
% Each group corresponds to a specific experimental condition.
% Multiple runs in the same group will be pooled together.
run_groups = struct(...
    'Baseline_EyesOpen', {'R01'}, ...
    'Baseline_EyesClosed', {'R02'}, ...
    'Real_LeftRight_Fist', {'R03','R07','R11'}, ...
    'Imagined_LeftRight_Fist', {'R04','R08','R12'}, ...
    'Real_Both_FistFeet', {'R05','R09','R13'}, ...
    'Imagined_Both_FistFeet', {'R06','R10','R14'} ...
);

% Labels describing what T0/T1/T2 correspond to in the experiment
event_labels = struct(...
    'T0', 'Rest (no movement)', ...
    'T1', 'Onset of left hand / both fists movement (real or imagined)', ...
    'T2', 'Onset of right hand / both feet movement (real or imagined)' ...
);

group_names = fieldnames(run_groups);

% ---------- FIND SUBJECTS ----------
subjects = dir(fullfile(processed_root, 'S*'));
subjects = {subjects([subjects.isdir]).name};
fprintf('Found %d subjects to process.\n\n', length(subjects));

% ---------- INITIALIZE EEGLAB ----------
[ALLEEG, EEG, CURRENTSET] = eeglab('nogui');

% ---------- MAIN LOOP OVER SUBJECTS ----------
for s = 1:length(subjects)
    subject = subjects{s};
    fprintf('\nProcessing subject %s ...\n', subject);
    subj_folder = fullfile(processed_root, subject);

    ERP_struct = struct();

    % ----- LOOP OVER TASK GROUPS -----
    for g = 1:length(group_names)
        grp_name = group_names{g};
        runs = run_groups.(grp_name);

        % Ensure runs is always a cell array
        if ischar(runs)
            runs = {runs};
        end

        % Matrices to collect all epochs across runs for this group
        all_T0 = [];
        all_T1 = [];
        all_T2 = [];

        % ---- LOAD AND COLLECT EPOCHS ----
        for r = 1:length(runs)
            run = runs{r};

            % Build expected filename for processed .set file
            setfile = fullfile(subj_folder, sprintf('%s_%s_processed.set', subject, run));
            if ~isfile(setfile)
                fprintf('Missing %s\n', setfile);
                continue;
            end

            % Load epoched, cleaned EEG
            EEG = pop_loadset('filename', sprintf('%s_%s_processed.set', subject, run), ...
                              'filepath', subj_folder);

            % Loop through trials and sort them by event type
            for i = 1:EEG.trials
                etype = EEG.epoch(i).eventtype;

                % eventtype-nested: convert to string
                if iscell(etype), etype = etype{1}; end

                % Append epoch data to corresponding event accumulator
                switch etype
                    case 'T0'
                        all_T0 = cat(3, all_T0, EEG.data(:,:,i));   % concatenate along 3rd dim
                    case 'T1'
                        all_T1 = cat(3, all_T1, EEG.data(:,:,i));
                    case 'T2'
                        all_T2 = cat(3, all_T2, EEG.data(:,:,i));
                end
            end
        end

        % ---- COMPUTE AVERAGE ERPs ----
        % Skip group if absolutely no epochs found
        if isempty(all_T0) && isempty(all_T1) && isempty(all_T2)
            fprintf('No epochs found for %s %s\n', subject, grp_name);
            continue;
        end

        % Compute trial-averaged ERP per event type
        ERP_struct.(grp_name).T0_mean = mean(all_T0, 3, 'omitnan');
        ERP_struct.(grp_name).T1_mean = mean(all_T1, 3, 'omitnan');
        ERP_struct.(grp_name).T2_mean = mean(all_T2, 3, 'omitnan');

        % Count epochs in each category
        ERP_struct.(grp_name).T0_n = size(all_T0,3);
        ERP_struct.(grp_name).T1_n = size(all_T1,3);
        ERP_struct.(grp_name).T2_n = size(all_T2,3);

        % Store human-readable event descriptions
        ERP_struct.(grp_name).T0_info = event_labels.T0;
        ERP_struct.(grp_name).T1_info = event_labels.T1;
        ERP_struct.(grp_name).T2_info = event_labels.T2;

        fprintf('%s: %d T0 (Rest), %d T1 (Left/Both Fist), %d T2 (Right/Both Feet)\n', ...
                grp_name, ERP_struct.(grp_name).T0_n, ERP_struct.(grp_name).T1_n, ERP_struct.(grp_name).T2_n);
    end

    % ---- SAVE PER-SUBJECT ERP FILE ----
    save(fullfile(results_folder, sprintf('%s_ERP.mat', subject)), 'ERP_struct', '-v7.3');
    fprintf('Saved ERP data for %s\n', subject);
end

fprintf('\nAll subjects processed. ERP files saved in %s\n', results_folder);

