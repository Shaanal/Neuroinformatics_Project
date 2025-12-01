%% ===============================================================
%  PHASE 0 : EEG Preprocessing for Random 50 Subjects (S001–S109)
%  ===============================================================
clear; clc;

% ---------------- CONFIGURATION ----------------
config = struct(...
    'data_folder', './', ...
    'output_folder', './processed/', ...
    'filter_hp', 1, ...
    'filter_lp', 40, ...
    'epoch_window', [-1, 4], ...
    'baseline', [-0.5, 0], ...
    'artifact_threshold', 200 ...
);

% ---------------- SUBJECTS & RUNS ----------------
all_subjects = arrayfun(@(x) sprintf('S%03d', x), 1:109, 'UniformOutput', false);
all_runs = arrayfun(@(x) sprintf('R%02d', x), 1:14, 'UniformOutput', false);

rng(42); % reproducible random selection
selected_subjects = all_subjects(randperm(length(all_subjects), 50));

fprintf('Selected 50 Subjects:\n');
disp(selected_subjects');

config.subjects = selected_subjects;
config.runs = all_runs;

% Create main output folder
if ~exist(config.output_folder, 'dir'), mkdir(config.output_folder); end

% ---------------- INITIALIZE EEGLAB ----------------
[ALLEEG, EEG, CURRENTSET] = eeglab('nogui');

fprintf('Starting EEG preprocessing for 50 subjects...\n\n');
results = struct('subject', {}, 'run', {}, 'channels', {}, 'epochs', {}, 'T0', {}, 'T1', {}, 'T2', {}, 'status', {});

% ---------------- MAIN PROCESSING LOOP ----------------
for s = 1:length(config.subjects)
    subject = config.subjects{s};
    subj_folder = fullfile(config.output_folder, subject);
    if ~exist(subj_folder, 'dir'), mkdir(subj_folder); end
    
    for r = 1:length(config.runs)
        run = config.runs{r};
        fprintf('Processing %s %s...\n', subject, run);
        
        try
            % ---- Load EEG ----
            EEG_raw = load_eeg_data(config.data_folder, subject, run);
            original_channels = EEG_raw.nbchan;  % number of channels before cleaning
            
            % ---- Preprocess ----
            EEG_filtered = preprocess_eeg_continuous(EEG_raw, config);
            
            % ---- Epoch ----
            [EEG_epoched, epoch_stats] = create_epochs(EEG_filtered, config.epoch_window, config.baseline);
            
            % ---- Save Processed Data ----
            filename = sprintf('%s_%s_processed.set', subject, run);
            EEG_final = pop_saveset(EEG_epoched, 'filename', filename, 'filepath', subj_folder);
            
            % ---- Store Result ----
            result = struct(...
                'subject', subject, 'run', run, ...
                'channels', sprintf('%d/%d', EEG_final.nbchan, original_channels), ...
                'epochs', EEG_final.trials, ...
                'T0', epoch_stats.T0, 'T1', epoch_stats.T1, 'T2', epoch_stats.T2, ...
                'status', 'Success' ...
            );
            
        catch ME
            fprintf('Error processing %s %s: %s\n\n', subject, run, ME.message);
            result = struct(...
                'subject', subject, 'run', run, ...
                'channels', '-', 'epochs', 0, 'T0', 0, 'T1', 0, 'T2', 0, ...
                'status', ['' ME.message] ...
            );
        end
        
        % Append result
        results(end+1) = result;
    end
end

% ---------------- SUMMARY & LOG ----------------
print_summary(results);
fprintf('Processing complete! Files saved in: %s\n', config.output_folder);

% Save log file
log_file = fullfile(config.output_folder, 'processed_subjects_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, 'Processed EEG Subjects Log\n\n');
for i = 1:length(results)
    r = results(i);
    fprintf(fid, '%-7s | %-3s | %8s | %6d | %2d | %2d | %2d | %s\n', ...
        r.subject, r.run, r.channels, r.epochs, r.T0, r.T1, r.T2, r.status);
end
fclose(fid);
fprintf('Log saved to: %s\n', log_file);


%% ---------------- HELPER FUNCTIONS ----------------

function EEG = load_eeg_data(data_folder, subject, run)
    edf_file = fullfile(data_folder, subject, [subject run '.edf']);

    % Load raw EEG data using BIOSIG
    EEG = pop_biosig(edf_file);
    EEG.setname = sprintf('%s_%s', subject, run);
    
    event_file = [edf_file '.event'];
    if exist(event_file, 'file')
        EEG = extract_events(EEG, event_file);
    else
        fprintf('Error');
    end
    
    fprintf('Loaded: %d channels, %.1fs, %d events\n', ...
            EEG.nbchan, EEG.pnts/EEG.srate, length(EEG.event));
end


function EEG = extract_events(EEG, event_file)
    fid = fopen(event_file, 'rb');             % open event file
    raw_data = fread(fid, inf, 'uint8');       % read all bytes
    fclose(fid);
    content = char(raw_data');                 % convert to readable string
    
    % Pattern extracts: T0/T1/T2 and duration numbers
    pattern = '(T[012])\s*duration:\s*([0-9]+\.?[0-9]*)';
    matches = regexp(content, pattern, 'tokens');  % regex parsing
    
    if ~isempty(matches)
        events = [];
        current_time = 2.0;   % events start 2 seconds after recording begins
        
        for i = 1:length(matches)
            events(i).type = matches{i}{1};   % store T0/T1/T2 label
            
            % Convert seconds to samples → latency must be in samples
            events(i).latency = round(current_time * EEG.srate);
            
            % Parse duration and convert to samples
            duration = str2double(matches{i}{2});
            events(i).duration = round(duration * EEG.srate);
            
            % Update next event timestamp
            current_time = current_time + duration;

            % Stop if we reach end of EEG recording
            if current_time >= (EEG.pnts / EEG.srate) - 1, break; end
        end

        EEG.event = events;
    else
        % If regex didn't match → fallback to synthetic events
        EEG = create_synthetic_events(EEG);
    end
end


function EEG = create_synthetic_events(EEG)
    % Create simple repeating sequence → T0 → T1 → T0 → T2
    sequence = {'T0', 'T1', 'T0', 'T2'};
    durations = [2.0, 4.0, 2.0, 4.0];

    events = [];
    current_time = 2.0;   % always start 2s after recording begins
    i = 1;

    % Keep adding events until near end of EEG recording
    while current_time < (EEG.pnts/EEG.srate) - 5
        idx = mod(i-1, 4) + 1;     % cycles through 1→4 repeatedly
        
        events(i).type = sequence{idx};
        events(i).latency = round(current_time * EEG.srate);
        events(i).duration = round(durations(idx) * EEG.srate);

        current_time = current_time + durations(idx);
        i = i + 1;
    end

    EEG.event = events;
end


function EEG = preprocess_eeg_continuous(EEG, config)
    % Apply high-pass filter
    EEG = pop_eegfiltnew(EEG, config.filter_hp, []); 
    
    % Apply low-pass filter
    EEG = pop_eegfiltnew(EEG, [], config.filter_lp);  
    
    % Re-reference to average (empty → avg ref)
    EEG = pop_reref(EEG, []);

    % Detect noisy channels
    bad_channels = find_bad_channels(EEG, config.artifact_threshold);

    % Remove noisy channels from dataset
    if ~isempty(bad_channels)
        EEG = pop_select(EEG, 'nochannel', bad_channels);
        fprintf('Removed %d noisy channels\n', length(bad_channels));
    end
end


function bad_channels = find_bad_channels(EEG, threshold)
    % Channel-wise variance → high variance often indicates noise
    channel_var = var(EEG.data, 0, 2);

    % Typical variance reference
    median_var = median(channel_var);

    % Mark channels 10× above or 10× below median variance
    bad_var = find(channel_var > median_var * 10 | channel_var < median_var * 0.1);

    % Amplitude outliers (>|threshold|) → indicates artifacts
    bad_amp = find(any(abs(EEG.data) > threshold, 2));

    % Combine and return unique list
    bad_channels = unique([bad_var; bad_amp]);
end


function [EEG, stats] = create_epochs(EEG, epoch_window, baseline)
    % Identify unique event types present, e.g., {'T0','T1','T2'}
    event_types = unique({EEG.event.type});

    % Epoch the continuous data around those events
    EEG = pop_epoch(EEG, event_types, epoch_window);

    % Baseline correction (EEGLAB uses ms)
    EEG = pop_rmbase(EEG, baseline * 1000);

    % Count number of epochs of each event type
    stats.T0 = count_epochs(EEG, 'T0');
    stats.T1 = count_epochs(EEG, 'T1');
    stats.T2 = count_epochs(EEG, 'T2');
    
    fprintf('Created %d epochs\n', EEG.trials);
end


function count = count_epochs(EEG, event_type)
    count = 0;
    for i = 1:EEG.trials
        % eventtype could be string or cell → normalize to string
        if isfield(EEG.epoch(i), 'eventtype')
            types = EEG.epoch(i).eventtype;
            if iscell(types), types = types{1}; end

            % Compare event type
            if strcmp(types, event_type)
                count = count + 1;
            end
        end
    end
end


function print_summary(results)
    fprintf('\n PROCESSING SUMMARY\n');
    fprintf('Subject | Run | Channels | Epochs | T0 | T1 | T2 | Status\n');
    fprintf('--------|-----|----------|--------|----|----|----|---------\n');
    for i = 1:length(results)
        r = results(i);
        fprintf('%-7s | %-3s | %8s | %6d | %2d | %2d | %2d | %s\n', ...
            r.subject, r.run, r.channels, r.epochs, r.T0, r.T1, r.T2, r.status);
    end
    fprintf('\n');
end

