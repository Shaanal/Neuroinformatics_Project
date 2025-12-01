%% ===============================================================
%  PHASE 0 : VERIFY PROCESSED EEG FILES & CREATE SUMMARY METADATA
%  ===============================================================
clear; clc;

% ---------- CONFIGURATION ----------
processed_root = './processed/';
output_file = './results/per_subject_summary.csv';

% Create results folder if it doesn't exist
if ~exist('./results', 'dir'), mkdir('./results'); end

% ---------- FIND ALL PROCESSED SET FILES ----------
% recursively search inside subject folders for files ending in *_processed.set
files = dir(fullfile(processed_root, 'S*', '*_processed.set'));
fprintf('Found %d processed EEG files.\n\n', length(files));

% ---------- PREPARE SUMMARY TABLE ----------
% Pre-allocate table with fixed variable types for consistent output
summary = table('Size', [0 10], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double','double','string'}, ...
    'VariableNames', {'Subject','Run','SRate','OrigCh','KeptCh','Trials','T0','T1','T2','Status'});

% ---------- INITIALIZE EEGLAB ----------
[ALLEEG, EEG, CURRENTSET] = eeglab('nogui');

% ---------- LOOP THROUGH FILES ----------
for i = 1:length(files)
    try
        filepath = fullfile(files(i).folder, files(i).name);
        fprintf('Checking %s ...\n', filepath);

        % Load the .set file for inspection
        EEG = pop_loadset('filename', files(i).name, 'filepath', files(i).folder);

        % Extract subject name from folder path  
        % (works for both processed/Sxxx and processed/Sxxx/)
        subj = extractBetween(files(i).folder, 'processed/', '/');
        if isempty(subj)
            subj = extractAfter(files(i).folder, 'processed/');
        end
        subject = string(subj);

        % Extract run label (e.g., R001) using regex on filename
        runinfo = regexp(files(i).name, '(R\d+)', 'tokens', 'once');
        run = string(runinfo{1});

        % collect metadata
        srate = EEG.srate;          % sampling rate
        keptCh = EEG.nbchan;        % final number of channels
        origCh = EEG.nbchan;        % assuming original channels unknown here
        trials = EEG.trials;        % number of epochs in data

        % ---------- count events per type ----------
        % EEG.epoch.eventtype may be nested cells → flatten if needed
        types = {EEG.epoch.eventtype};
        if iscell(types{1})
            % extract the first (main) event type from nested cell
            types = cellfun(@(x) x{1}, types, 'UniformOutput', false);
        end

        % Count occurrences per condition
        T0 = sum(strcmp(types, 'T0'));
        T1 = sum(strcmp(types, 'T1'));
        T2 = sum(strcmp(types, 'T2'));

        % Append new row into the summary table
        summary = [summary; {subject, run, srate, origCh, keptCh, trials, T0, T1, T2, "OK"}];

    catch ME
        warning('Could not read %s: %s', files(i).name, ME.message);

        % Attempt subject & run extraction even when load fails
        subj = extractBetween(files(i).folder, 'processed/', '/');
        if isempty(subj)
            subj = extractAfter(files(i).folder, 'processed/');
        end
        subject = string(subj);

        runinfo = regexp(files(i).name, '(R\d+)', 'tokens', 'once');
        if ~isempty(runinfo)
            run = string(runinfo{1});
        else
            run = "";
        end

        % Insert failure row with NaNs
        summary = [summary; {subject, run, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "LoadFailed"}];
    end
end

% ---------- SAVE SUMMARY ----------
writetable(summary, output_file);
fprintf('\nSummary saved to %s\n', output_file);
disp(head(summary));
fprintf('\nQuality overview:\n');
fprintf('  Total files checked : %d\n', height(summary));
fprintf('  Successful loads    : %d\n', sum(summary.Status == "OK"));
fprintf('  Failed loads        : %d\n', sum(summary.Status ~= "OK"));
fprintf('  Mean #trials        : %.1f\n', mean(summary.Trials,'omitnan'));
fprintf('  Mean #channels      : %.1f\n', mean(summary.KeptCh,'omitnan'));
fprintf('\nDone.\n');

