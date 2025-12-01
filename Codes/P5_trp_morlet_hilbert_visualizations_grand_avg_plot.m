%% ======================================================================================
%  PHASE 5 : Time Frequency Analysis: Morlet + Hilbert, per-subject plots, Executed vs Imagined, Grand Averages
% ====================================================================================== 
clear; clc;

%% ---------------- USER / PATH CONFIG ----------------
processed_root = './processed/';
results_root   = './results/';
subject_outdir = fullfile(results_root, 'tfr_subjects');
hilbert_outdir = fullfile(results_root, 'hilbert_subjects');
grand_avg_dir  = fullfile(results_root, 'grand_averages');

% Create output folders if missing
if ~exist(results_root,'dir'), mkdir(results_root); end
if ~exist(subject_outdir,'dir'), mkdir(subject_outdir); end
if ~exist(hilbert_outdir,'dir'), mkdir(hilbert_outdir); end
if ~exist(grand_avg_dir,'dir'), mkdir(grand_avg_dir); end

maxSubjects = 250;   % upper bound for processing

%% ---------------- TIME-FREQUENCY CONFIG ----------------
freqs = 4:1:40;              % frequencies for both Morlet & Hilbert TFRs
nFreqs = numel(freqs);
ncycles_min = 3;             % minimum cycles at low frequencies
ncycles_max = 8;             % maximum cycles at high frequencies

baseline_sec = [-0.5 0];     % baseline for ERSP normalization
motor_chans = {'C3','Cz','C4'}; 
plot_time_limits = [];       % empty = use full epoch

% Run groups
executed_runs = [3,7,11];
imagined_runs = [4,8,12];

DO_WAVELET = true; 
DO_HILBERT = true;
save_figs  = true;
verbose    = true;

%% ---------------- FIND FILES ----------------
files = dir(fullfile(processed_root, 'S*', '*_processed.set'));
if isempty(files), error('No processed .set files found under %s', processed_root); end
if numel(files) > maxSubjects, files = files(1:maxSubjects); end
nFiles = numel(files);
if verbose, fprintf('Found %d processed files (limited to %d).\n', nFiles, maxSubjects); end

%% ---------------- EEGLAB INIT ----------------
[ALLEEG, EEG, CURRENTSET] = eeglab('nogui');

%% ---------------- GRAND STRUCT INIT ----------------
% grand.executed.C3 / grand.imagined.C3 etc.
conds = {'executed','imagined'};
for cc = 1:numel(conds)
    for mc = 1:numel(motor_chans)
        grand.(conds{cc}).(motor_chans{mc}).power = [];         % Morlet ERSP
        grand.(conds{cc}).(motor_chans{mc}).itpc  = [];         % Morlet ITPC
        grand.(conds{cc}).(motor_chans{mc}).power_hilbert = []; % Hilbert power
        grand.(conds{cc}).(motor_chans{mc}).nsubj = 0;
        grand.(conds{cc}).(motor_chans{mc}).times_ref = [];
    end
end

%% =============== PROCESS EACH SUBJECT FILE (MAIN LOOP) ===============
for fi = 1:nFiles

    filepath = fullfile(files(fi).folder, files(fi).name);
    if verbose, fprintf('\n--- Processing: %s ---\n', filepath); end

    try
        EEG = pop_loadset('filename', files(fi).name, 'filepath', files(fi).folder);
    catch ME
        warning('Could not load %s: %s', files(fi).name, ME.message);
        continue;
    end

    % Ensure data is epoched (chan × time × trials)
    if ndims(EEG.data) < 3
        warning('File %s not epoched. Skipping.', files(fi).name);
        continue;
    end

    %% ------- Determine executed vs imagined based on run number -------
    [~, fname] = fileparts(files(fi).name); 
    run_match = regexp(fname, 'R(\d+)', 'tokens', 'once');
    if isempty(run_match)
        warning('No run number found in %s.', files(fi).name);
        continue;
    end

    run_num = str2double(run_match{1});
    if ismember(run_num, executed_runs)
        condition = 'executed';
    elseif ismember(run_num, imagined_runs)
        condition = 'imagined';
    else
        if verbose, fprintf('Run %d not relevant. Skipping.\n', run_num); end
        continue;
    end

    %% ------- Basic metadata -------
    chan_labels = {EEG.chanlocs.labels};
    srate = EEG.srate;
    times_sec = EEG.times / 1000;
    nt = numel(times_sec);
    ntrials = EEG.trials;

    %% ------- Baseline index lookup -------
    bidx = find(times_sec >= baseline_sec(1) & times_sec <= baseline_sec(2));
    if isempty(bidx)
        bidx = 1:round(nt * 0.1); % fallback if baseline outside epoch
    end

    %% ------- Time window for plotting -------
    if isempty(plot_time_limits)
        pidx = 1:nt;
    else
        pidx = find(times_sec >= plot_time_limits(1) & times_sec <= plot_time_limits(2));
        if isempty(pidx), pidx = 1:nt; end
    end

    subj_name = regexprep(fname, '_processed$', '', 'ignorecase');

    %% ------- Auto-detect motor channels -------
    chan_idx = nan(1, numel(motor_chans));
    matched_labels = cell(1, numel(motor_chans));
    for c = 1:numel(motor_chans)
        tok = motor_chans{c};
        idx = find(contains(chan_labels, tok, 'IgnoreCase', true));

        % fallback: match beginnings
        if isempty(idx)
            idx = find(~cellfun(@isempty, regexpi(chan_labels, ['^' tok])));
        end

        if ~isempty(idx)
            chan_idx(c) = idx(1);
            matched_labels{c} = chan_labels{idx(1)};
        end
    end

    if all(isnan(chan_idx))
        warning('No motor channels found in %s. Skipping.', subj_name);
        continue;
    end

    valid_chans = find(~isnan(chan_idx));

    %% ===================================================================
    %                     MORLET WAVELET TFR
    % ===================================================================
    subj_tfr = struct();

    if DO_WAVELET
        if verbose
            fprintf('  Morlet TFR: %d channels × %d freqs × %d trials\n', ...
                numel(valid_chans), nFreqs, ntrials);
        end

        for v = 1:numel(valid_chans)

            ch = chan_idx(valid_chans(v));
            data = double(squeeze(EEG.data(ch,:,:))'); % time × trials

            % Allocate power & phase across trials
            pow_all   = zeros(nFreqs, nt, ntrials);
            phase_all = zeros(nFreqs, nt, ntrials);

            for fiq = 1:nFreqs
                f = freqs(fiq);

                % Linearly increasing cycles for higher frequencies
                ncycles = ncycles_min + (ncycles_max - ncycles_min) * ...
                        ((f - freqs(1)) / (freqs(end) - freqs(1)));

                % SD of wavelet Gaussian
                s = ncycles / (2*pi*f);

                % Construct Morlet wavelet on ±3.5σ window
                wavtime = -3.5*s : 1/srate : 3.5*s;
                mwave = exp(2i*pi*f.*wavtime) .* exp(-wavtime.^2./(2*s^2));

                % FFT-based convolution length
                nWave = numel(mwave);
                nConv = nWave + nt - 1;

                waveX = fft(mwave, nConv);  % freq-domain wavelet

                for tr = 1:ntrials
                    sigX = fft(data(:,tr)', nConv);  % trial FFT
                    conv_res = ifft(sigX .* waveX);  % convolution
                    conv_res = conv_res(floor(nWave/2)+1 : floor(nWave/2)+nt);

                    pow_all(fiq,:,tr)   = abs(conv_res).^2;
                    phase_all(fiq,:,tr) = angle(conv_res);
                end
            end

            % Average across trials → freq × time
            mean_pow = mean(pow_all, 3);

            % Baseline mean per frequency
            mean_base = mean(mean_pow(:, bidx), 2);
            mean_base(mean_base==0) = eps;

            % ERSP in dB
            ersp_db = 10 .* log10(bsxfun(@rdivide, mean_pow, mean_base));

            % ITPC = |mean(exp(iθ))|
            itpc = abs(mean(exp(1i * phase_all), 3));

            token_name = motor_chans{valid_chans(v)};
            safe_field = matlab.lang.makeValidName(token_name);

            % Store per-subject
            subj_tfr.(safe_field).freqs = freqs;
            subj_tfr.(safe_field).times = times_sec;
            subj_tfr.(safe_field).power_db = ersp_db;
            subj_tfr.(safe_field).itpc = itpc;

            % Append to grand
            G = grand.(condition).(token_name);
            if isempty(G.times_ref), G.times_ref = times_sec; end
            if isempty(G.power)
                G.power = ersp_db(:, pidx);
            else
                G.power(:,:,end+1) = ersp_db(:, pidx);
            end
            if isempty(G.itpc)
                G.itpc = itpc(:, pidx);
            else
                G.itpc(:,:,end+1) = itpc(:, pidx);
            end
            G.nsubj = G.nsubj + 1;
            grand.(condition).(token_name) = G;
        end

        % Save subject-level wavelet .mat
        save(fullfile(subject_outdir, sprintf('%s_%s_wavelet.mat', subj_name, condition)), ...
             'subj_tfr', 'subj_name', 'chan_idx', 'motor_chans', '-v7.3');
    end

    %% ===================================================================
    %                               HILBERT
    % ===================================================================
    subj_hilb = struct();

    if DO_HILBERT
        if verbose, fprintf('  Hilbert: %d channels\n', numel(valid_chans)); end

        for v = 1:numel(valid_chans)

            ch = chan_idx(valid_chans(v));
            data = double(squeeze(EEG.data(ch,:,:))); % samples × trials

            hilb_pow = zeros(nFreqs, nt);  % freq × time

            %% Hilbert per frequency band
            for fiq = 1:nFreqs
                f = freqs(fiq);

                % Bandwidth = 30% of frequency
                bw = max(2, round(f * 0.3));

                bpFilt = designfilt('bandpassfir','FilterOrder',64, ...
                    'CutoffFrequency1', max(0.5, f - bw/2), ...
                    'CutoffFrequency2', f + bw/2, ...
                    'SampleRate', srate);

                pow_trials = zeros(nt, ntrials);

                for tr = 1:ntrials
                    x = filtfilt(bpFilt, data(:,tr));   % zero-phase bandpass
                    a = abs(hilbert(x));                % analytic amplitude
                    pow_trials(:,tr) = a.^2;            % power per trial
                end

                hilb_pow(fiq,:) = mean(pow_trials,2)';  % avg across trials
            end

            % Baseline-correct in dB
            mean_base = mean(hilb_pow(:, bidx), 2);
            mean_base(mean_base==0) = eps;
            hilb_db = 10 .* log10(bsxfun(@rdivide, hilb_pow, mean_base));

            token_name = motor_chans{valid_chans(v)};
            safe_field = matlab.lang.makeValidName(token_name);

            subj_hilb.(safe_field).freqs = freqs;
            subj_hilb.(safe_field).times = times_sec;
            subj_hilb.(safe_field).power_db = hilb_db;

            % Append to grand Hilbert matrix
            G = grand.(condition).(token_name);
            if isempty(G.times_ref), G.times_ref = times_sec; end
            if isempty(G.power_hilbert)
                G.power_hilbert = hilb_db(:, pidx);
            else
                G.power_hilbert(:,:,end+1) = hilb_db(:, pidx);
            end
            grand.(condition).(token_name) = G;
        end

        % Save subject-level Hilbert .mat
        save(fullfile(hilbert_outdir, sprintf('%s_%s_hilbert.mat', subj_name, condition)), ...
             'subj_hilb', 'subj_name', 'chan_idx', 'motor_chans', '-v7.3');
    end

end % file loop

%% ===================================================================
%                  GRAND AVERAGE & DIFF PLOTS
% ===================================================================
fprintf('\n--- Computing grand-averages and comparisons ---\n');

for mc = 1:numel(motor_chans)
    token = motor_chans{mc};
    Gexec = grand.executed.(token);
    Gimg  = grand.imagined.(token);

    % ---- Morlet Grand Averages ----
    if ~isempty(Gexec.power)
        mean_exec = mean(Gexec.power, 3, 'omitnan');
        plot_save_TFR(mean_exec, freqs, Gexec.times_ref, sprintf('grand_wavelet_executed_%s', token), grand_avg_dir);
    end

    if ~isempty(Gimg.power)
        mean_img = mean(Gimg.power, 3, 'omitnan');
        plot_save_TFR(mean_img, freqs, Gimg.times_ref, sprintf('grand_wavelet_imagined_%s', token), grand_avg_dir);
    end

    if ~isempty(Gexec.power) && ~isempty(Gimg.power)
        plot_save_TFR(mean_exec - mean_img, freqs, Gexec.times_ref, sprintf('grand_wavelet_diff_%s', token), grand_avg_dir);
    end

    % ---- Hilbert Grand Averages ----
    if ~isempty(Gexec.power_hilbert)
        mean_exec_h = mean(Gexec.power_hilbert, 3, 'omitnan');
        plot_save_TFR(mean_exec_h, freqs, Gexec.times_ref, sprintf('grand_hilbert_executed_%s', token), grand_avg_dir);
    end

    if ~isempty(Gimg.power_hilbert)
        mean_img_h = mean(Gimg.power_hilbert, 3, 'omitnan');
        plot_save_TFR(mean_img_h, freqs, Gimg.times_ref, sprintf('grand_hilbert_imagined_%s', token), grand_avg_dir);
    end

    if ~isempty(Gexec.power_hilbert) && ~isempty(Gimg.power_hilbert)
        plot_save_TFR(mean_exec_h - mean_img_h, freqs, Gexec.times_ref, sprintf('grand_hilbert_diff_%s', token), grand_avg_dir);
    end
end

fprintf('\nAll done. Per-subject and grand-average plots saved under:\n  %s\n  %s\n  %s\n', ...
    subject_outdir, hilbert_outdir, grand_avg_dir);

%% ===================================================================
%                       HELPER FUNCTION
% ===================================================================
function plot_save_TFR(mat, freqs, times, fname, outdir)
    % Plot and save a single TFR figure (freq × time)
    h = figure('units','normalized','position',[0.2 0.2 0.6 0.45],'visible','off');
    imagesc(times, freqs, mat);
    axis xy;
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title(strrep(fname,'_',' '));
    c = colorbar; ylabel(c,'dB');
    colormap('jet');
    saveas(h, fullfile(outdir, [fname '.png']));
    close(h);
end

