%% ========================================================================
%   PHASE 4: GROUP ERSP ANALYSIS (20 SUBJECTS)
% ========================================================================

clear; clc; close all;
[ALLEEG, EEG, CURRENTSET] = eeglab;

%% ========================== CONFIG ==========================
processed_root = './processed/';
results_root   = './group_results_ersp/';
if ~exist(results_root,'dir'), mkdir(results_root); end

N_pick = 20;     % number of subjects to sample

baseline_ms = [-2000 -500];    % baseline for ERSP (ms)
task_window = [0 1500];        % time window (ms) to average ERSP over task
freq_range  = [4 40];          % frequency range for time-frequency decomposition
cycles      = [3 0.5];         % wavelet cycles parameter (as in your S001 call)

alpha_band  = [8 12];
beta_band   = [13 30];

subdirs = dir(fullfile(processed_root, 'S*'));
subdirs = subdirs([subdirs.isdir]);
subjects = {subdirs.name};

% Ensuring only subject folders matching pattern S### are considered
subjects = subjects(~cellfun(@isempty, regexp(subjects, '^S\d+$')));

fprintf('\nFound %d valid subjects.\n', length(subjects));

if length(subjects) < N_pick
    N_pick = length(subjects);
end

rng(42); % reproducible random sampling
sel_idx = randperm(length(subjects), N_pick);
picked_subjects = subjects(sel_idx);

fprintf('\nRandomly selected %d subjects:\n', N_pick);
disp(picked_subjects');

%% ========== STORAGE FOR GROUP ANALYSIS ================
group_ids = {};
real_alpha = [];    imag_alpha = [];
real_beta  = [];    imag_beta  = [];

%% ====================== PROCESS EVERY SUBJECT ======================

for s = 1:N_pick
    
    subj = picked_subjects{s};
    subj_path = fullfile(processed_root, subj);
    
    fprintf('\n\n==========================================\n');
    fprintf('             Subject: %s\n', subj);
    fprintf('==========================================\n');
    
    % Find processed .set files for this subject
    run_files = dir(fullfile(subj_path, '*_processed.set'));
    if isempty(run_files)
        fprintf('No SET files for %s. Skipping.\n', subj);
        continue;
    end
    
    % Per-subject accumulators for runs
    subj_real_alpha = []; 
    subj_imag_alpha = [];
    subj_real_beta  = [];
    subj_imag_beta  = [];
    
    for i = 1:length(run_files)
        
        fname = run_files(i).name;
        % Extract run number (R##) from filename
        tok = regexp(fname, 'R(\d+)', 'tokens','once');
        if isempty(tok), continue; end
        R = str2double(tok{1});
        
        % Determine whether this run belongs to Real or Imag conditions
        if ismember(R, [3 7 11])
            condition = 'Real';
        elseif ismember(R, [4 8 12])
            condition = 'Imag';
        else
            % Skip runs that are not part of the motor tasks of interest
            continue;
        end
        
        fprintf('Loading %s (%s)\n', fname, condition);
        
        EEG = pop_loadset('filename',fname,'filepath',subj_path);
        
        % Dynamically find C3 channel 
        labels = {EEG.chanlocs.labels};
        chanC3 = find(strcmpi(labels,'C3') | contains(labels,'C3'));
        
        if isempty(chanC3)
            % If C3 not found, skip this run (we focus on C3 for group ERSP)
            fprintf('No C3 channel found — skipping run.\n');
            continue;
        end
        
        % Extract single-channel data for C3 (channels × times × trials -> times × trials)
        % squeeze to get (times × trials)
        dataC3 = squeeze(EEG.data(chanC3,:,:));
        
        % Computing ERSP 
        try
            [erspC3, ~, ~, times, freqs] = newtimef( ...
                dataC3, EEG.pnts, ...
                [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles, ...
                'baseline', baseline_ms, ...
                'freqs', freq_range, ...
                'plotersp','off','plotitc','off','verbose','off');
        catch
            fprintf('ERSP failed — skipping run.\n');
            continue;
        end
        
        % Create boolean indices for alpha/beta/time windows
        aidx = freqs>=alpha_band(1) & freqs<=alpha_band(2);
        bidx = freqs>=beta_band(1)  & freqs<=beta_band(2);
        tidx = times>=task_window(1) & times<=task_window(2);
        
        % If any index is empty, the requested band or time window isn't present -> skip
        if ~any(aidx) || ~any(bidx) || ~any(tidx)
            fprintf('Invalid index (no matching freq/time). Skipping run.\n');
            continue;
        end
        
        % Average ERSP over selected freq × time window → scalar summary per run
        alpha_val = mean(erspC3(aidx, tidx), 'all', 'omitnan');
        beta_val  = mean(erspC3(bidx, tidx), 'all', 'omitnan');
        
        % Accumulate per-condition values for this subject
        if strcmp(condition,'Real')
            subj_real_alpha(end+1) = alpha_val;
            subj_real_beta(end+1)  = beta_val;
        else
            subj_imag_alpha(end+1) = alpha_val;
            subj_imag_beta(end+1)  = beta_val;
        end
    end
    
    % Require both Real and Imag data for the subject to be useful in group analysis
    if isempty(subj_real_alpha) || isempty(subj_imag_alpha)
        fprintf('Subject missing Real or Imag data — skipping.\n');
        continue;
    end
    
    % Average across runs for this subject and add to group arrays
    group_ids{end+1} = subj;
    real_alpha(end+1) = mean(subj_real_alpha);
    imag_alpha(end+1) = mean(subj_imag_alpha);
    real_beta(end+1)  = mean(subj_real_beta);
    imag_beta(end+1)  = mean(subj_imag_beta);
    
    fprintf('%s: Alpha(Real=%.2f | Imag=%.2f), Beta(Real=%.2f | Imag=%.2f)\n', ...
        subj, real_alpha(end), imag_alpha(end), real_beta(end), imag_beta(end));
end

%% ==================== GROUP ANALYSIS =========================

n = length(real_alpha);
fprintf('\n\n==========================================\n');
fprintf('         GROUP RESULTS (%d subjects)\n', n);
fprintf('==========================================\n');

if n < 3
    fprintf('Not enough subjects for t-test (need ≥ 3).\n');
    p_alpha = NaN;
    p_beta  = NaN;
else
    % Paired t-tests: Real vs Imag for alpha and beta band ERSP
    [~,p_alpha,~,stats_alpha] = ttest(real_alpha, imag_alpha);
    [~,p_beta ,~,stats_beta ] = ttest(real_beta , imag_beta );
    
    fprintf('Alpha: t(%d)=%.2f, p=%.4f\n', stats_alpha.df, stats_alpha.tstat, p_alpha);
    fprintf('Beta : t(%d)=%.2f, p=%.4f\n', stats_beta.df, stats_beta.tstat, p_beta);
end

%% ========================= SAVE CSV ==========================
% Save per-subject group summaries for downstream reporting
T = table(group_ids', real_alpha', imag_alpha', real_beta', imag_beta',...
    'VariableNames', {'Subject','RealAlpha','ImagAlpha','RealBeta','ImagBeta'});
writetable(T, fullfile(results_root,'Group_ERSP_Results.csv'));

%% ========================= BARPLOTS ==========================
fig = figure('Position',[200 200 900 400],'Color','w');

subplot(1,2,1)
% Bar plot of mean alpha ERSP (Real vs Imag)
bar([mean(real_alpha) mean(imag_alpha)],'FaceColor',[0.6 0.8 1]); hold on;
if n>1
    % Plot SEM errorbars when more than one subject
    errorbar(1:2,[mean(real_alpha) mean(imag_alpha)], ...
             [std(real_alpha)/sqrt(n), std(imag_alpha)/sqrt(n)], 'k.','LineWidth',1.5);
end
set(gca,'XTickLabel',{'Real','Imag'}); ylabel('Alpha ERSP (dB)');
title('Alpha Band'); grid on;

subplot(1,2,2)
% Bar plot of mean beta ERSP (Real vs Imag)
bar([mean(real_beta) mean(imag_beta)],'FaceColor',[1 0.8 0.6]); hold on;
if n>1
    errorbar(1:2,[mean(real_beta) mean(imag_beta)], ...
             [std(real_beta)/sqrt(n), std(imag_beta)/sqrt(n)], 'k.','LineWidth',1.5);
end
set(gca,'XTickLabel',{'Real','Imag'}); ylabel('Beta ERSP (dB)');
title('Beta Band'); grid on;

saveas(fig, fullfile(results_root,'Group_ERSP_Barplots.png'));
close(fig);

%% ========================= BOXPLOTS ==========================
% Boxplots to visualize distribution across subjects (only when >=2)
if n>=2
    fig2 = figure('Position',[200 200 900 400],'Color','w');
    
    subplot(1,2,1)
    boxplot([real_alpha' imag_alpha'], {'Real','Imag'}); ylabel('Alpha ERSP (dB)');
    title('Alpha Distribution'); grid on;
    
    subplot(1,2,2)
    boxplot([real_beta' imag_beta'], {'Real','Imag'}); ylabel('Beta ERSP (dB)');
    title('Beta Distribution'); grid on;
    
    saveas(fig2, fullfile(results_root,'Group_ERSP_Boxplots.png'));
    close(fig2);
end

fprintf('\n Group analysis complete. Results in: %s\n', results_root);

