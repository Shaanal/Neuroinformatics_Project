%% ===============================================================
%  PHASE 2 : ERP Visualization + Feature Extraction + Statistics
%  ===============================================================
clear; clc;

% ---------- CONFIGURATION ----------
erp_folder = './results/erp/';
fig_folder = './results/figures/';
if ~exist(fig_folder,'dir'), mkdir(fig_folder); end

channels_to_plot = [9 11 13];   % approx. C3, Cz, C4 (motor cortex)
time_window = [0 1];            % window for peak extraction (post-cue)
Fs = 160;                       % sampling rate (used for latency conversion)
erp_files = dir(fullfile(erp_folder,'S*_ERP.mat'));

all_real = [];
all_imag = [];
subj_ids = {};
valid_subj_ids = {};  % subjects where both Real + Imagined were processed

fprintf('Loaded %d ERP files.\n\n', length(erp_files));

%% ---------- LOOP OVER SUBJECTS ----------
for i = 1:length(erp_files)
    fname = fullfile(erp_folder, erp_files(i).name);
    load(fname, 'ERP_struct');
    
    % Extract subject ID (e.g., S001)
    subj = extractBefore(erp_files(i).name,'_ERP.mat');
    subj_ids{end+1} = subj;

    % ---- Display what fields exist (for debugging/inconsistencies) ----
    fn = fieldnames(ERP_struct);
    disp(['Fields in ' subj ': ' strjoin(fn', ', ')]);

    % ---- Detect Real & Imagined blocks ----
    hasReal = isfield(ERP_struct,'Real_LeftRight_Fist') || ...
              isfield(ERP_struct,'Real_Hands');

    hasImag = isfield(ERP_struct,'Imagined_LeftRight_Fist') || ...
              isfield(ERP_struct,'Imagined_Hands');

    if ~(hasReal && hasImag)
        fprintf('Skipping %s (missing Real or Imagined fields)\n', subj);
        continue; % cannot analyze this subject
    end

    % ---- Assign Real/Imag structures based on field name present ----
    if isfield(ERP_struct,'Real_LeftRight_Fist')
        real = ERP_struct.Real_LeftRight_Fist;
    else
        real = ERP_struct.Real_Hands;
    end

    if isfield(ERP_struct,'Imagined_LeftRight_Fist')
        imag = ERP_struct.Imagined_LeftRight_Fist;
    else
        imag = ERP_struct.Imagined_Hands;
    end

    % ---- Average movement-related ERPs (T1 & T2 correspond to responses) ----
    real_avg = mean(cat(3, real.T1_mean, real.T2_mean), 3, 'omitnan');
    imag_avg = mean(cat(3, imag.T1_mean, imag.T2_mean), 3, 'omitnan');

    % ---- Compute difference waves by subtracting the rest (T0 baseline) ----
    real_diff = real_avg - real.T0_mean;
    imag_diff = imag_avg - imag.T0_mean;

    % ---- Time axis (from preprocessing phase: epoch = [-1, 4] sec) ----
    t = linspace(-1,4,size(real_avg,2));

    % ---- Create ERP plot for Real vs Imagined ----
    fig = figure('Name',subj,'Position',[200 200 900 500],'Visible','off');
    hold on;

    % Verify requested channels exist (after cleaning, some channels may be removed)
    valid_channels = channels_to_plot(channels_to_plot <= size(real_diff,1));
    if isempty(valid_channels)
        valid_channels = size(real_diff,1); % fallback to last existing channel
        fprintf('%s has only %d channels. Using last available channel.\n',...
                subj, size(real_diff,1));
    end

    % ---- Plot Real and Imagined ERPs for selected channels ----
    for ch = valid_channels
        plot(t, real_diff(ch,:), 'r', 'LineWidth',1.2);
        plot(t, imag_diff(ch,:), 'b', 'LineWidth',1.2);
    end

    xlabel('Time (s)'); ylabel('Amplitude (µV)');
    title(sprintf('%s: ERP (Real vs Imagined) C3–Cz–C4', subj));
    legend({'Real','Imagined'}); 
    grid on; 
    xlim([-0.5 2]);

    saveas(fig, fullfile(fig_folder,[subj '_ERP_RealVsImagined.png']));
    close(fig);

    % ---- FEATURE EXTRACTION ----
    % We extract peak amplitude in the 0–1s window after cue.
    idx = t>=time_window(1) & t<=time_window(2);

    % Find peak amplitude within the window for each channel
    [peak_real,lat_real] = max(abs(real_diff(:,idx)),[],2);
    [peak_imag,lat_imag] = max(abs(imag_diff(:,idx)),[],2);

    % Convert latency index to seconds (t was defined from -1s)
    lat_real = (lat_real + find(idx,1)-1)/Fs - 1;
    lat_imag = (lat_imag + find(idx,1)-1)/Fs - 1;

    % Use the mean peak value across C3–Cz–C4
    all_real = [all_real; mean(peak_real(valid_channels))];
    all_imag = [all_imag; mean(peak_imag(valid_channels))];

    valid_subj_ids{end+1} = subj; % mark subject as valid
end

%% ---------- GROUP-LEVEL STATISTICS ----------
validN = length(all_real);
if validN < 3
    error('Too few subjects with both Real & Imagined for statistical comparison.');
end

% Paired t-test between real and imagined ERPs across subjects
[~,p_ttest,~,stats] = ttest(all_real, all_imag);

fprintf('\nGroup Paired t-test Real vs Imagined (mean peak ERP amplitude):\n');
fprintf('  Subjects included: %d\n', validN);
fprintf('  t(%d)=%.2f,  p=%.4f\n', stats.df, stats.tstat, p_ttest);

%% ---------- GROUP BAR PLOT ----------
fig2 = figure('Position',[300 300 600 400],'Visible','off');
bar([mean(all_real) mean(all_imag)], 'FaceColor','flat');
hold on;

% Add error bars (SEM)
errorbar(1:2, [mean(all_real) mean(all_imag)], ...
         [std(all_real)/sqrt(validN), std(all_imag)/sqrt(validN)], ...
         '.k','LineWidth',1.2);

set(gca,'XTickLabel',{'Real','Imagined'});
ylabel('Peak ERP amplitude (µV)');
title(sprintf('Group Comparison (p = %.4f)',p_ttest));
grid on;

saveas(fig2, fullfile(fig_folder,'Group_ERP_Barplot.png'));
close(fig2);

%% ---------- SAVE FEATURE SUMMARY ----------
summary_tbl = table(valid_subj_ids', all_real, all_imag, ...
    'VariableNames', {'Subject','Real_PeakAmp','Imagined_PeakAmp'});

writetable(summary_tbl,'./results/erp/ERP_Features.csv');

fprintf('Saved ERP feature summary to results/erp/ERP_Features.csv\n');
fprintf('Figures saved in %s\n', fig_folder);

