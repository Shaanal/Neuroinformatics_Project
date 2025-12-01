%% ===============================================================
%  PHASE 3 : SPECTRAL ANALYSIS (PSD & BANDPOWER)
%  ===============================================================
clear; clc;

% ---------- CONFIGURATION ----------
erp_folder = './results/erp/';
psd_folder = './results/psd/';
if ~exist(psd_folder,'dir'), mkdir(psd_folder); end

Fs = 160;                           % sampling rate (Hz)
channels_to_plot = [9 11 13];      % approximate C3, Cz, C4
mu_band = [8 13];                  % μ (mu) band
beta_band = [13 30];               % β (beta) band

erp_files = dir(fullfile(erp_folder,'S*_ERP.mat'));
fprintf('Loaded %d ERP files for PSD analysis.\n', length(erp_files));

% ---------- Initialize storage ----------
all_mu_real = []; all_mu_imag = [];
all_beta_real = []; all_beta_imag = [];
valid_subjects = {};

%% ---------- LOOP OVER SUBJECTS ----------
for i = 1:length(erp_files)
    fname = fullfile(erp_folder, erp_files(i).name);
    load(fname, 'ERP_struct');
    subj = extractBefore(erp_files(i).name,'_ERP.mat');

    % ---- check and assign Real/Imagined conditions ----
    if isfield(ERP_struct,'Real_LeftRight_Fist')
        real = ERP_struct.Real_LeftRight_Fist;
    elseif isfield(ERP_struct,'Real_Hands')
        real = ERP_struct.Real_Hands;
    else
        fprintf('Skipping %s (no Real condition)\n', subj);
        continue;
    end

    if isfield(ERP_struct,'Imagined_LeftRight_Fist')
        imag = ERP_struct.Imagined_LeftRight_Fist;
    elseif isfield(ERP_struct,'Imagined_Hands')
        imag = ERP_struct.Imagined_Hands;
    else
        fprintf('Skipping %s (no Imagined condition)\n', subj);
        continue;
    end

    % ---- average movement-related ERPs (T1 and T2 represent movement events) ----
    real_avg = mean(cat(3, real.T1_mean, real.T2_mean), 3, 'omitnan');
    imag_avg = mean(cat(3, imag.T1_mean, imag.T2_mean), 3, 'omitnan');

    % ---- handle cases where Rest (T0) may have fewer channels ----
    min_ch_real = min(size(real_avg,1), size(real.T0_mean,1));
    min_ch_imag = min(size(imag_avg,1), size(imag.T0_mean,1));

    % ---- subtract baseline (T0) for real/imag conditions ----
    real_sub = real_avg(1:min_ch_real,:) - real.T0_mean(1:min_ch_real,:);
    imag_sub = imag_avg(1:min_ch_imag,:) - imag.T0_mean(1:min_ch_imag,:);

    % ---- find motor channels that still exist post-cleaning ----
    valid_real_ch = channels_to_plot(channels_to_plot <= min_ch_real);
    valid_imag_ch = channels_to_plot(channels_to_plot <= min_ch_imag);

    % ---- average motor cortex channels for PSD input ----
    real_sig = mean(real_sub(valid_real_ch,:), 1, 'omitnan');
    imag_sig = mean(imag_sub(valid_imag_ch,:), 1, 'omitnan');

    % ---- Subjects with only NaN data are unusable ----
    if all(isnan(real_sig)) || all(isnan(imag_sig))
        fprintf('Skipping %s (signal all NaN)\n', subj);
        continue;
    end

    % Replace occasional NaNs after cleaning/interpolation
    real_sig(isnan(real_sig)) = 0;
    imag_sig(isnan(imag_sig)) = 0;

    % ---- compute PSD via Welch’s method ----
    % hamming(512) = window size; 256 = 50% overlap; 1024 FFT points
    [Pxx_real,freqs] = pwelch(real_sig, hamming(512), 256, 1024, Fs);
    [Pxx_imag,~]     = pwelch(imag_sig, hamming(512), 256, 1024, Fs);

    % ---- compute power in μ and β frequency bands ----
    mu_real   = bandpower(Pxx_real, freqs, mu_band, 'psd');
    mu_imag   = bandpower(Pxx_imag, freqs, mu_band, 'psd');
    beta_real = bandpower(Pxx_real, freqs, beta_band, 'psd');
    beta_imag = bandpower(Pxx_imag, freqs, beta_band, 'psd');

    % Store subject-level bandpower values
    all_mu_real   = [all_mu_real; mu_real];
    all_mu_imag   = [all_mu_imag; mu_imag];
    all_beta_real = [all_beta_real; beta_real];
    all_beta_imag = [all_beta_imag; beta_imag];
    valid_subjects{end+1} = subj;

    % ---- per-subject PSD plot ----
    fig = figure('Visible','off','Position',[200 200 900 400]);
    plot(freqs,10*log10(Pxx_real),'r','LineWidth',1.5); hold on;
    plot(freqs,10*log10(Pxx_imag),'b','LineWidth',1.5);

    xlabel('Frequency (Hz)'); 
    ylabel('Power (dB/Hz)');
    title(sprintf('%s PSD – Real vs Imagined',subj));
    legend({'Real','Imagined'}); 
    grid on; 
    xlim([0 40]);

    saveas(fig, fullfile(psd_folder,[subj '_PSD.png']));
    close(fig);
end

%% ---------- GROUP-LEVEL STATISTICS ----------
fprintf('\nPerforming group statistics...\n');
validN = length(valid_subjects);

if validN < 3
    error('Too few valid subjects for statistical analysis.');
end

% Paired t-tests: Real vs Imagined bandpower differences
[~,p_mu,~,stats_mu]     = ttest(all_mu_real, all_mu_imag);
[~,p_beta,~,stats_beta] = ttest(all_beta_real, all_beta_imag);

fprintf('μ-band  (8–13 Hz):  t(%d)=%.2f,  p=%.4f\n', stats_mu.df, stats_mu.tstat, p_mu);
fprintf('β-band (13–30 Hz): t(%d)=%.2f,  p=%.4f\n', stats_beta.df, stats_beta.tstat, p_beta);

% ---- grouped bar plots for Real vs Imagined ----
fig2 = figure('Visible','off','Position',[300 300 700 300]);

subplot(1,2,1);
bar([mean(all_mu_real) mean(all_mu_imag)], 'FaceColor','flat'); hold on;
errorbar(1:2,[mean(all_mu_real) mean(all_mu_imag)], ...
    [std(all_mu_real)/sqrt(validN), std(all_mu_imag)/sqrt(validN)], '.k');
set(gca,'XTickLabel',{'Real','Imagined'});
ylabel('μ Power (µV²)');
title(sprintf('μ-band Power (p=%.4f)',p_mu)); 
grid on;

subplot(1,2,2);
bar([mean(all_beta_real) mean(all_beta_imag)], 'FaceColor','flat'); hold on;
errorbar(1:2,[mean(all_beta_real) mean(all_beta_imag)], ...
    [std(all_beta_real)/sqrt(validN), std(all_beta_imag)/sqrt(validN)], '.k');
set(gca,'XTickLabel',{'Real','Imagined'});
ylabel('β Power (µV²)');
title(sprintf('β-band Power (p=%.4f)',p_beta)); 
grid on;

saveas(fig2, fullfile(psd_folder,'Group_Bandpower_Barplots.png'));
close(fig2);

%% ---------- SAVE SUMMARY ----------
summary_tbl = table(valid_subjects', all_mu_real, all_mu_imag, all_beta_real, all_beta_imag, ...
    'VariableNames', {'Subject','Mu_Real','Mu_Imagined','Beta_Real','Beta_Imagined'});

writetable(summary_tbl, fullfile(psd_folder,'PSD_Bandpower_Summary.csv'));

fprintf('Saved bandpower summary and figures to %s\n', psd_folder);

