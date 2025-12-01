%% ===============================================================
%  PHASE 5_3: FINAL GROUP STATISTICS FOR TFR + HILBERT + ITPC
% ===============================================================

clear; clc;

%% ---------------- CONFIG ----------------
root = "./results/tfr_new_results/";       % parent folder for TFR results
tfr_path = fullfile(root, "tfr_subjects"); % per-subject wavelet files
hilbert_path = fullfile(root, "hilbert_subjects"); % per-subject Hilbert files
out_path = fullfile(root, "grand_averages");       % output folder for summary

if ~exist(out_path,'dir'), mkdir(out_path); end

alpha_band = [8 12];   % classical mu-alpha
beta_band  = [13 30];  % motor beta

fprintf("\n==============================================\n");
fprintf("      FINAL GROUP TFR / HILBERT STATS\n");
fprintf("==============================================\n");

%% ---------------- FIND ALL FILES ----------------
all_files = dir(fullfile(tfr_path, "*.mat")); % list all TFR .mat files

fprintf("\nFound %d TFR MAT files.\n", numel(all_files));

if numel(all_files)==0
    error("No .mat files found! Run preprocessing first.");
end

% Show filenames to help debug naming issues
fprintf("\nAll MAT files found:\n");
for i = 1:numel(all_files)
    fprintf("  %d: %s\n", i, all_files(i).name);
end

%% ---- Extract unique subject IDs using regex ----
names = {all_files.name};
tokens = regexp(names, "^(S\d+)_", "tokens");   % extract S### prefix
subjects = cell(0);

% Collect only valid subject IDs
for i = 1:numel(tokens)
    if ~isempty(tokens{i})
        subjects{end+1} = tokens{i}{1}{1};
    end
end

subjects = unique(subjects);  % remove duplicates
fprintf("\nUnique subjects found: %d\n", numel(subjects));
disp(subjects');

%% ------------- STORAGE FOR GROUP DATA -------------
summary = struct([]);   % struct-array for final CSV

%% ==========================================================
%                PROCESS EACH SUBJECT
% ==========================================================
for s = 1:numel(subjects)

    subj = subjects{s};
    fprintf("\n------------------------------------------\n");
    fprintf("Processing subject: %s\n", subj);

    % File patterns for Real and Imag TFR files (subject_R##_executed_tfr.mat)
    pat_real = sprintf("%s_R*_executed_tfr.mat", subj);
    pat_imag = sprintf("%s_R*_imagined_tfr.mat", subj);

    % Find matching files
    file_real = dir(fullfile(tfr_path, pat_real));
    file_imag = dir(fullfile(tfr_path, pat_imag));
    
    % Debug output to confirm correct file detection
    if isempty(file_real)
        fprintf("  No Real files found matching: %s\n", pat_real);
    else
        fprintf("  Found Real: %s\n", file_real(1).name);
    end
    if isempty(file_imag)
        fprintf("  No Imagined files found matching: %s\n", pat_imag);
    else
        fprintf("  Found Imagined: %s\n", file_imag(1).name);
    end

    % Cannot compute features if one condition is missing
    if isempty(file_real) || isempty(file_imag)
        fprintf("Skipping %s — missing Real or Imagined.\n", subj);
        continue;
    end

    % Load first matching Real and Imag files
    realmat = load(fullfile(file_real(1).folder, file_real(1).name));
    imagmat = load(fullfile(file_imag(1).folder, file_imag(1).name));

    % Each MAT contains subj_tfr.C3 (or other channels)
    R = realmat.subj_tfr;
    I = imagmat.subj_tfr;

    % Require at least C3 for analysis
    if ~isfield(R,'C3') || ~isfield(I,'C3')
        fprintf("Skipping %s — C3 missing.\n", subj);
        continue;
    end

    freqs = R.C3.freqs;   % frequency vector
    times = R.C3.times;   % time vector

    % Frequency band index masks
    aidx = freqs>=alpha_band(1) & freqs<=alpha_band(2);
    bidx = freqs>=beta_band(1)  & freqs<=beta_band(2);

    % Wavelet power (already baseline-corrected in dB)
    RA = mean(R.C3.power_db(aidx,:), "all", "omitnan");
    RB = mean(R.C3.power_db(bidx,:), "all", "omitnan");
    IA = mean(I.C3.power_db(aidx,:), "all", "omitnan");
    IB = mean(I.C3.power_db(bidx,:), "all", "omitnan");

    % Default Hilbert values in case files are missing
    HA = NaN; HB = NaN; HIA = NaN; HIB = NaN;

    % Correct patterns for Hilbert MAT files
    h_real = dir(fullfile(hilbert_path, sprintf("%s_R*_executed_hilbert.mat", subj)));
    h_imag = dir(fullfile(hilbert_path, sprintf("%s_R*_imagined_hilbert.mat", subj)));

    % Hilbert extraction if both files exist
    if ~isempty(h_real) && ~isempty(h_imag)
        H_R = load(fullfile(h_real(1).folder, h_real(1).name));
        H_I = load(fullfile(h_imag(1).folder, h_imag(1).name));

        % Check channel inside subj_hilb struct
        if isfield(H_R.subj_hilb,'C3')
            HA  = mean(H_R.subj_hilb.C3.power_db(aidx,:), "all", "omitnan");
            HB  = mean(H_R.subj_hilb.C3.power_db(bidx,:), "all", "omitnan");
            HIA = mean(H_I.subj_hilb.C3.power_db(aidx,:), "all", "omitnan");
            HIB = mean(H_I.subj_hilb.C3.power_db(bidx,:), "all", "omitnan");
        end
    end

    % Add all extracted values into one row
    summary(end+1).Subject = subj;
    summary(end).RealAlpha = RA;
    summary(end).ImagAlpha = IA;
    summary(end).RealBeta  = RB;
    summary(end).ImagBeta  = IB;

    summary(end).HilbRealAlpha = HA;
    summary(end).HilbImagAlpha = HIA;
    summary(end).HilbRealBeta  = HB;
    summary(end).HilbImagBeta  = HIB;

end

%% ============================================================
%                 MAKE TABLE + SAVE CSV
% ============================================================
if isempty(summary)
    error("No subjects had both Real + Imagined data.");
end

T = struct2table(summary);                  % convert to table
csvfile = fullfile(out_path, "Group_TFR_Features.csv");
writetable(T, csvfile);

fprintf("\n==========================================\n");
fprintf("Saved GROUP CSV: %s\n", csvfile);
fprintf("==========================================\n");

%% ============================================================
%            GROUP-LEVEL STATISTICS (Paired t-tests)
% ============================================================
RA = T.RealAlpha;  IA = T.ImagAlpha;
RB = T.RealBeta;   IB = T.ImagBeta;

fprintf("\n------- Wavelet Band Power Statistics (C3) -------\n");

% Need >=3 subjects for meaningful t-test
if height(T) >= 3
    [~,pA,~,sA] = ttest(RA, IA);
    [~,pB,~,sB] = ttest(RB, IB);

    fprintf("Alpha: t(%d)=%.3f, p=%.5f\n", sA.df, sA.tstat, pA);
    fprintf("Beta : t(%d)=%.3f, p=%.5f\n", sB.df, sB.tstat, pB);
else
    fprintf("Not enough subjects for t-tests.\n");
end

%% Hilbert Stats
fprintf("\n------- Hilbert Band Power Statistics (C3) -------\n");

HA  = T.HilbRealAlpha;
HIA = T.HilbImagAlpha;
HB  = T.HilbRealBeta;
HIB = T.HilbImagBeta;

% Only compute if at least 3 non-NaN Hilbert values exist
if height(T) >= 3 && sum(~isnan(HA)) >= 3
    [~,pHA,~,sHA] = ttest(HA, HIA);
    [~,pHB,~,sHB] = ttest(HB, HIB);

    fprintf("Hilbert Alpha: t(%d)=%.3f, p=%.5f\n", sHA.df, sHA.tstat, pHA);
    fprintf("Hilbert Beta : t(%d)=%.3f, p=%.5f\n", sHB.df, sHB.tstat, pHB);
else
    fprintf("Not enough Hilbert subjects.\n");
end

%% ============================================================
%                   PLOT GROUP RESULTS
% ============================================================
fig = figure('Position',[300 200 900 400]);

% Wavelet alpha
subplot(1,2,1)
bar([mean(RA) mean(IA)]); hold on
errorbar(1:2, [mean(RA) mean(IA)], ...
         [std(RA)/sqrt(height(T)) std(IA)/sqrt(height(T))],'k.');
set(gca,'XTickLabel',{'Real','Imag'}); ylabel("Alpha Power (dB)")
title("Wavelet Alpha")

% Wavelet beta
subplot(1,2,2)
bar([mean(RB) mean(IB)]); hold on
errorbar(1:2, [mean(RB) mean(IB)], ...
         [std(RB)/sqrt(height(T)) std(IB)/sqrt(height(T))],'k.');
set(gca,'XTickLabel',{'Real','Imag'}); ylabel("Beta Power (dB)")
title("Wavelet Beta")

saveas(fig, fullfile(out_path,"Wavelet_Barplots.png"));

fig2 = figure('Position',[300 200 900 400]);

% Boxplot distributions
subplot(1,2,1)
boxplot([RA IA], 'Labels', {'Real','Imag'}); ylabel("Alpha (dB)")
title("Wavelet Alpha Distribution")

subplot(1,2,2)
boxplot([RB IB], 'Labels', {'Real','Imag'}); ylabel("Beta (dB)")
title("Wavelet Beta Distribution")

saveas(fig2, fullfile(out_path,"Wavelet_Boxplots.png"));

fprintf("\nPlots saved.\n");

