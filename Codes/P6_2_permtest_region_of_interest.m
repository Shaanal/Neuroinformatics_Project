%% ===============================================================
%  PHASE 6_2: ROI-BASED FOCUSED STATISTICAL ANALYSIS
% ===============================================================

clear; clc;

%% ---------------- CONFIG ----------------
root = "./results/tfr_new_results/";
tfr_path = fullfile(root, "tfr_subjects");
out_path = fullfile(root, "roi_stats");
if ~exist(out_path,'dir'), mkdir(out_path); end

% ROIs defined based on classical sensorimotor ERD/ERS timing and frequencies
ROIs = struct();

ROIs(1).name = 'Alpha_Movement';
ROIs(1).freq = [8 12];
ROIs(1).time = [0 1.5];
ROIs(1).hypothesis = 'Real movement should show stronger alpha suppression';

ROIs(2).name = 'Beta_Movement';
ROIs(2).freq = [13 30];
ROIs(2).time = [0 1.5];
ROIs(2).hypothesis = 'Real movement should show stronger beta ERD';

ROIs(3).name = 'Beta_Rebound';
ROIs(3).freq = [13 30];
ROIs(3).time = [1.5 3.5];
ROIs(3).hypothesis = 'Post-movement beta rebound (ERS) expected';

ROIs(4).name = 'Alpha_Preparation';
ROIs(4).freq = [8 12];
ROIs(4).time = [-0.5 0];
ROIs(4).hypothesis = 'Pre-movement motor preparation';

% Permutation settings (non-parametric robustness)
n_permutations = 1000;
alpha_level = 0.05;

fprintf("\n==============================================\n");
fprintf("   ROI-BASED FOCUSED ANALYSIS\n");
fprintf("==============================================\n");
fprintf("Bonferroni α = %.4f\n", alpha_level/length(ROIs));

%% ---------------- LOAD ALL SUBJECT DATA ----------------
all_files = dir(fullfile(tfr_path, "*.mat"));
names = {all_files.name};
tokens = regexp(names, "^(S\d+)_", "tokens");

subjects = {};
for i = 1:numel(tokens)
    if ~isempty(tokens{i}), subjects{end+1} = tokens{i}{1}{1}; end
end
subjects = unique(subjects);

real_data = [];
imag_data = [];
freqs = [];
times = [];

fprintf("\nLoading subjects: ");
for s = 1:numel(subjects)
    subj = subjects{s};

    % match run-level TFR files
    fR = dir(fullfile(tfr_path, sprintf("%s_R*_executed_tfr.mat", subj)));
    fI = dir(fullfile(tfr_path, sprintf("%s_R*_imagined_tfr.mat", subj)));
    if isempty(fR) || isempty(fI), continue; end

    R = load(fullfile(fR(1).folder, fR(1).name));
    I = load(fullfile(fI(1).folder, fI(1).name));

    % Only C3 used → strongest MI signature
    if ~isfield(R.subj_tfr,'C3') || ~isfield(I.subj_tfr,'C3'), continue; end

    if isempty(freqs)
        freqs = R.subj_tfr.C3.freqs;
        times = R.subj_tfr.C3.times;
    end

    % Store entire TFR (freq × time)
    real_data(end+1,:,:) = R.subj_tfr.C3.power_db;
    imag_data(end+1,:,:) = I.subj_tfr.C3.power_db;

    fprintf(".");
end
fprintf("\n");

n_subjects = size(real_data,1);
fprintf("Loaded %d subjects\n", n_subjects);
if n_subjects < 3, error("Need ≥ 3 subjects."); end

%% ============================================================
%                   EXTRACT ROI VALUES
% ============================================================
roi_real = zeros(n_subjects, length(ROIs));
roi_imag = zeros(n_subjects, length(ROIs));

fprintf("\nExtracting ROI values...\n");

for r = 1:length(ROIs)

    % ROI → convert freq/time windows into indices
    fidx = freqs >= ROIs(r).freq(1) & freqs <= ROIs(r).freq(2);
    tidx = times >= ROIs(r).time(1) & times <= ROIs(r).time(2);

    fprintf("  ROI %d (%s): using %d freq x %d time points\n", ...
        r, ROIs(r).name, sum(fidx), sum(tidx));

    % Mean across ROI area for each subject
    for s = 1:n_subjects
        roi_real(s,r) = mean(real_data(s,fidx,tidx), 'all');
        roi_imag(s,r) = mean(imag_data(s,fidx,tidx), 'all');
    end
end

%% ============================================================
%          PARAMETRIC T-TESTS PER ROI
% ============================================================
fprintf("\n==============================================\n");
fprintf("       PARAMETRIC TESTS (REAL vs IMAG)\n");
fprintf("==============================================\n");

results_parametric = struct();

for r = 1:length(ROIs)
    fprintf("\nROI %d: %s\n", r, ROIs(r).name);

    [~, p, ~, stats] = ttest(roi_real(:,r), roi_imag(:,r));

    diff_vals = roi_real(:,r) - roi_imag(:,r);
    coh_d = mean(diff_vals) / std(diff_vals);

    results_parametric(r).roi_name = ROIs(r).name;
    results_parametric(r).t = stats.tstat;
    results_parametric(r).df = stats.df;
    results_parametric(r).p = p;
    results_parametric(r).cohens_d = coh_d;
    results_parametric(r).mean_real = mean(roi_real(:,r));
    results_parametric(r).mean_imag = mean(roi_imag(:,r));

    fprintf("  Real: %.2f ± %.2f\n", mean(roi_real(:,r)), std(roi_real(:,r)));
    fprintf("  Imag: %.2f ± %.2f\n", mean(roi_imag(:,r)), std(roi_imag(:,r)));
    fprintf("  t(%d)=%.3f, p=%.4f\n", stats.df, stats.tstat, p);
    fprintf("  Cohen's d = %.3f\n", coh_d);
end

%% ============================================================
%          NON-PARAMETRIC PERMUTATION TEST
% ============================================================
fprintf("\n==============================================\n");
fprintf("            PERMUTATION TESTS\n");
fprintf("==============================================\n");

results_permutation = struct();

for r = 1:length(ROIs)
    fprintf("\nROI %d: %s — %d permutations\n", r, ROIs(r).name, n_permutations);

    obs_diff = mean(roi_real(:,r) - roi_imag(:,r));
    null_diffs = zeros(n_permutations,1);

    % sign-flip permutation test (paired design)
    for p = 1:n_permutations
        flip = (rand(n_subjects,1) > 0.5) * 2 - 1;
        null_diffs(p) = mean((roi_real(:,r)-roi_imag(:,r)) .* flip);
    end

    p_perm = mean(abs(null_diffs) >= abs(obs_diff));

    results_permutation(r).roi_name = ROIs(r).name;
    results_permutation(r).observed_diff = obs_diff;
    results_permutation(r).p_perm = p_perm;

    fprintf("  Observed = %.3f dB, p_perm = %.4f\n", obs_diff, p_perm);
end

%% ============================================================
%         VISUALIZATIONS (Barplots, Null distributions)
% ============================================================
fprintf("\nGenerating figures...\n");

% ------- Barplot Comparison -------
fig1 = figure('Position',[100 100 1200 800]);

for r = 1:length(ROIs)
    subplot(2,2,r);

    mR = results_parametric(r).mean_real;
    mI = results_parametric(r).mean_imag;
    sR = std(roi_real(:,r))/sqrt(n_subjects);
    sI = std(roi_imag(:,r))/sqrt(n_subjects);

    bar([mR, mI]); hold on;
    errorbar(1:2,[mR,mI],[sR,sI],'k.');

    set(gca,'XTickLabel',{'Real','Imag'}); ylabel('Power (dB)');
    title(strrep(ROIs(r).name,'_',' '));
end

saveas(fig1, fullfile(out_path,'ROI_Comparison.png'));

% ------- Permutation Null Distributions -------
fig2 = figure('Position',[100 100 1200 800]);

for r = 1:length(ROIs)
    subplot(2,2,r);

    hist(results_permutation(r).null_diffs,50); hold on;
    xline(results_permutation(r).observed_diff,'r','LineWidth',2);

    title(sprintf('%s (p=%.4f)', ...
        strrep(ROIs(r).name,'_',' '), results_permutation(r).p_perm));
end

saveas(fig2, fullfile(out_path,'Permutation_Nulls.png'));

% ------- Effect sizes -------
fig3 = figure('Position',[100 100 800 500]);
dvals = [results_parametric.cohens_d];
barh(dvals); xlabel("Cohen's d");
set(gca,'YTick',1:length(ROIs),'YTickLabel',{ROIs.name});
saveas(fig3, fullfile(out_path,'Effect_Sizes.png'));

%% ============================================================
%               SAVE RESULTS
% ============================================================
summary_table = struct2table(results_parametric);
summary_table.p_perm = [results_permutation.p_perm]';
summary_table.p_bonf = summary_table.p * length(ROIs);
writetable(summary_table, fullfile(out_path,'ROI_Statistics.csv'));

fprintf("\nSaved all ROI statistics & figures to: %s\n", out_path);
fprintf("============================================================\n");

