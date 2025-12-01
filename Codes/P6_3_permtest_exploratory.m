%% ===============================================================
%  PHASE 6_3: CLUSTER-BASED PERMUTATION (STRICT + EXPLORATORY)
% ===============================================================

clear; clc;

%% ---- CONFIG ----
root = "./results/tfr_new_results/";
tfr_path = fullfile(root, "tfr_subjects");
out_path = fullfile(root, "permutation_stats_exploratory");
if ~exist(out_path,'dir'), mkdir(out_path); end

n_perm = 1000;                     % Number of permutations for FWER control
cluster_form_p = 0.05;            % Pixelwise threshold for forming clusters
                                  % (NOT the corrected threshold—used only for selecting suprathreshold pixels)

% Exploratory maps:
relaxed_p1 = 0.05; min_cluster1 = 20;   % Uncorrected threshold + minimum spatial extent
relaxed_p2 = 0.01; min_cluster2 = 10;

fprintf("\n==============================================\n");
fprintf("   CLUSTER PERMUTATION + EXPLORATORY MAPS\n");
fprintf("==============================================\n");

%% ---- LOAD SUBJECTS ----
all_files = dir(fullfile(tfr_path,"*.mat"));
names = {all_files.name};
tokens = regexp(names,"^(S\d+)_","tokens");

subjects = {};
for i=1:numel(tokens)
    if ~isempty(tokens{i}), subjects{end+1} = tokens{i}{1}{1}; end
end
subjects = unique(subjects);

fprintf("\nFound %d subjects: ", numel(subjects));

%% ---- COLLECT REAL & IMAGINED TFR FOR C3 ----
% We use C3 because motor imagery ERD/ERS effects are strongest here.
real_data = [];
imag_data = [];
freqs = []; times = [];

for s = 1:numel(subjects)
    subj = subjects{s};

    rfile = dir(fullfile(tfr_path, sprintf("%s_R*_executed_tfr.mat", subj)));
    ifile = dir(fullfile(tfr_path, sprintf("%s_R*_imagined_tfr.mat", subj)));

    if isempty(rfile) || isempty(ifile), continue; end

    R = load(fullfile(rfile(1).folder, rfile(1).name));
    I = load(fullfile(ifile(1).folder, ifile(1).name));

    % Ensure both have Morlet results
    if ~isfield(R.subj_tfr,'C3') || ~isfield(I.subj_tfr,'C3'), continue; end

    if isempty(freqs)
        freqs = R.subj_tfr.C3.freqs;     % all subjects share same freq/time axes
        times = R.subj_tfr.C3.times;
    end

    real_data(end+1,:,:) = R.subj_tfr.C3.power_db;   % (subjects × freq × time)
    imag_data(end+1,:,:) = I.subj_tfr.C3.power_db;

    fprintf(".");
end
fprintf("\n");

nS = size(real_data,1);
nF = size(real_data,2);
nT = size(real_data,3);

fprintf("Loaded %d subjects, %d freqs, %d times\n", nS, nF, nT);
if nS < 3, error("Need ≥3 subjects."); end

%% ---- OBSERVED T-MAP ----
% Classic paired t-test per pixel (freq × time). No multiple comparisons yet.
fprintf("\nComputing observed t-statistics...\n");

diff = real_data - imag_data;  % paired differences
mean_diff = squeeze(mean(diff, 1));
std_diff  = squeeze(std(diff, 0, 1));
se_diff   = std_diff / sqrt(nS);
observed_t = mean_diff ./ se_diff;

% Two-tailed p-values (not used for correction, but for exploratory maps)
df = nS - 1;
observed_p = 2 * tcdf(-abs(observed_t), df);

fprintf("Done.\n");

%% ===============================================================
%   STRICT CLUSTER-BASED PERMUTATION TEST (CONTROLS FWER)
% ===============================================================

% Pre-clustering threshold:
% Convert p<0.05 to a t-value. This is NOT corrected; used only to form clusters.
t_crit = tinv(1 - cluster_form_p/2, df);

fprintf("\nRunning %d permutations (strict FWER control)...\n", n_perm);

null_max_cluster = zeros(n_perm,1);

tic
for p = 1:n_perm
    if mod(p,50)==0
        fprintf("  %d/%d (%.1f%%)\n", p, n_perm, p/n_perm*100);
    end

    % Random sign flips implement the null hypothesis:
    % Real − Imag has no systematic direction across subjects.
    sign_flip = (rand(nS,1) > 0.5)*2 - 1;      % ±1 with equal probability
    sign_flip = reshape(sign_flip,[nS 1 1]);

    shuffled = diff .* sign_flip;              % permuted paired differences

    % Compute permuted t-map
    mean_s = squeeze(mean(shuffled,1));
    std_s  = squeeze(std(shuffled,0,1));
    se_s   = std_s / sqrt(nS);
    perm_t = mean_s ./ se_s;

    % Identify suprathreshold pixels
    suprathresh = abs(perm_t) > t_crit;

    % Measure cluster mass (sum of |t| inside each cluster)
    % This is the standard cluster statistic for FWER correction.
    if exist('bwlabel','file')
        [lbl, nCl] = bwlabel(suprathresh);
        if nCl > 0
            cluster_sums = zeros(nCl,1);
            for c = 1:nCl
                cluster_mask = (lbl == c);
                cluster_sums(c) = sum(abs(perm_t(cluster_mask)));
            end
            null_max_cluster(p) = max(cluster_sums);
        else
            null_max_cluster(p) = 0;
        end
    else
        % Manual clustering when Image Processing Toolbox is missing
        info = manual_cluster_sum(suprathresh, perm_t);
        null_max_cluster(p) = info.max_sum;
    end
end
toc

cluster_threshold = prctile(null_max_cluster, 95);  % FWER=0.05
fprintf("\nFWER cluster threshold = %.3f\n", cluster_threshold);

%% ---- APPLY STRICT CORRECTION ----
strict_sig = false(size(observed_t));

% Clusters in observed data that exceed permutation threshold
suprathresh_obs = abs(observed_t) > t_crit;

if exist('bwlabel','file')
    [lbl_obs, nCl_obs] = bwlabel(suprathresh_obs);
    for c = 1:nCl_obs
        mask = (lbl_obs == c);
        cluster_sum = sum(abs(observed_t(mask)));
        if cluster_sum > cluster_threshold
            strict_sig = strict_sig | mask;
        end
    end
else
    cl = manual_cluster_analysis(suprathresh_obs, observed_t, cluster_threshold);
    strict_sig = cl.sig_mask;
end

fprintf("Strict-corrected significant pixels: %d\n", sum(strict_sig(:)));

%% ===============================================================
%                EXPLORATORY (UNCORRECTED) MAPS
% ===============================================================

% These maps do NOT control FWER.
% They help visualize weak/early/late effects that don't survive correction.
fprintf("\nCreating exploratory maps...\n");

% Exploratory 1: moderate uncorrected threshold + cluster extent filter
t1 = tinv(1 - relaxed_p1/2, df);
mask1 = abs(observed_t) > t1;
mask1 = filter_clusters(mask1, min_cluster1);
fprintf("Exploratory1 (p<%.2f, cluster≥%d): %d pixels\n", relaxed_p1, min_cluster1, sum(mask1(:)));

% Exploratory 2: stricter uncorrected p + smaller minimum cluster size
t2 = tinv(1 - relaxed_p2/2, df);
mask2 = abs(observed_t) > t2;
mask2 = filter_clusters(mask2, min_cluster2);
fprintf("Exploratory2 (p<%.2f, cluster≥%d): %d pixels\n", relaxed_p2, min_cluster2, sum(mask2(:)));

%% ===============================================================
%                         SAVE RESULTS
% ===============================================================

results = struct();
results.observed_t = observed_t;
results.observed_p = observed_p;
results.strict_sig = strict_sig;
results.exploratory1 = mask1;
results.exploratory2 = mask2;
results.freqs = freqs;
results.times = times;
results.null_max_cluster = null_max_cluster;
results.cluster_threshold = cluster_threshold;
results.n_subjects = nS;

save(fullfile(out_path,"cluster_results.mat"), 'results');

%% ===============================================================
%                  VISUALIZATION FIGURES
% ===============================================================

fprintf("Creating visualizations...\n");

% ---------- Figure 1: strict vs exploratory ----------
fig1 = figure("Position",[50 50 1600 500]);

subplot(1,3,1)
imagesc(times,freqs,strict_sig);
axis xy; colorbar;
title(sprintf("Strict FWER corrected (%d px)", sum(strict_sig(:))));
colormap(gca,[1 1 1; 1 0 0]);

subplot(1,3,2)
imagesc(times,freqs,mask1);
axis xy; colorbar;
title(sprintf("Exploratory p<%.2f (≥%d px)", relaxed_p1, min_cluster1));
colormap(gca,[1 1 1; 0 0 1]);

subplot(1,3,3)
imagesc(times,freqs,mask2);
axis xy; colorbar;
title(sprintf("Exploratory p<%.2f (≥%d px)", relaxed_p2, min_cluster2));
colormap(gca,[1 1 1; 0 0.5 0]);

saveas(fig1, fullfile(out_path,"Exploratory_Cluster_Comparison.png"));

% ---------- Figure 2: T-statistics with overlays ----------
fig2 = figure("Position",[50 50 1600 500]);

subplot(1,3,1)
imagesc(times,freqs,observed_t);
axis xy; colorbar;
title("Observed t-statistics");
colormap(gca, jet);

subplot(1,3,2)
imagesc(times,freqs,observed_t); hold on;
contour(times,freqs,mask1,1,'k','LineWidth',1.5);
axis xy; colorbar;
title("Exploratory 1 Overlay");

subplot(1,3,3)
imagesc(times,freqs,observed_t); hold on;
contour(times,freqs,mask2,1,'w','LineWidth',1.5);
axis xy; colorbar;
title("Exploratory 2 Overlay");

saveas(fig2, fullfile(out_path,"T_Stats_With_Overlays.png"));

% ---------- Figure 3: Null distribution ----------
fig3 = figure("Position",[50 50 800 400]);
histogram(null_max_cluster,50);
hold on;
xline(cluster_threshold,'r--','LineWidth',2);
xlabel("Max cluster mass (null)");
ylabel("Count");
title("Permutation Null Distribution");
saveas(fig3, fullfile(out_path,"Null_Distribution.png"));

%% ===============================================================
%                            SUMMARY
% ===============================================================

fprintf("\n==============================================\n");
fprintf("              RESULTS SUMMARY\n");
fprintf("==============================================\n");
fprintf("Subjects: %d\n", nS);
fprintf("Strict corrected significant pixels: %d\n", sum(strict_sig(:)));
fprintf("Exploratory 1 pixels: %d\n", sum(mask1(:)));
fprintf("Exploratory 2 pixels: %d\n", sum(mask2(:)));
fprintf("Results saved in: %s\n", out_path);
fprintf("==============================================\n");

%% ———————————————— HELPER FUNCTIONS ————————————————

function mask_out = filter_clusters(mask_in, min_size)
% Removes clusters smaller than min_size (extent thresholding)
...
end

function info = manual_cluster_sum(mask, t_values)
% Computes cluster mass without Image Processing Toolbox
...
end

function info = manual_cluster_analysis(mask, t_values, threshold)
% Manual strict significance test without IPT
...
end

