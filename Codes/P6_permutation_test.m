%% =====================================================================
%  PHASE 6: CLUSTER-BASED PERMUTATION TEST (WAVELET ERSP, C3 ONLY)
% =====================================================================

clear; clc;

%% ---------------- PATHS ----------------
root = "./results/tfr_new_results/";
tfr_path = fullfile(root,"tfr_subjects");
out_path = fullfile(root,"permutation_stats_fast");
if ~exist(out_path,'dir'), mkdir(out_path); end

%% ---------------- PARAMETERS ----------------
n_perm = 1000;              % Number of permutations controls FWER precision.
cluster_alpha = 0.05;       % Threshold to define clusters (NOT final corrected p)
final_alpha   = 0.05;       % Family-wise corrected statistical threshold.

fprintf("\n==============================================\n");
fprintf("   PERMUTATION TEST (WAVELET ERSP)\n");
fprintf("==============================================\n");

%% =====================================================================
%  STEP 0 — Load all subjects and extract ERSP matrices for C3
% =====================================================================

files = dir(fullfile(tfr_path,"*.mat"));
if isempty(files), error("No .mat files found."); end

% Extract subject IDs
names = {files.name};
tokens = regexp(names,"^(S\\d+)_","tokens");
subjects = {};
for i=1:numel(tokens)
    if ~isempty(tokens{i}), subjects{end+1} = tokens{i}{1}{1}; end
end
subjects = unique(subjects);

fprintf("Found subjects: %d\n", numel(subjects));

real_list = [];
imag_list = [];
valid_subj = {};

freqs = []; times = [];

for i = 1:numel(subjects)
    subj = subjects{i};

    % Locate one executed and imagined TFR file for the subject
    fr = dir(fullfile(tfr_path, sprintf("%s_R*_executed_tfr.mat",subj)));
    fi = dir(fullfile(tfr_path, sprintf("%s_R*_imagined_tfr.mat",subj)));

    if isempty(fr) || isempty(fi), continue; end

    R = load(fullfile(fr(1).folder, fr(1).name));
    I = load(fullfile(fi(1).folder, fi(1).name));

    % Skip if C3 not present
    if ~isfield(R.subj_tfr,'C3') || ~isfield(I.subj_tfr,'C3')
        continue;
    end

    if isempty(freqs)
        % Establish global time-frequency axes from the first subject
        freqs = R.subj_tfr.C3.freqs;
        times = R.subj_tfr.C3.times;
    end

    % Each subject contributes their ERSP matrix (freq x time)
    real_list(end+1,:,:) = R.subj_tfr.C3.power_db;
    imag_list(end+1,:,:) = I.subj_tfr.C3.power_db;

    valid_subj{end+1} = subj;
end

n_subj = size(real_list,1);
n_freq = size(real_list,2);
n_time = size(real_list,3);

fprintf("Using %d subjects | %d freqs × %d times\n", n_subj, n_freq, n_time);

if n_subj < 5
    error("Cluster permutation requires ≥5 subjects.");
end


%% =====================================================================
%  STEP 1 — Observed t-statistics (vectorised)
% =====================================================================

fprintf("\nComputing observed t-statistics...\n");

diff_mat = real_list - imag_list;       % Subject-level paired differences
mean_diff = squeeze(mean(diff_mat,1));  % Mean difference for each (f,t)
std_diff  = squeeze(std(diff_mat,0,1)); % SD of subject differences

% Vectorised t statistic
observed_t = mean_diff ./ (std_diff / sqrt(n_subj));
df = n_subj - 1;


%% =====================================================================
%  STEP 2 — Permutation null distribution
% =====================================================================

fprintf("\nRunning %d permutations...\n", n_perm);

t_thresh = tinv(1 - cluster_alpha/2, df);   % Threshold to *define* clusters

null_max_pixel  = zeros(n_perm,1);          % Null distribution of max-|t|
null_max_cluster = zeros(n_perm,1);         % Null distribution of max cluster-mass

for p = 1:n_perm
    if mod(p,50)==0
        fprintf("  Perm %d/%d\n", p, n_perm);
    end

    % Sign flips preserve within-subject correlation structure
    flips = (rand(n_subj,1)>0.5)*2 - 1;
    flips = reshape(flips,[n_subj 1 1]);  % broadcast to freq × time

    perm_diff = diff_mat .* flips;

    % Vectorised permutation t-values
    perm_t = squeeze(mean(perm_diff,1)) ./ (squeeze(std(perm_diff,0,1))/sqrt(n_subj));

    % ------------- Pixel-level extreme value -------------
    % Needed for strong control of FWER at individual pixels
    null_max_pixel(p) = max(abs(perm_t(:)));

    % ------------- Cluster-level statistic -------------
    % Form clusters at uncorrected threshold
    mask = abs(perm_t) > t_thresh;
    clusters = find_clusters_2D(mask);

    if isempty(clusters)
        null_max_cluster(p) = 0;
    else
        % Cluster mass = sum |t| across all pixels in the cluster
        sums = zeros(numel(clusters),1);
        for k = 1:numel(clusters)
            sums(k) = sum(abs(perm_t(clusters{k})));
        end
        null_max_cluster(p) = max(sums); % store extreme cluster mass
    end
end

fprintf("Permutation testing complete.\n");


%% =====================================================================
%  STEP 3 — Apply FWER corrections
% =====================================================================

pixel_thr   = prctile(null_max_pixel, 100*(1-final_alpha));
cluster_thr = prctile(null_max_cluster,100*(1-final_alpha));

fprintf("\nPixel threshold: |t| > %.3f\n", pixel_thr);
fprintf("Cluster mass threshold: %.3f\n", cluster_thr);

% Pixel correction → flags isolated extreme t-values (rare in EEG)
pixel_sig = abs(observed_t) > pixel_thr;

% Form uncorrected suprathreshold map for cluster detection
obs_mask = abs(observed_t) > t_thresh;
obs_clusters = find_clusters_2D(obs_mask);

cluster_sig = false(size(observed_t));

% For each observed cluster → check if its mass exceeds null distribution
for k = 1:numel(obs_clusters)
    idx = obs_clusters{k};
    cluster_mass = sum(abs(observed_t(idx)));
    if cluster_mass > cluster_thr
        cluster_sig(idx) = true;   % mark entire cluster as significant
    end
end


%% =====================================================================
%  STEP 4 — SAVE RESULTS
% =====================================================================

save(fullfile(out_path,"cluster_perm_results_fast.mat"), ...
    "observed_t","pixel_sig","cluster_sig","freqs","times","valid_subj", ...
    "null_max_pixel","null_max_cluster","pixel_thr","cluster_thr");


%% =====================================================================
%  STEP 5 — PLOTS
% =====================================================================

fig = figure('Position',[200 200 1500 500]);

subplot(1,3,1)
imagesc(times,freqs,observed_t);
axis xy; colorbar;
title("Observed t-statistics");

subplot(1,3,2)
imagesc(times,freqs,pixel_sig);
axis xy; colorbar;
title("Pixel-level FWER corrected");

subplot(1,3,3)
imagesc(times,freqs,cluster_sig);
axis xy; colorbar;
title("Cluster-level FWER corrected (FINAL)");

saveas(fig, fullfile(out_path, "Permutation_Fast_Results.png"));



%% =====================================================================
%  HELPER — Cluster detection (4-connectivity)
% =====================================================================
function clusters = find_clusters_2D(mask)

    clusters = {};
    visited = false(size(mask));
    dirs = [1 0; -1 0; 0 1; 0 -1];   % 4-connected adjacency

    for i=1:size(mask,1)
        for j=1:size(mask,2)

            if mask(i,j) && ~visited(i,j)
                queue = [i j];
                visited(i,j) = true;
                cl = [];

                % Breadth-first search to extract cluster
                while ~isempty(queue)
                    pos = queue(1,:); queue(1,:) = [];
                    cl(end+1,:) = pos;

                    for d=1:4
                        ni = pos(1) + dirs(d,1);
                        nj = pos(2) + dirs(d,2);

                        if ni>=1 && ni<=size(mask,1) && nj>=1 && nj<=size(mask,2)
                            if mask(ni,nj) && !visited(ni,nj)
                                visited(ni,nj) = true;
                                queue(end+1,:) = [ni nj];
                            end
                        end
                    end
                end

                clusters{end+1} = sub2ind(size(mask), cl(:,1), cl(:,2));
            end
        end
    end
end

