%% ===============================================================
%  PHASE 2_2 : GRAND-AVERAGE ERP PLOT (Across Subjects) 
%  ===============================================================
clear; clc;

erp_folder = './results/erp/';
fig_folder = './results/figures/';
if ~exist(fig_folder,'dir'), mkdir(fig_folder); end

erp_files = dir(fullfile(erp_folder,'S*_ERP.mat'));

% --- Define common epoch range and time base ---
Fs_target = 160;          % target sampling rate for interpolation
epoch_range = [-1, 4];    % seconds
t_common = linspace(epoch_range(1), epoch_range(2), diff(epoch_range)*Fs_target); 
% Creates a common time axis for all subjects (ensures same length)

channels_to_plot = [9 11 13]; % C3, Cz, C4 (motor cortex)

% Matrices storing interpolated ERPs for each subject
all_real = [];
all_imag = [];

fprintf('Loaded %d ERP files for group grand average.\n', length(erp_files));

for i = 1:length(erp_files)
    fname = fullfile(erp_folder, erp_files(i).name);
    load(fname, 'ERP_struct');   % load the subject-specific ERP struct

    % --- Locate fields ---
    if isfield(ERP_struct,'Real_LeftRight_Fist')
        real = ERP_struct.Real_LeftRight_Fist;
    elseif isfield(ERP_struct,'Real_Hands')
        real = ERP_struct.Real_Hands;
    else
        fprintf('Skipping %s (no Real task)\n', fname); 
        continue;
    end

    if isfield(ERP_struct,'Imagined_LeftRight_Fist')
        imag = ERP_struct.Imagined_LeftRight_Fist;
    elseif isfield(ERP_struct,'Imagined_Hands')
        imag = ERP_struct.Imagined_Hands;
    else
        fprintf('Skipping %s (no Imagined task)\n', fname); 
        continue;
    end

    % --- Average T1 & T2 and subtract T0 baseline condition ---
    % T1_mean and T2_mean are (channels × samples); cat(3,...) stacks them
    real_avg = mean(cat(3, real.T1_mean, real.T2_mean), 3, 'omitnan');
    imag_avg = mean(cat(3, imag.T1_mean, imag.T2_mean), 3, 'omitnan');

    % Difference waves: (movement - rest)
    real_diff = real_avg - real.T0_mean;
    imag_diff = imag_avg - imag.T0_mean;

    % --- Construct subject-specific time axis based on its sample count ---
    subj_len = size(real_diff,2);                    % number of samples
    subj_Fs = subj_len / diff(epoch_range);          % infer sampling rate
    t_subj = linspace(epoch_range(1), epoch_range(2), subj_len);
    valid_channels = channels_to_plot(channels_to_plot <= size(real_diff,1));

    % --- Average ERP across C3–Cz–C4 for the subject ---
    real_mean = mean(real_diff(valid_channels,:), 1, 'omitnan');
    imag_mean = mean(imag_diff(valid_channels,:), 1, 'omitnan');

    % --- Interpolate subject ERP onto the common time axis ---
    real_interp = interp1(t_subj, real_mean, t_common, 'linear', 'extrap');
    imag_interp = interp1(t_subj, imag_mean, t_common, 'linear', 'extrap');

    % Store into group matrices
    all_real(i,:) = real_interp;
    all_imag(i,:) = imag_interp;
end

% ---- Compute group average and SEM ----
grand_real = mean(all_real,1,'omitnan');
grand_imag = mean(all_imag,1,'omitnan');

% SEM = std / sqrt(N)
sem_real = std(all_real,0,1,'omitnan') / sqrt(size(all_real,1));
sem_imag = std(all_imag,0,1,'omitnan') / sqrt(size(all_imag,1));

% ---- Plot ----
fig = figure('Position',[200 200 900 500]);
hold on;

% Draw shaded error bands
fill([t_common fliplr(t_common)], ...
     [grand_real+sem_real fliplr(grand_real-sem_real)], ...
     [1 0.6 0.6], 'EdgeColor','none');        % red-ish for real

fill([t_common fliplr(t_common)], ...
     [grand_imag+sem_imag fliplr(grand_imag-sem_imag)], ...
     [0.6 0.6 1], 'EdgeColor','none');        % blue-ish for imagined

% Plot mean ERPs
plot(t_common, grand_real, 'r', 'LineWidth',2);
plot(t_common, grand_imag, 'b', 'LineWidth',2);

xlim([-0.5 2]);
xlabel('Time (s)'); 
ylabel('Amplitude (µV)');
title('Grand-Average ERP (Real vs Imagined) - C3–Cz–C4');
legend({'Real ± SEM','Imagined ± SEM'}, 'Location','best');
grid on;

saveas(fig, fullfile(fig_folder,'Group_GrandAverage_ERP.png'));
fprintf('Saved grand-average ERP plot to results/figures/Group_GrandAverage_ERP.png\n');

