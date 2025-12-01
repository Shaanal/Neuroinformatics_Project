%% =========================================================================
%   PHASE 5_2: GROUP MOTOR IMAGERY ANALYSIS (20 subjects)
% =========================================================================

clear; clc; close all;
[ALLEEG, EEG, CURRENTSET] = eeglab;

%% =============================== CONFIG ==================================
processed_root = './processed/';
results_root   = './group_results/';
if ~exist(results_root,'dir'), mkdir(results_root); end

N_subjects = 20;                     % total subjects to analyze

% Run numbers corresponding to executed vs. imagined MI
executed_runs = [3 7 11];
imagined_runs = [4 8 12];

motor_tokens = {'C3','Cz','C4'};     % auto-detection tokens for motor electrodes
freqs = 4:40;                        % frequency range for ERSP
cycles = [3 0.5];                    % Morlet wavelet parameters
baseline_ms = [-1000 0];            % baseline for ERSP / Hilbert
task_window = [300 1200];           % ms window for ERD/ERS extraction

% Frequency bands for Hilbert method
bands = struct;
bands.theta = [4 7];
bands.alpha = [8 12];
bands.beta  = [13 30];

%% ======================== SUBJECT SELECTION ==============================
S = dir(fullfile(processed_root,'S*'));
S = S([S.isdir]);
subjects = {S.name};

rng(42);                                  % reproducible random selection
subjects = subjects(randperm(length(subjects)));
subjects = subjects(1:min(N_subjects,length(subjects)));

fprintf("\n=== Selected %d subjects ===\n", length(subjects));
disp(subjects');


%% ======================= STORAGE STRUCTURES ===============================
GA = struct;  % stores group-level wavelet + hilbert results
conditions = {'executed','imagined'};

% Create fields for each condition and each motor channel
for cond = conditions
    c = cond{1};
    GA.(c).wave = struct;
    GA.(c).hil  = struct;

    for m = 1:length(motor_tokens)
        tk = motor_tokens{m};

        % Wavelet ERSP + ITPC
        GA.(c).wave.(tk).pow = [];   % F × T × subjects
        GA.(c).wave.(tk).itpc = [];
        GA.(c).wave.(tk).times = [];
        GA.(c).wave.(tk).n = 0;

        % Hilbert-band envelopes
        GA.(c).hil.(tk).theta = [];
        GA.(c).hil.(tk).alpha = [];
        GA.(c).hil.(tk).beta  = [];
        GA.(c).hil.(tk).times = [];
        GA.(c).hil.(tk).n = 0;
    end
end

group_table = [];   % stores subject-level ERD/ERS features


%% =========================================================================
%                         SUBJECT LOOP
% =========================================================================
for si = 1:length(subjects)

    subj = subjects{si};
    subj_path = fullfile(processed_root, subj);

    fprintf("\n=================================================\n");
    fprintf("Processing Subject: %s\n", subj);
    fprintf("=================================================\n");

    % Stores subject-averaged ERSP/Hilbert for all runs
    subj_exec_wave = {};
    subj_imag_wave = {};
    subj_exec_hil  = {};
    subj_imag_hil  = {};

    %% ================= LOAD ALL MI RUNS =====================
    run_files = dir(fullfile(subj_path, '*_processed.set'));
    if isempty(run_files)
        fprintf("No processed files. Skipping.\n");
        continue;
    end

    % Per-run containers
    exec_wave_runs = {}; exec_hil_runs = {};
    imag_wave_runs = {}; imag_hil_runs = {};

    for r = 1:length(run_files)

        fname = run_files(r).name;

        % Extract run index (Rxx)
        tok = regexp(fname,'R(\d+)','tokens','once');
        if isempty(tok), continue; end
        run_num = str2double(tok{1});

        % Map run to executed/imagined condition
        if ismember(run_num, executed_runs)
            cond = 'executed';
        elseif ismember(run_num, imagined_runs)
            cond = 'imagined';
        else
            continue; % skip irrelevant runs
        end

        fprintf("   Loading %s (%s)\n", fname, cond);
        EEG = pop_loadset('filename',fname,'filepath',subj_path);

        % Skip if too few trials for reliable TF/Hilbert analysis
        if EEG.trials < 5
            fprintf("Not enough trials\n");
            continue;
        end

        %% ========== AUTO-DETECT MOTOR CHANNELS ==========
        labels = {EEG.chanlocs.labels};
        chan_idx = nan(1,length(motor_tokens));

        % Allow partial-match detection of C3/Cz/C4
        for m = 1:length(motor_tokens)
            tk = motor_tokens{m};
            idx = find(contains(labels, tk,'IgnoreCase',true));
            if isempty(idx)
                idx = find(strcmpi(labels, tk));
            end
            if ~isempty(idx)
                chan_idx(m) = idx(1);
            end
        end

        if all(isnan(chan_idx))
            fprintf("None of C3/Cz/C4 found. Skipping run.\n");
            continue;
        end

        %% ----------------- WAVELET (Morlet) ------------------
        for m = 1:length(motor_tokens)

            if isnan(chan_idx(m)), continue; end

            ch = chan_idx(m);
            data = squeeze(EEG.data(ch,:,:));  % (time × trials)

            % ERSP (dB) from newtimef
            try
                [ersp,~,~,t_ms,f] = newtimef( ...
                    data, EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, ...
                    cycles, 'freqs',freqs, 'plotersp','off','plotitc','off');
            catch
                fprintf("Wavelet failed\n");
                continue;
            end

            W = struct();
            W.pow   = ersp;
            W.freqs = f;
            W.times = t_ms;

            % ----------------- ITPC ---------------------------
            try
                [~,itpc] = newtimef(data,EEG.pnts,[EEG.xmin EEG.xmax]*1000,...
                    EEG.srate, cycles, 'freqs',freqs, ...
                    'plotersp','off','plotitc','on');
            catch
                itpc = nan(size(ersp));
            end

            W.itpc = itpc;

            % Append run
            if strcmp(cond,'executed')
                exec_wave_runs{end+1,m} = W;
            else
                imag_wave_runs{end+1,m} = W;
            end
        end

        %% ----------------- HILBERT ---------------------------
        [~,nt,ntr] = size(EEG.data);        % samples × channels × trials
        t_sec = EEG.times/1000;             % time axis in seconds

        % Logical index for baseline sample points
        bidx = t_sec >= baseline_ms(1)/1000 & t_sec <= baseline_ms(2)/1000;

        for m = 1:length(motor_tokens)

            if isnan(chan_idx(m)), continue; end
            ch = chan_idx(m);

            % Extract data (samples × trials)
            sig = squeeze(EEG.data(ch,:,:));

            H = struct();

            % Process each frequency band separately
            for fn = fieldnames(bands)'
                bandname = fn{1};
                fr = bands.(bandname);

                % IIR bandpass per trial (Hilbert requires single-band analytic signal)
                bp = designfilt('bandpassiir','FilterOrder',4,...
                    'HalfPowerFrequency1',fr[1],...
                    'HalfPowerFrequency2',fr[2],...
                    'SampleRate',EEG.srate);

                amp = zeros(nt,ntr);

                for tr=1:ntr
                    x = filtfilt(bp,double(sig(:,tr))); % zero-phase filtering
                    A = abs(hilbert(x));               % analytic amplitude
                    amp(:,tr)=A;
                end

                % Baseline normalization in log(dB)
                base = mean(amp(bidx,:),1);
                amp_db = 10*log10(bsxfun(@rdivide,amp,base));

                H.(bandname) = mean(amp_db,2); % avg across trials
            end

            H.times = t_sec;

            if strcmp(cond,'executed')
                exec_hil_runs{end+1,m} = H;
            else
                imag_hil_runs{end+1,m} = H;
            end
        end

    end % run loop

    %% =================== MERGE RUNS INTO SUBJECT MEAN ====================
    if isempty(exec_wave_runs) || isempty(imag_wave_runs)
        fprintf("Subject %s missing Real/Imag → skipping\n",subj);
        continue;
    end

    % subject-level result row
    subj_row = struct;
    subj_row.Subject = {subj};

    for m = 1:length(motor_tokens)
        tk = motor_tokens{m};

        %% Merge WAVELET ERSP
        exec_runs = exec_wave_runs(:,m);
        imag_runs = imag_wave_runs(:,m);

        exec_runs = exec_runs(~cellfun(@isempty,exec_runs));
        imag_runs = imag_runs(~cellfun(@isempty,imag_runs));

        if isempty(exec_runs) || isempty(imag_runs)
            continue;
        end

        % Ensure all runs match in size (F × T)
        ref = exec_runs{1}.pow;
        F = size(ref,1); T = size(ref,2);

        EX = zeros(F,T,length(exec_runs));
        IM = zeros(F,T,length(imag_runs));

        for k=1:length(exec_runs), EX(:,:,k)=exec_runs{k}.pow; end
        for k=1:length(imag_runs), IM(:,:,k)=imag_runs{k}.pow; end

        exec_subj = mean(EX,3);
        imag_subj = mean(IM,3);

        % Append ERSP to group arrays
        GA.executed.wave.(tk).pow(:,:,end+1)=exec_subj;
        GA.imagined.wave.(tk).pow(:,:,end+1)=imag_subj;

        % Compute subject-mean ITPC
        GA.executed.wave.(tk).itpc(:,:,end+1)=mean(cat(3,exec_runs{:}.itpc),3);
        GA.imagined.wave.(tk).itpc(:,:,end+1)=mean(cat(3,imag_runs{:}.itpc),3);

        GA.executed.wave.(tk).times = exec_runs{1}.times;
        GA.imagined.wave.(tk).times = imag_runs{1}.times;

        GA.executed.wave.(tk).n = GA.executed.wave.(tk).n+1;
        GA.imagined.wave.(tk).n = GA.imagined.wave.(tk).n+1;

        %% ================== MERGE HILBERT =====================
        execH = exec_hil_runs(:,m);
        imagH = imag_hil_runs(:,m);

        execH = execH(~cellfun(@isempty,execH));
        imagH = imagH(~cellfun(@isempty,imagH));

        if isempty(execH) || isempty(imagH)
            continue;
        end

        t = execH{1}.times;

        % Construct subject-level matrices for theta/alpha/beta
        E = []; I = [];

        for fn = fieldnames(bands)'
            b = fn{1};
            E = [E mean(cat(2,execH{:}.(b)),2)]; 
            I = [I mean(cat(2,imagH{:}.(b)),2)]; 
        end

        % Append per-band time-courses to group structure
        GA.executed.hil.(tk).theta(:,end+1)=E(:,1);
        GA.executed.hil.(tk).alpha(:,end+1)=E(:,2);
        GA.executed.hil.(tk).beta(:,end+1) =E(:,3);

        GA.imagined.hil.(tk).theta(:,end+1)=I(:,1);
        GA.imagined.hil.(tk).alpha(:,end+1)=I(:,2);
        GA.imagined.hil.(tk).beta(:,end+1) =I(:,3);

        GA.executed.hil.(tk).times = t;
        GA.imagined.hil.(tk).times = t;

        GA.executed.hil.(tk).n = GA.executed.hil.(tk).n+1;
        GA.imagined.hil.(tk).n = GA.imagined.hil.(tk).n+1;

        %% Extract subject-level ERD/ERS features (from Hilbert alpha)
        tidx = t*1000>=task_window(1) & t*1000<=task_window(2);

        subj_row.(tk+"_alpha_exec") = mean(E(tidx,2),'omitnan');
        subj_row.(tk+"_alpha_imag") = mean(I(tidx,2),'omitnan');

    end % channels

    % Append subject row to table
    group_table = [group_table; struct2table(subj_row)];

end % subject loop

%% =========================================================================
%                     GROUP STATISTICS + PLOTS
% =========================================================================

stats_dir = fullfile(results_root,'stats');
if ~exist(stats_dir,'dir'), mkdir(stats_dir); end

disp(" ");
disp("=== GROUP RESULTS ===");
writetable(group_table, fullfile(stats_dir,'group_features.csv'));

% Example: stats for C3 alpha ERD (Real vs Imag)
if isfield(group_table,'C3_alpha_exec') && isfield(group_table,'C3_alpha_imag')
    exec = group_table.C3_alpha_exec;
    imag = group_table.C3_alpha_imag;

    [~,p,~,s] = ttest(exec,imag);

    fprintf("\nC3 Alpha ERD:\n");
    fprintf("Exec = %.2f ± %.2f\n", mean(exec), std(exec));
    fprintf("Imag = %.2f ± %.2f\n", mean(imag), std(imag));
    fprintf("t(%d) = %.2f, p = %.4f\n", s.df, s.tstat, p);

    figure('Color','white');
    bar([mean(exec) mean(imag)]);
    hold on;
    errorbar(1:2,[mean(exec) mean(imag)],...
             [std(exec)/sqrt(length(exec)) std(imag)/sqrt(length(imag))],'.k');
    set(gca,'XTickLabel',{'Real','Imag'});
    title(sprintf("C3 Alpha ERD (p=%.4f)",p));
    ylabel('ERD (dB)');
    saveas(gcf,fullfile(stats_dir,'C3_alpha_ERD_barplot.png'));
    close;
end

fprintf("\nResults saved to: %s\n",results_root);

