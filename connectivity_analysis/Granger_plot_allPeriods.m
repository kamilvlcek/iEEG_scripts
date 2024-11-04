%%% aggregates Granger data for all task periods and conditions ('same' and 'diff' trials) from all patients
%%% computes the mean Granger values across channel pairs and subjects
%%% generates plots of the average Granger spectra for both directions
%%% highlights the maximal difference between the two directions 
% --------------------
% by Sofiia Moraresku
% September 2024
% --------------------

% Set up paths and variables 
outputPath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\Granger_Sept2024';

% Define the periods and corresponding data files
% periods = {
%     struct('name', 'encoding', 'Granger_data_s', 'Granger_VTC-IPL_1.9s encoding_same_trials_2024-09.mat', 'Granger_data_d', 'Granger_VTC-IPL_1.9s encoding_diff_trials_2024-09.mat'),
%     struct('name', 'delay_first', 'Granger_data_s', 'Granger_VTC-IPL_first 1.9s delay_same_trials_2024-09.mat', 'Granger_data_d', 'Granger_VTC-IPL_first 1.9s delay_diff_trials_2024-09.mat'),
%     struct('name', 'delay_last', 'Granger_data_s', 'Granger_VTC-IPL_last 1.9s delay_same_trials_2024-09.mat', 'Granger_data_d', 'Granger_VTC-IPL_last 1.9s delay_diff_trials_2024-09.mat'),
%     struct('name', 'baseline', 'Granger_data_s', 'Granger_VTC-IPL_1.9s baseline_same_trials_PLV_sign_pairs_2024-10.mat', 'Granger_data_d', 'Granger_VTC-IPL_1.9s baseline_diff_trials_PLV_sign_pairs_2024-10.mat')
% };

periods = {
    struct('name', 'encoding', 'Granger_data_s', 'Granger_Hip-IPL_1.9s encoding_same_trials_2024-09.mat', 'Granger_data_d', 'Granger_Hip-IPL_1.9s encoding_diff_trials_2024-09.mat'),
    struct('name', 'delay_first', 'Granger_data_s', 'Granger_Hip-IPL_first 1.9s delay_same_trials_2024-09.mat', 'Granger_data_d', 'Granger_Hip-IPL_first 1.9s delay_diff_trials_2024-09.mat'),
    struct('name', 'delay_last', 'Granger_data_s', 'Granger_Hip-IPL_last 1.9s delay_same_trials_2024-09.mat', 'Granger_data_d', 'Granger_Hip-IPL_last 1.9s delay_diff_trials_2024-09.mat'),
    struct('name', 'baseline', 'Granger_data_s', 'Granger_Hip-IPL_1.9s baseline_same_trials_PLV_sign_pairs_2024-10.mat', 'Granger_data_d', 'Granger_Hip-IPL_1.9s baseline_diff_trials_PLV_sign_pairs_2024-10.mat')
};

Gr_folder = 'Granger_permut_stat';
setup = setup_memact(1);
basedir = setup.basedir;
subfolder = setup.subfolder;

% Get the list of patients
[pacienti] = pacienti_memact(); % Set up pacienti_memact by selecting which patients to analyze

%% Initialize variables to store Granger values across all conditions and periods
grangerCombined = [];
freq_combined = [];

% Loop over each period and condition to collect Granger values
for i = 1:length(periods)
    period = periods{i};
    Granger_data_s = period.Granger_data_s;
    Granger_data_d = period.Granger_data_d;
    
    % Extract ROI info from the file name
    ROI1 = regexp(Granger_data_s, '(?<=Granger_)[^-]+', 'match', 'once');
    ROI2 = regexp(Granger_data_s, '(?<=-)[^-_]+', 'match', 'once');
    
    % Conditions
    conditions = {'same', 'diff'};
    
    for c = 1:length(conditions)
        condition = conditions{c};
        granger_mean_all = []; % Store mean Granger across subjects for this condition and period
        
        ipat = 1;
        
        for p = 1:numel(pacienti)
            if ~pacienti(p).todo
                continue; % Skip if the subject is not marked 'todo'
            end
            
            % Get patient folder
            patientFolder = pacienti(p).folder;
            filePath = [basedir patientFolder '\' subfolder '\' Gr_folder '\'];
            
            % Load Granger data
            if strcmp(condition, 'same')
                filename = Granger_data_s;
            else
                filename = Granger_data_d;
            end
            fullFileName = fullfile(filePath, filename);
            if isfile(fullFileName)
                data = load(fullFileName);
                grangerPeriod = data.grangerPeriod;
                grangerspctrm = grangerPeriod.grangerspctrm; % [N_rows x N_freqs]
                freq = grangerPeriod.freq;
                
                % Number of channel pairs
                N_rows = size(grangerspctrm, 1);
                N_pairs = N_rows / 2;
                N_freqs = length(freq);
                
                % Reshape to [N_pairs x 2 x N_freqs]
                grangerspctrm_reshaped = reshape(grangerspctrm, [2, N_pairs, N_freqs]);
                grangerspctrm_reshaped = permute(grangerspctrm_reshaped, [2,1,3]); % [N_pairs x 2 x N_freqs]
                
                % Compute mean over channel pairs for both directions
                mean_granger = squeeze(mean(grangerspctrm_reshaped, 1)); % [2 x N_freqs]
                
                % Store the results for each subject separately for both directions
                granger_mean_all(:,:,ipat) = mean_granger; % [2 x N_freqs x N_subjects]
                ipat = ipat + 1;
                
            else
                warning('File %s not found for patient %s', filename, patientFolder);
            end
        end
        
        % Store Granger data from this period/condition in a combined matrix
%         grangerCombined.(condition).mean_granger_all(:,:,i) = mean(granger_mean_all, 1); % [2 x N_freqs]
        grangerCombined.(condition).mean_granger_all(:,:,:,i) = granger_mean_all; % [2 x N_freqs x N_subjects x N_periods]
        freq_combined = freq; % Store frequency vector
    end
end

%% Combine across all periods and conditions
combined_granger_mean = [];
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Mean across all periods for each subject and direction
    granger_combined_mean_periods = squeeze(mean(grangerCombined.(condition).mean_granger_all, 4)); % Average over periods [2 x N_freqs x N_subjects]
    
    % Now average across subjects for each direction separately
    granger_combined_mean_condition = squeeze(mean(granger_combined_mean_periods, 3)); % [2 x N_freqs] (2 directions)
    
    % Add this condition's data to the combined matrix
    combined_granger_mean(:,:,c) = granger_combined_mean_condition; % [2 x N_freqs x N_conditions]
end

% Now average across conditions
final_granger_mean = squeeze(mean(combined_granger_mean, 3)); % [2 x N_freqs] (2 directions)

%% Plotting the average Granger spectra for both directions across all conditions and periods
figure('Name', 'Average Granger Spectrum Across All Conditions and Periods');
hold on;
plot(freq_combined, final_granger_mean(1,:), 'b', 'LineWidth', 2, 'DisplayName', 'Direction 1 (IPL -> VTC)');
plot(freq_combined, final_granger_mean(2,:), 'r', 'LineWidth', 2, 'DisplayName', 'Direction 2 (VTC -> IPL)');
xlabel('Frequency (Hz)');
ylabel('Mean Granger causality');
title('Average Granger Spectrum Across All Conditions and Periods');
legend;
grid on;

%% Plotting the difference between the two directions
difference_granger = final_granger_mean(2,:) - final_granger_mean(1,:); % VTC->IPL  - IPL->VTC; Hip->IPL  - IPL->Hip
[max_diff, idx_max_diff] = max(abs(difference_granger));
freq_max_diff = freq_combined(idx_max_diff);

figure('Name', 'Difference Between Directions');
plot(freq_combined, difference_granger, 'k', 'LineWidth', 2);
hold on;
plot(freq_max_diff, difference_granger(idx_max_diff), 'ro', 'MarkerSize', 8);
hold on;
% Add a dashed line at y=0 for reference
line(get(gca, 'xlim'), [0 0], 'Color', 'black', 'LineStyle', '--');
xlim([1 20])
xticks(1:1:20)
xlabel('Frequency (Hz)');
ylabel('Difference in Granger causality');
% title('Difference Between Directions (VTC->IPL) Across All Conditions and Periods');
title('Difference Between Directions (Hip->IPL) Across All Conditions and Periods');
legend(sprintf('Max Difference at %.2f Hz', freq_max_diff));
% grid on;



