% --------------------
% a script to run a Linear mixed effects model (LMEM) on Granger net flow for 2 conditions (same and diff) and 4 time periods - fixed effects, 
% while accounting for random effects - ch pair and patientID
% by Sofiia Moraresku
% September 2024
% --------------------

% Set up paths and variables 
outputPath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\Granger_Sept2024';

% Granger data for different conditions and periods
% pairs significant in PLV analysis
Granger_data_s_encoding = 'Granger_VTC-IPL_1.9s encoding_same_trials_2024-09.mat';
Granger_data_d_encoding = 'Granger_VTC-IPL_1.9s encoding_diff_trials_2024-09.mat';
Granger_data_s_delay_first = 'Granger_VTC-IPL_first 1.9s delay_same_trials_2024-09.mat';
Granger_data_d_delay_first = 'Granger_VTC-IPL_first 1.9s delay_diff_trials_2024-09.mat';
Granger_data_s_delay_last = 'Granger_VTC-IPL_last 1.9s delay_same_trials_2024-09.mat';
Granger_data_d_delay_last = 'Granger_VTC-IPL_last 1.9s delay_diff_trials_2024-09.mat';
Granger_data_s_bs = 'Granger_VTC-IPL_1.9s baseline_same_trials_PLV_sign_pairs_2024-10.mat';
Granger_data_d_bs = 'Granger_VTC-IPL_1.9s baseline_diff_trials_PLV_sign_pairs_2024-10.mat';

% Granger_data_s_encoding =  'Granger_Hip-IPL_1.9s encoding_same_trials_2024-09.mat';
% Granger_data_d_encoding =  'Granger_Hip-IPL_1.9s encoding_diff_trials_2024-09.mat';
% Granger_data_s_delay_first =  'Granger_Hip-IPL_first 1.9s delay_same_trials_2024-09.mat';
% Granger_data_d_delay_first =  'Granger_Hip-IPL_first 1.9s delay_diff_trials_2024-09.mat';
% Granger_data_s_delay_last =  'Granger_Hip-IPL_last 1.9s delay_same_trials_2024-09.mat';
% Granger_data_d_delay_last =  'Granger_Hip-IPL_last 1.9s delay_diff_trials_2024-09.mat';
% Granger_data_s_bs = 'Granger_Hip-IPL_1.9s baseline_same_trials_PLV_sign_pairs_2024-10.mat';
% Granger_data_d_bs = 'Granger_Hip-IPL_1.9s baseline_diff_trials_PLV_sign_pairs_2024-10.mat';

Gr_folder = 'Granger_permut_stat';
setup = setup_memact(1);
basedir = setup.basedir;
subfolder = setup.subfolder;

% Extract ROI info from the file name
ROI1 = regexp(Granger_data_s_encoding, '(?<=Granger_)[^-]+', 'match', 'once');
ROI2 = regexp(Granger_data_s_encoding, '(?<=-)[^-_]+', 'match', 'once');

% Get the list of patients
[pacienti] = pacienti_memact();

% Conditions and Periods
conditions = {'same', 'diff'};
periods = {'baseline', 'encoding', 'delay1', 'delay2'};
Granger_files_s = {Granger_data_s_bs, Granger_data_s_encoding, Granger_data_s_delay_first, Granger_data_s_delay_last};
Granger_files_d = {Granger_data_d_bs, Granger_data_d_encoding, Granger_data_d_delay_first, Granger_data_d_delay_last};

% Define frequency band to analyze
% freq_band = [9 11];
% freq_band = [2 4];
freq_band = [3 5];

%% Initialize a structure to store data tables for each condition
T_conditions = struct();

% Loop over conditions
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Initialize variables for saving results for this condition
    MeanNetGranger = []; % ?Granger values
    Participant = {};
    ChannelPair = {};
    Period = {};
    
    % Loop over periods
    for p_idx = 1:length(periods)
        period = periods{p_idx};
        
        % Loop over subjects
        for p = 1:numel(pacienti)
            if ~pacienti(p).todo
                continue;
            end
            
            % Get patient folder name
            patientFolder = pacienti(p).folder;
            % Construct the file path
            filePath = [basedir patientFolder '\' subfolder '\' Gr_folder '\'];
            
            % Determine file name based on condition and period
            if strcmp(condition, 'same')
                filename = Granger_files_s{p_idx};
            else
                filename = Granger_files_d{p_idx};
            end
            fullFileName = fullfile(filePath, filename);
            
            if isfile(fullFileName)
                % Load data
                data = load(fullFileName);
                grangerPeriod = data.grangerPeriod;
                grangerspctrm = grangerPeriod.grangerspctrm;
                labelcmb = grangerPeriod.labelcmb;
                freq = grangerPeriod.freq;
                freq_idx = find(freq >= freq_band(1) & freq <= freq_band(2));
                N_rows = size(grangerspctrm, 1);
                N_pairs = N_rows / 2; % Number of channel pairs
                N_freqs = length(freq);
                % Reshape to [N_pairs x 2 x N_freqs]
                grangerspctrm_reshaped = reshape(grangerspctrm, [2, N_pairs, N_freqs]);
                grangerspctrm_reshaped = permute(grangerspctrm_reshaped, [2, 1, 3]); % [N_pairs x 2 x N_freqs]
                mean_granger_freqband = mean(grangerspctrm_reshaped(:,:,freq_idx), 3); % [N_pairs x 2]
                
                for chnp = 1:N_pairs
                    chlabel1 = split(labelcmb{2*(chnp-1)+1,1}, '[');
                    chlabel2 = split(labelcmb{2*(chnp-1)+2,1}, '[');
                    pair_label = sprintf('%s-%s', chlabel1{1}, chlabel2{1});
                    
                    ichn_pair = data.ROI_chanpairs(chnp, :); % number of channel in chn pair (N pairs x 2)
                    % Determine the direction order based on ROI labels
                    if strcmp(data.dataPeriod.channelInfo(ichn_pair(2)).ROI, ROI1) % first direction in granger spectrum
                        % Direction 1: ROI1 -> ROI2
                        % Direction 2: ROI2 -> ROI1
                        net_granger = mean_granger_freqband(chnp,1) - mean_granger_freqband(chnp,2);
                    elseif strcmp(data.dataPeriod.channelInfo(ichn_pair(2)).ROI, ROI2)
                        % Direction 1: ROI2 -> ROI1
                        % Direction 2: ROI1 -> ROI2
                        net_granger = mean_granger_freqband(chnp,2) - mean_granger_freqband(chnp,1);
                    end
                    
                    MeanNetGranger(end+1) = net_granger;
                    Participant{end+1} = patientFolder;
                    ChannelPair{end+1} = pair_label;
                    Period{end+1} = period;
                end
            else
                warning('File %s not found for patient %s', filename, patientFolder);
            end
        end
    end
    
    % Create table for this condition and store it in the structure
    T = table(MeanNetGranger', Participant', ChannelPair', Period', ...
        'VariableNames', {'MeanNetGranger', 'Participant', 'ChannelPair', 'Period'});
    T.Participant = categorical(T.Participant);
    T.ChannelPair = categorical(T.ChannelPair);
    % Ensure Period is a categorical variable with the correct levels
    T.Period = categorical(T.Period, {'baseline', 'encoding', 'delay1', 'delay2'}, 'Ordinal', false);
    % Store the table in the structure
    T_conditions.(condition) = T;
end


%% fit LMEM for each condition separately
% Initialize variables for saving results
results_summary = {};  % To store results for saving in Excel

for ic = 1:length(conditions)
    condition = conditions{ic};
    T = T_conditions.(condition);
    
    % Fit linear mixed model with Period as fixed effect for this condition
    lme = fitlme(T, 'MeanNetGranger ~ -1+Period + (1|Participant) + (1|Participant:ChannelPair)', 'DummyVarCoding', 'full'); % without intercept, each period is tested against zero
    
    % Display the results
    fprintf('\nLinear Mixed Model Results for Condition: %s\n', condition);
    disp(lme);
    
    % Define the number of comparisons (number of models)
    numComparisons = 2;
    
    % original alpha level
    alpha = 0.05;
    
    % Adjust the alpha level using Bonferroni correction
    alpha_bonferroni = alpha / numComparisons;
    
    % Extract fixed effects
    fixed_Effects = lme.Coefficients;
    
    % Append model results to the summary for saving
    results_summary{end+1} = struct(...
        'Condition', condition, ...
        'Formula', lme.Formula, ...
        'Statistics', fixed_Effects, ...
        'BonferroniAlpha', alpha_bonferroni ...
        );
end

%% Save all results to an Excel file
% results_filename = fullfile(outputPath, 'Granger_LMEM_results_3-5Hz_VTC_IPL_4periods.xlsx');
results_filename = fullfile(outputPath, 'Granger_LMEM_results_9-11Hz_Hip_IPL_4periods.xlsx');
% results_filename = fullfile(outputPath, 'Granger_LMEM_results_9-12Hz_Hip_IPL.xlsx');
for i = 1:length(results_summary)
    cond = results_summary{i}.Condition;
    fixedEffects = results_summary{i}.Statistics;
    alpha_bonferroni = results_summary{i}.BonferroniAlpha;
    
    % Add Bonferroni correction and significance
    fixedEffects.Significance = fixedEffects.pValue < alpha_bonferroni;
    
    % Convert fixed effects to table format
    fixedEffectsTable = dataset2table(fixedEffects);  % Convert to table
    
    % Write results for this condition to Excel
    writetable(fixedEffectsTable, results_filename, 'Sheet', cond);
end

%% plot mean +-sem of granger net flow across periods and conditions
% Assuming T_conditions contains the tables for each condition ('same', 'diff')
conditions = {'same', 'diff'};
% periods = {'encoding', 'delay1', 'delay2', 'recall'};
periods = {'baseline', 'encoding', 'delay1', 'delay2'};
colors = [1 0.5 0.5; 0.5 0.5 1]; % Red for 'same', blue for 'diff'
shiftWidth = 0.1;

% Initialize figure
figure;
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 20, 15]); % [x_position, y_position, width, height]
hold on;

% Iterate over conditions
for c = 1:length(conditions)
    condition = conditions{c};
    T = T_conditions.(condition);
    
    % Initialize arrays for storing mean and SEM values across all periods
    meanValues = zeros(1, numel(periods));
    semValues = zeros(1, numel(periods));
    
    % Iterate over periods
    for p_idx = 1:length(periods)
        period = periods{p_idx};
        
        % Extract the Granger data for the current period
        periodData = T.MeanNetGranger(T.Period == period);
        
        % Compute mean and SEM for this period
        meanValues(p_idx) = mean(periodData);
        semValues(p_idx) = std(periodData) / sqrt(length(periodData)); % SEM = std / sqrt(N)
    end
    
    % Shift x positions for conditions to avoid overlap
    xPositions = (1:length(periods)) + (c - 1.5) * shiftWidth;
    
    % Plot the mean with error bars for each condition
    errorbar(xPositions, meanValues, semValues, 'o', 'Color', colors(c,:),...
        'DisplayName', conditions{c}, 'LineWidth', 1.5, 'MarkerSize', 7, 'CapSize', 7, 'MarkerFaceColor', colors(c,:));
end

% Customize the plot
ax = gca; % Get current axis
ax.XTick = 1:length(periods);
% ax.XTickLabel = {'Encoding', 'Delay 1', 'Delay 2', 'Recall'};
ax.XTickLabel = {'Baseline', 'Encoding', 'Delay 1', 'Delay 2'};
ax.FontSize = 14; % Increase font size for better readability
ax.LineWidth = 1.5; % Make the axes lines thicker
ax.XAxis.FontSize = 14; % Specific font size for X-axis
ax.YAxis.FontSize = 14; % Specific font size for Y-axis

% Rotate X-tick labels
xtickangle(20);

% Labeling
ylabel('Mean Net Granger Flow in 9-11 Hz', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Mean Net Granger Flow in 7-8 Hz', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time Period', 'FontSize', 14, 'FontWeight', 'bold');
% title('Mean ± SEM of Granger Net Flow for VTC -> IPL', 'FontSize', 14, 'FontWeight', 'bold');
title('Mean ± SEM of Granger Net Flow for Hip -> IPL', 'FontSize', 14, 'FontWeight', 'bold');

% Set axis limits
xlim([0.5, length(periods) + 0.5]);
ylim([-1.4*max(abs(meanValues)) 1.4*max(abs(meanValues))]); % Adjust limits based on your data
% yticks([-0.02:0.01:0.02])

% Add legend
legend(conditions, 'Location', 'Best', 'Interpreter', 'None', 'FontSize', 14);
hold on;
% Add dashed line at y=0
plot([0, length(periods) + 0.5], [0, 0], '--', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off');

% Save figure
outputPath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\Granger_Sept2024';
% figureName = 'meanGrangerNetFlow_3-5Hz_VTC_IPL_over_all_chnPairs_4periods.jpg';
figureName = 'meanGrangerNetFlow_9-11Hz_Hip_IPL_over_all_chnPairs_4periods.jpg';
% figureName = 'meanGrangerNetFlow_9-12Hz_Hip_IPL_over_all_chnPairs.jpg';
print(gcf, fullfile(outputPath, figureName), '-djpeg', '-r300');

% Save figure in .fig format
% savefig(fullfile(outputPath, 'meanGrangerNetFlow_3-5Hz_VTC_IPL_over_all_chnPairs_4periods.fig'));
savefig(fullfile(outputPath, 'meanGrangerNetFlow_9-11Hz_Hip_IPL_over_all_chnPairs_4periods.fig'));
% savefig(fullfile(outputPath, 'meanGrangerNetFlow_9-12Hz_Hip_IPL_over_all_chnPairs.fig'));

%% fit more complex on lmem with period and condition as fixed effects

% combine data from 2 tables
T_conditions.same.Condition = categorical(repmat({'same'}, height(T_conditions.same), 1));
T_conditions.diff.Condition = categorical(repmat({'diff'}, height(T_conditions.diff), 1));
T = [T_conditions.same; T_conditions.diff];

% Ensure categorical variables are correctly specified
T.Period = categorical(T.Period, {'baseline','encoding', 'delay1', 'delay2'}, 'Ordinal', false);
T.Condition = categorical(T.Condition, {'same', 'diff'}, 'Ordinal', false);

% fit the model
lme = fitlme(T, 'MeanNetGranger ~ Period*Condition + (1|Participant) + (1|Participant:ChannelPair)');
disp(lme);

% save to xls
fixedEff2 = lme.Coefficients;
% Convert fixed effects to table format
fixedEff2Table = dataset2table(fixedEff2);  
% Write results  to Excel
writetable(fixedEff2Table, results_filename, 'Sheet', 'combined');

%% to determine if the "diff" condition during encoding are significantly different from the baseline period under the "diff" condition
% Test the hypothesis that the sum of Period_encoding + Period_encoding:Condition_diff = 0
H = [0 1 0 0 0 1 0 0]; % coefficients to test, corresponding to Period_encoding and Period_encoding:Condition_diff
[pValue_encoding,F]  = coefTest(lme, H);
disp(['p-value for "diff" encoding vs. baseline (diff): ', num2str(pValue_encoding)]);
disp(['F-stat for "diff" encoding vs. baseline (diff): ', num2str(F)]);

%% Test for delay 1 in the "diff" condition vs. baseline in the "diff" condition

% Coefficients corresponding to:
% [Intercept, Period_encoding, Period_delay1, Period_delay2, Condition_diff, Period_encoding:Condition_diff, Period_delay1:Condition_diff, Period_delay2:Condition_diff]
% We are interested in Period_delay1 and Period_delay1:Condition_diff (positions 3 and 7 in the coefficient vector)

H_delay1 = [0 0 1 0 0 0 1 0];  % Test Period_delay1 + Period_delay1:Condition_diff = 0
[p_value_delay1,F] = coefTest(lme, H_delay1);
disp(['p-value for "diff" delay 1 vs. baseline (diff): ', num2str(p_value_delay1)]);
disp(['F-stat for "diff" delay 1 vs. baseline (diff): ', num2str(F)]);

%% Test for delay 2 in the "diff" condition vs. baseline in the "diff" condition

% We are interested in Period_delay2 and Period_delay2:Condition_diff (positions 4 and 8 in the coefficient vector)

H_delay2 = [0 0 0 1 0 0 0 1];  % Test Period_delay2 + Period_delay2:Condition_diff = 0
[p_value_delay2,F] = coefTest(lme, H_delay2);
disp(['p-value for "diff" delay 2 vs. baseline (diff): ', num2str(p_value_delay2)]);
disp(['F-stat for "diff" delay 2 vs. baseline (diff): ', num2str(F)]);

%% extract mean granger values for each direction

% Initialize a structure to store data tables for each condition
T_directions_conditions = struct();

% Loop over conditions
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Initialize variables for saving results for this condition
    MeanGranger = [];
    Direction1 = [];
    Direction2 = [];
    Participant = {};
    ChannelPair = {};
    Period = {};
    
    % Loop over periods
    for p_idx = 1:length(periods)
        period = periods{p_idx};
        
        % Loop over subjects
        for p = 1:numel(pacienti)
            if ~pacienti(p).todo
                continue;
            end
            
            % Get patient folder name
            patientFolder = pacienti(p).folder;
            % Construct the file path
            filePath = [basedir patientFolder '\' subfolder '\' Gr_folder '\'];
            
            % Determine file name based on condition and period
            if strcmp(condition, 'same')
                filename = Granger_files_s{p_idx};
            else
                filename = Granger_files_d{p_idx};
            end
            fullFileName = fullfile(filePath, filename);
            
            if isfile(fullFileName)
                % Load data
                data = load(fullFileName);
                grangerPeriod = data.grangerPeriod;
                grangerspctrm = grangerPeriod.grangerspctrm;
                labelcmb = grangerPeriod.labelcmb;
                freq = grangerPeriod.freq;
                freq_idx = find(freq >= freq_band(1) & freq <= freq_band(2));
                N_rows = size(grangerspctrm, 1);
                N_pairs = N_rows / 2; % Number of channel pairs
                N_freqs = length(freq);
                
                % Reshape to [N_pairs x 2 x N_freqs]
                grangerspctrm_reshaped = reshape(grangerspctrm, [2, N_pairs, N_freqs]);
                grangerspctrm_reshaped = permute(grangerspctrm_reshaped, [2, 1, 3]); % [N_pairs x 2 x N_freqs]
                mean_granger_freqband = mean(grangerspctrm_reshaped(:,:,freq_idx), 3); % [N_pairs x 2]
                
                for chnp = 1:N_pairs
                    chlabel1 = split(labelcmb{2*(chnp-1)+1,1}, '[');
                    chlabel2 = split(labelcmb{2*(chnp-1)+2,1}, '[');
                    pair_label = sprintf('%s-%s', chlabel1{1}, chlabel2{1});
                    
                    ichn_pair = data.ROI_chanpairs(chnp, :); % number of channel in chn pair (N pairs x 2)
                    
                    % Determine the direction order based on ROI labels
                    if strcmp(data.dataPeriod.channelInfo(ichn_pair(2)).ROI, ROI1) % first direction in granger spectrum
                        direction1 = [ROI1 '->' ROI2];
                        direction2 = [ROI2 '->' ROI1];
                    elseif strcmp(data.dataPeriod.channelInfo(ichn_pair(2)).ROI, ROI2)
                        direction1 = [ROI2 '->' ROI1];
                        direction2 = [ROI1 '->' ROI2];
                    end
                    
                    % values for each direction
                    dir1 = mean_granger_freqband(chnp,1);
                    dir2 = mean_granger_freqband(chnp,2);
                    
                    % Store data
                    MeanGranger = [MeanGranger; dir1, dir2];
                    Participant{end+1} = patientFolder;
                    ChannelPair{end+1} = pair_label;
                    Period{end+1} = period;
                    Direction1{end+1} = direction1;
                    Direction2{end+1} = direction2;
                end
            else
                warning('File %s not found for patient %s', filename, patientFolder);
            end
        end
    end
    
    % Create table for this condition and store it in the structure
    T_directions = table(MeanGranger(:,1), MeanGranger(:,2), Participant', ChannelPair', Period', Direction1', Direction2', ...
        'VariableNames', {'GrangerDir1', 'GrangerDir2', 'Participant', 'ChannelPair', 'Period', 'Direction1', 'Direction2'});
    T_directions.Participant = categorical(T_directions.Participant);
    T_directions.ChannelPair = categorical(T_directions.ChannelPair);
    T_directions.Period = categorical(T_directions.Period, {'baseline', 'encoding', 'delay1', 'delay2'}, 'Ordinal', false);
    % Store the table in the structure
    T_directions_conditions.(condition) = T_directions;
end
%% plot each direction
% Define colors for the plot
colors_same = [1 0.5 0.5; 1 0 0]; % Light red for Dir1, Dark red for Dir2
colors_diff = [0.5 0.5 1; 0 0 1]; % Light blue for Dir1, Dark blue for Dir2
periods = {'baseline', 'encoding', 'delay1', 'delay2'};
shiftWidth = 0.1;

% Initialize figure
figure;
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 35, 10]); % [x_position, y_position, width, height]

% Loop through each condition and create subplots
for c = 1:length(conditions)
    condition = conditions{c};
    T_directions = T_directions_conditions.(condition);
    
    % Determine the appropriate colors for the condition
    if strcmp(condition, 'same')
        colors = colors_same;
    else
        colors = colors_diff;
    end
    
    % Initialize mean and SEM arrays for both directions
    meanValuesDir1 = zeros(1, numel(periods));
    meanValuesDir2 = zeros(1, numel(periods));
    semValuesDir1 = zeros(1, numel(periods));
    semValuesDir2 = zeros(1, numel(periods));
    
    % Loop over each period
    for p_idx = 1:length(periods)
        period = periods{p_idx};
        
        % Extract Granger values for the current period
        periodDataDir1 = T_directions.GrangerDir1(T_directions.Period == period);
        periodDataDir2 = T_directions.GrangerDir2(T_directions.Period == period);
        
        % Compute mean and SEM for each direction
        meanValuesDir1(p_idx) = mean(periodDataDir1);
        meanValuesDir2(p_idx) = mean(periodDataDir2);
        semValuesDir1(p_idx) = std(periodDataDir1) / sqrt(length(periodDataDir1)); % SEM = std / sqrt(N)
        semValuesDir2(p_idx) = std(periodDataDir2) / sqrt(length(periodDataDir2)); % SEM = std / sqrt(N)
    end
    
    % Create subplot for each condition
    subplot(1, 2, c);
    hold on;
    
    % Shift x positions to avoid overlap of two directions
    xPositions = 1:length(periods);
    
    % Plot for Direction 1 (first direction)
    h1=errorbar(xPositions, meanValuesDir1, semValuesDir1, 'o', 'Color', colors(1,:),...
        'LineWidth', 1.5, 'MarkerSize', 7, 'CapSize', 7, 'MarkerFaceColor', colors(1,:),...
        'DisplayName', T_directions.Direction1{1});
    
    % Plot for Direction 2 (second direction)
    h2=errorbar(xPositions + shiftWidth, meanValuesDir2, semValuesDir2, 'o', 'Color', colors(2,:),...
        'LineWidth', 1.5, 'MarkerSize', 7, 'CapSize', 7, 'MarkerFaceColor', colors(2,:),...
        'DisplayName', T_directions.Direction2{1});
    
    % Customize the plot
    ax = gca; % Get current axis
    ax.XTick = 1:length(periods);
    ax.XTickLabel = {'Baseline', 'Encoding', 'Delay 1', 'Delay 2'};
    ax.FontSize = 14; % Increase font size for better readability
    ax.LineWidth = 1.5; % Make the axes lines thicker
    ax.XAxis.FontSize = 14; % Specific font size for X-axis
    ax.YAxis.FontSize = 14; % Specific font size for Y-axis
    
    % Rotate X-tick labels
    xtickangle(20);
    
    % Labeling
    ylabel('Mean Granger Value (9-11 Hz)', 'FontSize', 14, 'FontWeight', 'bold');
   
    
    % Title for each subplot
    if strcmp(condition, 'same')
        title('Same Condition', 'FontSize', 14, 'FontWeight', 'bold');
    else
        title('Different Condition', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Set axis limits
    xlim([0.5, length(periods) + 0.5]);
%     ylim([0.002 0.006])
    % Add legend
    legend([h1, h2], {T_directions.Direction1{1}, T_directions.Direction2{1}}, 'Location', 'Best', 'Interpreter', 'None', 'FontSize', 14);
    

end

% Save figure
outputPath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\Granger_Sept2024';
% figureName = 'Granger_Mean_Directions_Conditions_3-5Hz_VTC_IPL_4periods.jpg';
figureName = 'Granger_Mean_Directions_Conditions_9-11Hz_Hip_IPL_4periods.jpg';
print(gcf, fullfile(outputPath, figureName), '-djpeg', '-r300');

% Save figure in .fig format
% savefig(fullfile(outputPath, 'Granger_Mean_Directions_Conditions_3-5Hz_VTC_IPL_4periods.fig'));
savefig(fullfile(outputPath, 'Granger_Mean_Directions_Conditions_9-11Hz_Hip_IPL_4periods.fig'));

