% --------------------
% a script to do Linear mixed effect model (LMEM) for mean PLV (all trials together) for all time periods - bs, first and second half of the delay, recall - fixed effects,
% while accounting for random effects - ch pair and patientID
% by Sofiia Moraresku
% May 2024
% --------------------

% first, set up the indexes of freq bins to use for averaging PLV in the model (obtained earlier by plot_PLVdiff_variability)
% isignificant_bins = [1 2 3]; % delta-theta, 2-4Hz, Hip-IPL
% isignificant_bins = [6 7]; % alpha, 7-8Hz, Hip-IPL
isignificant_bins = [1:4]; % 2-5 Hz, VTC-IPL

% load an aggregated table with indices of all significant pairs for all three periods across patients
filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_VTC-IPL_aggregated_significant_chnPairs.mat';
% filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_Hip-IPL_aggregated_significant_chnPairs.mat';
load(filenameTable);

% path for the output files
outputPath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\LMEM_May2024';

% set up a path for the main data
PLV_data1 =  'PLV_VTC-IPL_first 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
PLV_data2 =  'PLV_VTC-IPL_last 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
PLV_data3 =  'PLV_VTC-IPL_0.5s recall_vs_bs_all_trials_200permut_2024-04.mat';
% PLV_data1 =  'PLV_Hip-IPL_first 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
% PLV_data2 =  'PLV_Hip-IPL_last 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
% PLV_data3 =  'PLV_Hip-IPL_0.5s recall_vs_bs_all_trials_200permut_2024-04.mat';
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % set up pacienti_memact by selecting which patients to analyze

% Extract important info from the file name
ROI1 = regexp(PLV_data1, '(?<=PLV_)[^-]+', 'match', 'once');
ROI2 = regexp(PLV_data1, '(?<=-)[^-_]+', 'match', 'once');

%% get average PLV over the selected freq range for all periods
% Initialize an empty table for storing the results
resultsTable = table([], [], [], [], 'VariableNames', {'patientId', 'chPair', 'timePeriod', 'meanPLV'});

% Iterate over each patient
for p = 1:numel(pacienti)
    
    if ~pacienti(p).todo
        continue; % Skip if the subject is not marked 'todo'
    end
    
    % Match the patient by name
    subjIndex = find(strcmp({aggregTable.patient}, pacienti(p).folder));
    if isempty(subjIndex)
        continue;
    end
    
    filePath = [basedir pacienti(p).folder '/' subfolder '/' PLV_folder '/'];
    
    if isfile([filePath PLV_data1])
        % Load the PLV data for all 3 periods
        first_delay_plv = load([filePath PLV_data1]);
        last_delay_plv = load([filePath PLV_data2]);
        recall_plv = load([filePath PLV_data3]);
        
        % Get indices of aggragated significant channel pairs for the current subject
        significantidxChan = aggregTable(subjIndex).idxChan_in_Pairs; % n signif pairs x 2 (each pair with chan indices)
        
        % For each significant pair extract raw PLV for each period
        for i = 1:size(significantidxChan,1)
            ipair = significantidxChan(i,:);
            
            % Calculate the mean PLV within the selected frequency range for this pair for all periods
            meanPLV_first_delay = mean(squeeze(first_delay_plv.PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            %             meanPLV_bs = mean(squeeze(first_delay_plv.PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_last_delay = mean(squeeze(last_delay_plv.PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_recall_plv = mean(squeeze(recall_plv.PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_bs = mean(squeeze(recall_plv.PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins))); % short 0.5s baseline
            
            chn_labels = [first_delay_plv.PLVCond1.label{ipair(1)} first_delay_plv.PLVCond1.label{ipair(2)}];  % original labels of chan in this pair
            
            % Append the results to the table
            newRow1 = {pacienti(p).folder, chn_labels, 'bs', meanPLV_bs}; % baseline
            newRow2 = {pacienti(p).folder, chn_labels, 'first_1.9s_delay', meanPLV_first_delay}; % first half of the delay
            newRow3 = {pacienti(p).folder, chn_labels, 'last_1.9s_delay', meanPLV_last_delay}; % last half of the delay
            newRow4 = {pacienti(p).folder, chn_labels, 'recall', meanPLV_recall_plv};
            resultsTable = [resultsTable; newRow1; newRow2; newRow3; newRow4];
        end
    end
end

resultsTable.timePeriod = categorical(resultsTable.timePeriod);
resultsTable.patientId = categorical(resultsTable.patientId);
resultsTable.chPair = categorical(resultsTable.chPair);

%% LMM
% Fit a linear mixed-effects model: time period(bs; delay) - fixed effect, patient and chan pair - random effects
model = fitlme(resultsTable, 'meanPLV ~ timePeriod + (1 | patientId) + (1 | patientId:chPair)'); % chapter 36.4.2 of iEEG book

% Display the model summary
disp(model); % the output can be copied to txt file

%% Test the significance of fixed effects coefficients
% The coefTest function in MATLAB conducts a linear hypothesis test of the model coefficients
% pVal = coefTest(model); % Test if all fixed-effects coefficients except for the intercept are 0.
pVal1 = coefTest(model,[0 1 -1 0]); % Test if there is any difference between the first and last delay
pVal2 = coefTest(model,[0 0 1 -1]); % Test if there is any difference between the last delay and recall
pVal3 = coefTest(model,[0 1 0 -1]); % Test if there is any difference between the first delay and recall

% Define the number of comparisons (post-hoc)
numComparisons = 3;

% original alpha level
alpha = 0.05;

% Adjust the alpha level using Bonferroni correction
alphaAdjusted_post_hoc = alpha / numComparisons;
significant_post_hoc = [pVal1; pVal2; pVal3] < alphaAdjusted_post_hoc;  % Apply Bonferroni correction

% % FDR correction
% [~, ~, adj_p_values] = fdr_bh([pVal1, pVal2, pVal3], 0.05, 'pdep');

%% export the results to the xls file

% Define the Excel file name
outputname = ['LMEM_output_meanPLV_' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '.xlsx'];
filename = [outputPath '\' outputname];

% Define the number of comparisons (models in total)
numComparisons = 3;
% original alpha level
alpha = 0.05;
% Adjust the alpha level using Bonferroni correction
alphaAdjusted = alpha / numComparisons;

% Manually create a table from the 'titleddataset'
coefficients = model.Coefficients;
significant = coefficients.pValue < alphaAdjusted;  % Apply Bonferroni correction
coefTable = table(coefficients.Name, coefficients.Estimate, coefficients.SE, coefficients.tStat, ...
    coefficients.DF,  coefficients.Lower, coefficients.Upper, coefficients.pValue, significant, repmat(alphaAdjusted, numel(significant),1),...
    'VariableNames', {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'Lower', 'Upper', 'pValue', 'significance_bonferroni_corr', 'alpha_bonferroni'});

% Write the fixed effects table to Excel
writetable(coefTable, filename, 'Sheet', 'Fixed Effects', 'Range', 'A1');

% Model Statistics
modelStats = {
    'Number of observations', size(model.Residuals.Raw, 1); % Adjusted to model.Residuals.Raw if applicable
    'AIC', model.ModelCriterion.AIC;
    'BIC', model.ModelCriterion.BIC;
    'LogLikelihood', model.LogLikelihood;
    };
modelStatsTable = cell2table(modelStats, 'VariableNames', {'Statistic', 'Value'});
writetable(modelStatsTable, filename, 'Sheet', 'Model Statistics', 'Range', 'A1');

% Model Formula
formula = {char(model.Formula)};
formulaTable = cell2table(formula, 'VariableNames', {'Formula'});
writetable(formulaTable, filename, 'Sheet', 'Model Formula', 'Range', 'A1');

% Prepare and save comparison results
comparisons = {
    'First 1.9s Delay vs Last 1.9s Delay';
    'Last 1.9s Delay vs Recall';
    'First 1.9s Delay vs Recall'
    };
comparisonResults = table(comparisons, [pVal1; pVal2; pVal3], significant_post_hoc, repmat(alphaAdjusted_post_hoc, numel(significant_post_hoc),1), ...
    'VariableNames', {'Comparison', 'Uncorrected_pValue', 'significance_bonferroni_corr', 'alpha_bonferroni'});

% Append the comparison results to the existing Excel file
writetable(comparisonResults, filename, 'Sheet', 'Comparison Results', 'Range', 'A1');

%% visualize the results
% Unique identifiers for patients and time periods
uniquePatients = unique(string(resultsTable.patientId), 'stable');
uniqueTimePeriods = unique(string(resultsTable.timePeriod));

% Initialize arrays for storing mean and SEM values
meanValues = zeros(length(uniquePatients), length(uniqueTimePeriods));
semValues = zeros(size(meanValues));

% Compute mean and SEM for each patient and time period
for i = 1:length(uniquePatients)
    for j = 1:length(uniqueTimePeriods)
        % Logical indices for the current patient and time period
        patientIdx = strcmp(string(resultsTable.patientId), uniquePatients(i));
        timePeriodIdx = string(resultsTable.timePeriod) == uniqueTimePeriods(j);
        
        % Extract PLV values for the current patient and time period
        plvData = resultsTable.meanPLV(patientIdx & timePeriodIdx);
        
        % Compute mean and SEM
        meanValues(i, j) = mean(plvData);
        semValues(i, j) = std(plvData) / sqrt(length(plvData));  % SEM = std / sqrt(N)
    end
end

% plot the mean and SEM for each patient
figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]); % Maximize the figure window
hold on;

% Define colors of markers
colors = lines(length(uniquePatients));

% Plot each patient's data
for i = 1:length(uniquePatients)
    % Plot the mean with error bars for each patient
    errorbar(1:numel(uniqueTimePeriods), meanValues(i, :), semValues(i, :), 'o-', ...
        'Color', colors(i, :), 'LineWidth', 1.5, 'MarkerSize', 6, 'CapSize', 8);
end

% Customize the plot
ax = gca; % Get current axis
ax.XTick = 1:numel(uniqueTimePeriods);
ax.XTickLabel = {'Baseline', 'First 1.9s Delay', 'Last 1.9s Delay', 'Recall'};
ax.FontSize = 16; % Increase font size for better readability
ax.LineWidth = 1.5; % Make the axes lines thicker
ax.XAxis.FontSize = 18; % Specific font size for X-axis
ax.YAxis.FontSize = 18; % Specific font size for Y-axis

% Labeling
ylabel(['Mean PLV in Freq Range: ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))  ' Hz'], ...
    'FontSize', 18);
xlabel('Time Period', 'FontSize', 18);
title(['Mean +- SEM of PLV in ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))...
    'Hz across chn pairs for each subject in ' ROI1 '-' ROI2], 'FontSize', 18);

% Set axis limits with some padding
xlim([0.5, numel(uniqueTimePeriods) + 0.5]);
ylim([min(meanValues(:) - semValues(:)) * 0.9, max(meanValues(:) + semValues(:)) * 1.05]);

% Add grid
% grid on;
% set(ax, 'GridLineStyle', '--'); % Customize grid lines

% Anti-aliasing for smoother lines
set(gcf, 'GraphicsSmoothing', 'on');

% Add legend
legend(uniquePatients, 'Location', 'Best', 'Interpreter', 'None', 'FontSize', 12);
hold off;

% Define the image name
figureName = ['meanPLV_' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '.jpg'];
filename2 = [outputPath '\' figureName];

% Save the figure as JPEG in high resolution
% print(gcf, filename2, '-djpeg', '-r300');


%% 22.05.2024 simpler plot with mean +- SEM across channel pairs
% Extract unique time periods
uniqueTimePeriods = unique(string(resultsTable.timePeriod), 'stable');
numPatients = numel(unique(string(resultsTable.patientId), 'stable'));
numTimePeriods = numel(uniqueTimePeriods);

% set the same y-limits
globalMin = 0.06;
globalMax = 0.18;

% Initialize figure
figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]); % Maximize the figure window
hold on;

% Initialize arrays for storing mean and SEM values across all time periods
meanValues = zeros(1, numTimePeriods);
semValues = zeros(1, numTimePeriods);

% Iterate over time periods
for timePeriodIdx = 1:numTimePeriods
    % Logical indices for the current time period
    timePeriodIdx_logical = string(resultsTable.timePeriod) == uniqueTimePeriods(timePeriodIdx);
    
    % Extract PLV values for the current time period
    plvData = resultsTable.meanPLV(timePeriodIdx_logical);
    
    % Compute mean and SEM
    meanValues(timePeriodIdx) = mean(plvData);
    semValues(timePeriodIdx) = std(plvData) / sqrt(length(plvData));  % SEM = std / sqrt(N)
end


% Plot the mean with error bars
errorbar(1:numTimePeriods, meanValues, semValues, 'o', 'Color', 'm',...
    'LineWidth', 2, 'MarkerSize', 10, 'CapSize', 10, 'MarkerFaceColor', 'm');

% Customize the plot
ax = gca; % Get current axis
ax.XTick = 1:numTimePeriods;
ax.XTickLabel = {'Baseline', 'First 1.9s Delay', 'Last 1.9s Delay', 'Recall'};
ax.FontSize = 12; % Increase font size for better readability
ax.LineWidth = 1.5; % Make the axes lines thicker
ax.XAxis.FontSize = 14; % Specific font size for X-axis
ax.YAxis.FontSize = 14; % Specific font size for Y-axis

% Labeling
ylabel(['Mean PLV in Freq Range: ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))  ' Hz'], ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time Period', 'FontSize', 14, 'FontWeight', 'bold');
title(['Mean +- SEM of PLV in ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))...
    'Hz in ' ROI1 '-' ROI2 ' over all chn pairs '], 'FontSize', 16, 'FontWeight', 'bold');

% Set axis limits with some padding
xlim([0.5, numTimePeriods + 0.5]);

% Set consistent y-axis limits
ylim([globalMin globalMax]);
yticks(globalMin:0.02:globalMax)

% Show grid
grid on;
set(ax, 'GridLineStyle', '--'); % Customize grid lines

% Release hold
hold off;

% Define the image name
figureName = ['meanPLV_' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '_over_all_chnPairs.jpg'];
filename2 = [outputPath '\' figureName];

% Save the figure as JPEG in high resolution
print(gcf, filename2, '-djpeg', '-r300');