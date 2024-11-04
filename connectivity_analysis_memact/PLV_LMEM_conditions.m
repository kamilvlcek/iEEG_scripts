% --------------------
% a script to run a Linear mixed effects model (LMEM) on mean PLV for 2 conditions (same and diff) and 4 time periods - fixed effects, 
% while accounting for random effects - ch pair and patientID
% by Sofiia Moraresku
% May 2024
% --------------------

% first, set up the indexes of freq bins to use for averaging PLV in the model (obtained earlier by plot_PLVdiff_variability)
isignificant_bins = [1 2 3]; % low theta, 2-4Hz, Hip-IPL
% isignificant_bins = [6 7]; % high theta, 7-8Hz, Hip-IPL
% isignificant_bins = [1:7]; % the whole range, Hip-IPL
% isignificant_bins = [1:4]; % low theta, 2-5 Hz, VTC-IPL

% load an aggregated table with indices of all significant pairs for all three periods across patients
% filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_VTC-IPL_aggregated_significant_chnPairs.mat';
% filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_VTC-IPL_aggregated_significant_chnPairs_Sept24.mat';
filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_Hip-IPL_aggregated_significant_chnPairs.mat';
% filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_Hip-IPL_aggregated_significant_chnPairs_Sept24.mat';
load(filenameTable);

% path for the output files
outputPath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\LMEM_May2024\conditions';

% set up a path for the main data
PLV_data1 =  'PLV_Hip-IPL_diff_vs_same_first 1.9s delay__trials_2024-05.mat';
PLV_data2 =  'PLV_Hip-IPL_diff_vs_same_last 1.9s delay__trials_2024-05.mat';
PLV_data3 =  'PLV_Hip-IPL_diff_vs_same_0.5s recall__trials_2024-05.mat';
PLV_data4 = 'PLV_Hip-IPL_diff_vs_same_0.5s bs__trials_2024-05.mat';
PLV_data5 = 'PLV_Hip-IPL_diff_vs_same_1.9s encoding__trials_2024-09.mat';
% PLV_data1 =  'PLV_VTC-IPL_diff_vs_same_first 1.9s delay__trials_2024-05.mat';
% PLV_data2 =  'PLV_VTC-IPL_diff_vs_same_last 1.9s delay__trials_2024-05.mat';
% PLV_data3 =  'PLV_VTC-IPL_diff_vs_same_0.5s recall__trials_2024-05.mat';
% PLV_data4 = 'PLV_VTC-IPL_diff_vs_same_0.5s bs__trials_2024-05.mat';
% PLV_data5 = 'PLV_VTC-IPL_diff_vs_same_1.9s encoding__trials_2024-09.mat';
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % set up pacienti_memact by selecting which patients to analyze

% Extract important info from the file name
ROI1 = regexp(PLV_data1, '(?<=PLV_)[^-]+', 'match', 'once');
ROI2 = regexp(PLV_data1, '(?<=-)[^-_]+', 'match', 'once');


%% get average PLV over the selected freq range for two conditions
% Initialize an empty table for storing the results
resultsTable = table([], [], [], [], [], 'VariableNames', {'patientId', 'chPair', 'timePeriod', 'condition', 'meanPLV'});

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
        % Load the file containing PLV data for two conditions for all periods
%         load([filePath PLV_data]);
        
        first_delay_plv = load([filePath PLV_data1]);
        last_delay_plv = load([filePath PLV_data2]);
        recall_plv = load([filePath PLV_data3]);
        bs_plv = load([filePath PLV_data4]);
        encoding_plv = load([filePath PLV_data5]);

        % Get indices of aggregated significant channel pairs for the current subject
        significantidxChan = aggregTable(subjIndex).idxChan_in_Pairs; % n signif pairs x 2 (each pair with chan indices)

        % For each significant pair extract raw PLV for each condition
        for i = 1:size(significantidxChan,1)
            ipair = significantidxChan(i,:);

            % Calculate the mean PLV within the selected frequency range for this pair for all conditions and periods
%             meanPLV_same = mean(squeeze(PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
%             meanPLV_diff = mean(squeeze(PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins))); 
            
            meanPLV_first_delay_same = mean(squeeze(first_delay_plv.PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_first_delay_diff = mean(squeeze(first_delay_plv.PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_last_delay_same = mean(squeeze(last_delay_plv.PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_last_delay_diff = mean(squeeze(last_delay_plv.PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_recall_same = mean(squeeze(recall_plv.PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_recall_diff = mean(squeeze(recall_plv.PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_bs_same = mean(squeeze(bs_plv.PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_bs_diff = mean(squeeze(bs_plv.PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_encod_same = mean(squeeze(encoding_plv.PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_encod_diff = mean(squeeze(encoding_plv.PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            
            chn_labels = [first_delay_plv.PLVCond1.label{ipair(1)} first_delay_plv.PLVCond1.label{ipair(2)}];  % original labels of chan in this pair

            % Append the results to the table
%             newRow1 = {pacienti(p).folder, chn_labels, 'same', meanPLV_same}; 
%             newRow2 = {pacienti(p).folder, chn_labels, 'diff', meanPLV_diff}; 
            newRow1 = {pacienti(p).folder, chn_labels, 'bs', 'same', meanPLV_bs_same}; % baseline
            newRow2 = {pacienti(p).folder, chn_labels, 'bs', 'diff', meanPLV_bs_diff}; % baseline
            newRow3 = {pacienti(p).folder, chn_labels, 'encoding', 'same', meanPLV_encod_same}; % encoding
            newRow4 = {pacienti(p).folder, chn_labels, 'encoding', 'diff', meanPLV_encod_diff}; % encoding            
            newRow5 = {pacienti(p).folder, chn_labels, 'first_1.9s_delay', 'same', meanPLV_first_delay_same}; % first half of the delay
            newRow6 = {pacienti(p).folder, chn_labels, 'first_1.9s_delay', 'diff', meanPLV_first_delay_diff}; % first half of the delay
            newRow7 = {pacienti(p).folder, chn_labels, 'last_1.9s_delay', 'same', meanPLV_last_delay_same}; % last half of the delay
            newRow8 = {pacienti(p).folder, chn_labels, 'last_1.9s_delay', 'diff', meanPLV_last_delay_diff}; % last half of the delay
            newRow9 = {pacienti(p).folder, chn_labels, 'recall', 'same', meanPLV_recall_same}; % recall
            newRow10 = {pacienti(p).folder, chn_labels, 'recall', 'diff', meanPLV_recall_diff}; % recall            
            resultsTable = [resultsTable; newRow1; newRow2; newRow3; newRow4; newRow5; newRow6; newRow7; newRow8; newRow9; newRow10];
        end
    end
end

resultsTable.timePeriod = categorical(resultsTable.timePeriod);
resultsTable.condition = categorical(resultsTable.condition);
resultsTable.patientId = categorical(resultsTable.patientId);
resultsTable.chPair = categorical(resultsTable.chPair);

%% LMM
% Fit a linear mixed-effects model: condition - fixed effect, patient and chan pair - random effects
model = fitlme(resultsTable, 'meanPLV ~ timePeriod * condition + (1 | patientId) + (1 | patientId:chPair)'); % chapter 36.4.2 of iEEG book

% Display the model summary
disp(model); % the output can be copied to txt file

% Define the number of comparisons (number of models)
numComparisons = 3;

% original alpha level
alpha = 0.05;

% Adjust the alpha level using Bonferroni correction
alphaAdjusted = alpha / numComparisons;

%% post-hoc encoding vs other periodss
% Test if there is a difference between 'encoding' and 'first_1.9s_delay'
pVal1 = coefTest(model, [0 1 -1 0 0 0 0 0 0 0]); 

% Test if there is a difference between 'encoding' and 'last_1.9s_delay'
pVal2 = coefTest(model, [0 1 0 -1 0 0 0 0 0 0]);

% Test if there is a difference between 'encoding' and 'recall'
pVal3 = coefTest(model, [0 1 0 0 -1 0 0 0 0 0]);

% Test if there is a difference between 'first_1.9s_delay' and 'last_1.9s_delay'
pVal4 = coefTest(model, [0 0 1 -1 0 0 0 0 0 0]);

% Number of comparisons (post-hoc testing)
numComparisons = 4;

% Original alpha level
alpha = 0.05;

% Adjust alpha using Bonferroni correction
alphaAdjusted_post_hoc = alpha / numComparisons;

% Check significance with Bonferroni correction
significant_post_hoc = [pVal1; pVal2; pVal3; pVal4] < alphaAdjusted_post_hoc;

% Display the results
disp('p-values for pairwise comparisons:');
disp([pVal1, pVal2, pVal3, pVal4]);

disp('Significant after Bonferroni correction:');
disp(significant_post_hoc);


%% export the results to the xls file

% Define the Excel file name
outputname = ['LMEM_output_meanPLV_' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '_allPeriods.xlsx'];
filename = [outputPath '\' outputname];

% Manually create a table from the 'titleddataset'
coefficients = model.Coefficients;
significant = coefficients.pValue < alphaAdjusted;  % Apply Bonferroni correction
coefTable = table(coefficients.Name, coefficients.Estimate, coefficients.SE, coefficients.tStat, ...
                  coefficients.DF, coefficients.Lower, coefficients.Upper, coefficients.pValue, significant, repmat(alphaAdjusted, numel(significant),1), ...
                  'VariableNames', {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'Lower', 'Upper', 'pValue','significance_bonferroni_corr', 'alpha_bonferroni'});

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

%% visualize the results (mean +- sem for each patient separately)
% Unique identifiers for patients and conditions
uniquePatients = unique(string(resultsTable.patientId), 'stable');
uniqueCondition = unique(string(resultsTable.condition), 'stable');

% Initialize arrays for storing mean and SEM values
meanValues = zeros(length(uniquePatients), length(uniqueCondition));
semValues = zeros(size(meanValues));

% Compute mean and SEM for each patient and condition
for i = 1:length(uniquePatients)
    for j = 1:length(uniqueCondition)
        % Logical indices for the current patient and condition
        patientIdx = strcmp(string(resultsTable.patientId), uniquePatients(i));
        conditionIdx = string(resultsTable.condition) == uniqueCondition(j);
        
        % Extract PLV values for the current patient and condition
        plvData = resultsTable.meanPLV(patientIdx & conditionIdx);
        
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
    errorbar(1:numel(uniqueCondition), meanValues(i, :), semValues(i, :), 'o-', ...
        'Color', colors(i, :), 'LineWidth', 1.5, 'MarkerSize', 6, 'CapSize', 8); 
end

% Customize the plot
ax = gca; % Get current axis
ax.XTick = 1:numel(uniqueCondition);
ax.XTickLabel = {'Same', 'Diff'};
ax.FontSize = 12; % Increase font size for better readability
ax.LineWidth = 1.5; % Make the axes lines thicker
ax.XAxis.FontSize = 14; % Specific font size for X-axis
ax.YAxis.FontSize = 14; % Specific font size for Y-axis

% Labeling
ylabel(['Mean PLV in Freq Range: ' num2str(PLVCond2.freq(isignificant_bins(1))) '-' num2str(PLVCond2.freq(isignificant_bins(end)))  ' Hz'], ...
    'FontSize', 14); 
xlabel('Condition', 'FontSize', 14); 
title(['Mean +- SEM of PLV in ' num2str(PLVCond2.freq(isignificant_bins(1))) '-' num2str(PLVCond2.freq(isignificant_bins(end)))...
    'Hz across chn pairs for each subject in ' ROI1 '-' ROI2 ' during ' period], 'FontSize', 14); 

% Set axis limits with some padding
xlim([0.5, numel(uniqueCondition) + 0.5]);
ylim([min(meanValues(:) - semValues(:)) * 0.9, max(meanValues(:) + semValues(:)) * 1.05]);

% Add grid
grid on;
set(ax, 'GridLineStyle', '--'); % Customize grid lines

% Anti-aliasing for smoother lines
set(gcf, 'GraphicsSmoothing', 'on'); 

% Add legend
legend(uniquePatients, 'Location', 'Best', 'Interpreter', 'None', 'FontSize', 10);
hold off;

% Define the image name
figureName = ['meanPLV_' num2str(PLVCond2.freq(isignificant_bins(1))) '-' num2str(PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '_' period '.jpg'];
filename2 = [outputPath '\' figureName];

% Save the figure as JPEG in high resolution
print(gcf, filename2, '-djpeg', '-r300'); 


%% visualize results - all conditions and time periods for each patient individually
% Extract unique time periods and conditions
uniqueTimePeriods = unique(string(resultsTable.timePeriod), 'stable');
uniqueConditions = unique(string(resultsTable.condition), 'stable');
uniquePatients = unique(string(resultsTable.patientId), 'stable');
numPatients = numel(uniquePatients);

% Initialize figure and subplots for each condition
figure;
for conditionIdx = 1:numel(uniqueConditions)
%     subplot(numel(uniqueConditions), 1, conditionIdx);
    figure;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]); % Maximize the figure window
    hold on;
    maxPLV = zeros(1, numPatients);
    % Plot each patient's data
    for patientIdx = 1:numPatients
        % Initialize arrays for storing mean and SEM values
        meanValues = zeros(1, numel(uniqueTimePeriods));
        semValues = zeros(1, numel(uniqueTimePeriods));

        % Compute mean and SEM for each time period
        for timePeriodIdx = 1:numel(uniqueTimePeriods)
            % Logical indices for the current patient, time period, and condition
            patientIdx_logical = string(resultsTable.patientId) == uniquePatients(patientIdx);
            timePeriodIdx_logical = string(resultsTable.timePeriod) == uniqueTimePeriods(timePeriodIdx);
            conditionIdx_logical = string(resultsTable.condition) == uniqueConditions(conditionIdx);
            
            % Extract PLV values for the current patient, time period, and condition
            plvData = resultsTable.meanPLV(patientIdx_logical & timePeriodIdx_logical & conditionIdx_logical);
            
            % Compute mean and SEM
            meanValues(timePeriodIdx) = mean(plvData);
            semValues(timePeriodIdx) = std(plvData) / sqrt(length(plvData));  % SEM = std / sqrt(N)
        end

        % Plot the mean with error bars for each patient
        errorbar(1:numel(uniqueTimePeriods), meanValues, semValues, 'o-', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'CapSize', 8); 
        maxPLV(patientIdx) = max(meanValues);
    end

    % Customize the plot
    ax = gca; % Get current axis
    ax.XTick = 1:numel(uniqueTimePeriods);
%     ax.XTickLabel = uniqueTimePeriods;
    ax.XTickLabel = {'Baseline', 'First 1.9s Delay', 'Last 1.9s Delay', 'Recall'};
    ax.FontSize = 12; % Increase font size for better readability
    ax.LineWidth = 1.5; % Make the axes lines thicker
    ax.XAxis.FontSize = 14; % Specific font size for X-axis
    ax.YAxis.FontSize = 14; % Specific font size for Y-axis
    
    % Labeling
    ylabel(['Mean PLV in Freq Range: ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))  ' Hz'], ...
        'FontSize', 14);
    xlabel('Time Period', 'FontSize', 14);
    title(['Mean +- SEM of PLV in ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))...
        'Hz in ' ROI1 '-' ROI2 ' for Condition: ' uniqueConditions{conditionIdx}], 'FontSize', 14);

    % Set axis limits with some padding
    xlim([0.5, numel(uniqueTimePeriods) + 0.5]);
%     ylim([min(meanValues(:) - semValues(:)) * 0.9, max(maxPLV) * 1.15]);
    ylim([0.09 0.23]);

    % Add grid
    grid on;
    set(ax, 'GridLineStyle', '--'); % Customize grid lines
    
    % Add legend
    legend(uniquePatients, 'Location', 'Best', 'Interpreter', 'None', 'FontSize', 10);
    hold off;
    
    % Define the image name
    figureName = ['meanPLV_' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '_condition_' uniqueConditions{conditionIdx} '.jpg'];
    filename2 = [outputPath '\' figureName];
    
    % Save the figure as JPEG in high resolution
    print(gcf, filename2, '-djpeg', '-r300');
    
end

%% 22.05.2024 simpler plot with mean +- SEM across channel pairs
% Extract unique time periods and conditions
uniqueTimePeriods = unique(string(resultsTable.timePeriod), 'stable');
uniqueConditions = unique(string(resultsTable.condition), 'stable');
numPatients = numel(unique(string(resultsTable.patientId), 'stable'));
numTimePeriods = numel(uniqueTimePeriods);
numConditions = numel(uniqueConditions);

% Define colors for conditions
colors = [
    1 0.5 0.5; % same condition (light red)
    0.5 0.5 1; % different condition (light blue)
    ];

% Width of the shift for conditions
shiftWidth = 0.1;

% set the same y-limits
globalMin = 0.1; 
globalMax = 0.21;

% Initialize figure
figure;
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 17, 13]); % [x_position, y_position, width, height]
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]); % Maximize the figure window
hold on;

% Initialize arrays to store maximum PLV values for each condition
maxPLV = zeros(1, numConditions);

% Iterate over conditions
for conditionIdx = 1:numConditions
    % Initialize arrays for storing mean and SEM values across all time periods
    meanValues = zeros(1, numTimePeriods);
    semValues = zeros(1, numTimePeriods);
    
    % Iterate over time periods
    for timePeriodIdx = 1:numTimePeriods
        % Logical indices for the current time period and condition
        timePeriodIdx_logical = string(resultsTable.timePeriod) == uniqueTimePeriods(timePeriodIdx);
        conditionIdx_logical = string(resultsTable.condition) == uniqueConditions(conditionIdx);

        % Extract PLV values for the current time period and condition
        plvData = resultsTable.meanPLV(timePeriodIdx_logical & conditionIdx_logical);
        
        % Compute mean and SEM
        meanValues(timePeriodIdx) = mean(plvData);
        semValues(timePeriodIdx) = std(plvData) / sqrt(length(plvData));  % SEM = std / sqrt(N)
    end

    % Shift the x positions for conditions to avoid overlap
    xPositions = (1:numTimePeriods) + (conditionIdx - 1.5) * shiftWidth;
    
    % Plot the mean with error bars for each condition
    errorbar(xPositions, meanValues, semValues, 'o', 'Color', colors(conditionIdx,:),...
        'DisplayName', char(uniqueConditions(conditionIdx)), 'LineWidth', 1.5, 'MarkerSize', 7, 'CapSize', 7, 'MarkerFaceColor', colors(conditionIdx,:));
    
    % Update max PLV for this condition
    maxPLV(conditionIdx) = max(meanValues);
end

% Customize the plot
ax = gca; % Get current axis
ax.XTick = 1:numTimePeriods;
ax.XTickLabel = {'Baseline', 'Encoding', 'Delay 1', 'Delay 2', 'Recall'};
ax.FontSize = 14; % Increase font size for better readability
ax.LineWidth = 1.5; % Make the axes lines thicker
ax.XAxis.FontSize = 16; % Specific font size for X-axis
ax.YAxis.FontSize = 16; % Specific font size for Y-axis

% Rotate X-tick labels to 45 degrees
xtickangle(20);

% Labeling
ylabel(['Mean PLV in Freq Range: ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))  ' Hz'], ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time Period', 'FontSize', 14, 'FontWeight', 'bold');
title(['Mean +- SEM of PLV in ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))...
    'Hz in ' ROI1 '-' ROI2 ' over all chn pairs '], 'FontSize', 16, 'FontWeight', 'bold');

% Set axis limits with some padding
xlim([0.5, numTimePeriods + 0.5]);
% ylim([min(meanValues(:) - semValues(:)) * 0.9, max(maxPLV) * 1.1]);
% Set consistent y-axis limits
ylim([globalMin globalMax]);
yticks(globalMin:0.03:globalMax)

% Add legend
legend(uniqueConditions, 'Location', 'northwest', 'Interpreter', 'None', 'FontSize', 14);

% % Show grid
% grid on;
% set(ax, 'GridLineStyle', '--'); % Customize grid lines

% Release hold
hold off;

% Define the image name
figureName = ['meanPLV_' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '_over_all_chnPairs.jpg'];
figureNameFull = [outputPath '\' figureName];

% Save the figure as JPEG in high resolution
print(gcf, figureNameFull, '-djpeg', '-r300');

% save also in matlab format
figureName2 = ['meanPLV_' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '_over_all_chnPairs.fig'];
savefig([outputPath '\' figureName2]);

%% violin plot 10.05.2024
% Extract unique time periods and conditions
uniqueTimePeriods = unique(string(resultsTable.timePeriod), 'stable');
uniqueConditions = unique(string(resultsTable.condition), 'stable');
numTimePeriods = numel(uniqueTimePeriods);
numConditions = numel(uniqueConditions);

% initialize cell for data for violin plot
dataCell = cell(1, numTimePeriods * numConditions);
i = 1;

% Iterate over time periods
for timePeriodIdx = 1:numTimePeriods
    
    % Iterate over conditions
    for conditionIdx = 1:numConditions
        % Logical indices for the current time period and condition
        timePeriodIdx_logical = string(resultsTable.timePeriod) == uniqueTimePeriods(timePeriodIdx);
        conditionIdx_logical = string(resultsTable.condition) == uniqueConditions(conditionIdx);
        
        % Extract PLV values for the current time period and condition
        plvData = resultsTable.meanPLV(timePeriodIdx_logical & conditionIdx_logical);
        
        % put to a cell array
        dataCell{i} = plvData;
        i = i + 1;
    end
    
end

% Define face colors for same and different conditions
base_colors = [
    1 0.5 0.5; % same condition (light red)
    0.5 0.5 1; % different condition (light blue)
];
facecolors = repmat(base_colors, 4, 1);

% Define x-axis locations for each period and condition
x = [1 2 4 5 7 8 10 11];  
figure;
% Call the violin function with combined data and x-axis locations
% violin(dataCell, 'x', x, 'facecolor', facecolors) % default kernel bandwidth for smoothing
violin(dataCell, 'x', x, 'facecolor', facecolors, 'bw', 0.01); % kernel bandwidth = 0.01 less smoothing

% Customize the plot
ax = gca; % Get current axis
ax.XTick = [1.5, 4.5, 7.5, 10.5];
ax.XTickLabel = {'Baseline', 'First 1.9s Delay', 'Last 1.9s Delay', 'Recall'};
ax.FontSize = 12; % Increase font size for better readability
ax.LineWidth = 1.5; % Make the axes lines thicker
ax.XAxis.FontSize = 14; % Specific font size for X-axis
ax.YAxis.FontSize = 14; % Specific font size for Y-axis

% Add legend
% legend(uniqueConditions, 'Location', 'Best', 'Interpreter', 'None', 'FontSize', 10);

text(min(ax.XLim)*1.5, max(ax.YLim)*0.9, 'same', 'Color', [1 0.5 0.5], 'FontSize', 14) 
text(min(ax.XLim)*1.5, max(ax.YLim)*0.85, 'different', 'Color', [0.5 0.5 1], 'FontSize', 14) 
hold off;

ylabel(['Mean PLV in Freq Range: ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))  ' Hz'], ...
    'FontSize', 14);
xlabel('Time Period', 'FontSize', 14);
title(['Violin Plot of Mean PLV in ' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(1))) '-' num2str(first_delay_plv.PLVCond2.freq(isignificant_bins(end)))...
    'Hz in ' ROI1 '-' ROI2], 'FontSize', 14);
