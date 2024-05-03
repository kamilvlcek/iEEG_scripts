% --------------------
% a script to run a Linear mixed effects model (LMEM) on mean PLV for 2 conditions - same and diff - fixed effects, 
% while accounting for random effects - ch pair and patientID
% by Sofiia Moraresku
% May 2024
% --------------------

% first, set up the indexes of freq bins to use for averaging PLV in the model (obtained earlier by plot_PLVdiff_variability)
% isignificant_bins = [1 2 3]; % delta-theta, 2-4Hz
isignificant_bins = [6 7]; % alpha, 7-8Hz
% isignificant_bins = [1:8]; % theta-alpha, 2-9Hz

% load an aggregated table with indices of all significant pairs for all three periods across patients
filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_Hip-IPL_aggregated_significant_chnPairs.mat';
load(filenameTable);

% path for the output files
outputPath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\LMEM_May2024\conditions';

% set up a path for the main data
PLV_data =  'PLV_Hip-IPL_diff_vs_same_first 1.9s delay__trials_2024-05.mat';
% PLV_data =  'PLV_Hip-IPL_diff_vs_same_0.5s recall__trials_2024-05.mat';
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % set up pacienti_memact by selecting which patients to analyze

% Extract important info from the file name
ROI1 = regexp(PLV_data, '(?<=PLV_)[^-]+', 'match', 'once');
ROI2 = regexp(PLV_data, '(?<=-)[^-_]+', 'match', 'once');
period = regexp(PLV_data, '(first|last) \d+\.\ds delay|(?<=_)\w+(?=__trials)', 'match', 'once');
% period = '0.5s recall';

%% get average PLV over the selected freq range for two conditions
% Initialize an empty table for storing the results
resultsTable = table([], [], [], [], 'VariableNames', {'patientId', 'chPair', 'condition', 'meanPLV'});

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

    if isfile([filePath PLV_data])
        % Load the file containing PLV data for two conditions
        load([filePath PLV_data]);

        % Get indices of aggregated significant channel pairs for the current subject
        significantidxChan = aggregTable(subjIndex).idxChan_in_Pairs; % n signif pairs x 2 (each pair with chan indices)

        % For each significant pair extract raw PLV for each condition
        for i = 1:size(significantidxChan,1)
            ipair = significantidxChan(i,:);

            % Calculate the mean PLV within the selected frequency range for this pair for all conditions
            meanPLV_same = mean(squeeze(PLVCond1.plvspctrm(ipair(1),ipair(2),isignificant_bins)));
            meanPLV_diff = mean(squeeze(PLVCond2.plvspctrm(ipair(1),ipair(2),isignificant_bins))); 
            
            chn_labels = [PLVCond1.label{ipair(1)} PLVCond1.label{ipair(2)}];  % original labels of chan in this pair

            % Append the results to the table
            newRow1 = {pacienti(p).folder, chn_labels, 'same', meanPLV_same}; 
            newRow2 = {pacienti(p).folder, chn_labels, 'diff', meanPLV_diff}; 
            resultsTable = [resultsTable; newRow1; newRow2];
        end
    end
end

resultsTable.condition = categorical(resultsTable.condition);
resultsTable.patientId = categorical(resultsTable.patientId);
resultsTable.chPair = categorical(resultsTable.chPair);

%% LMM
% Fit a linear mixed-effects model: condition - fixed effect, patient and chan pair - random effects
model = fitlme(resultsTable, 'meanPLV ~ condition + (1 | patientId) + (1 | patientId:chPair)'); % chapter 36.4.2 of iEEG book

% Display the model summary
disp(model); % the output can be copied to txt file

%% export the results to the xls file

% Define the Excel file name
outputname = ['LMEM_output_meanPLV_' num2str(PLVCond2.freq(isignificant_bins(1))) '-' num2str(PLVCond2.freq(isignificant_bins(end))) 'Hz_' ROI1 '-' ROI2 '_' period '.xlsx'];
filename = [outputPath '\' outputname];

% Manually create a table from the 'titleddataset'
coefficients = model.Coefficients;
coefTable = table(coefficients.Name, coefficients.Estimate, coefficients.SE, coefficients.tStat, ...
                  coefficients.DF, coefficients.pValue, coefficients.Lower, coefficients.Upper, ...
                  'VariableNames', {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'pValue', 'Lower', 'Upper'});

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

%% visualize the results
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
