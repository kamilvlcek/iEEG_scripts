% --------------------
% a script to find the freq bins of significant PLV, which first
% - aggregates significant ch pairs from 3 periods (first 2s delay, last 2s delay, recall), 
% - plots scatterplot and histogram of significant ch pairs in each freq bin for all patients, 
% - finds the significant freq bins across all patients by binomial test 
% by Sofiia Moraresku
% May 2024
% --------------------

%% aggregate data from 3 periods (delay first 2s, delay last 2s, recall 0.5s)

% path for the output files
outputPath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\LMEM_May2024';

% load 3 tables with indices of all significant pairs across patients
% filenameTable1 = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_VTC-IPL_first 2s delay_vs_bs_all_trials_200permut_2024-03_summaryNewSignif.mat';
% filenameTable2 = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_VTC-IPL_last 2s delay_vs_bs_all_trials_200permut_2024-03_summaryNewSignif.mat';
% filenameTable3 = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_VTC-IPL_0.5s recall_vs_bs_all_trials_200permut_2024-04_summaryNewSignif.mat';
filenameTable1 = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_Hip-IPL_first 1.9s delay_vs_bs_all_trials_200permut_2024-04_summaryNewSignif.mat';
filenameTable2 = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_Hip-IPL_last 1.9s delay_vs_bs_all_trials_200permut_2024-04_summaryNewSignif.mat';
filenameTable3 = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_Hip-IPL_0.5s recall_vs_bs_all_trials_200permut_2024-04_summaryNewSignif.mat';
first_delay = load(filenameTable1);
last_delay = load(filenameTable2);
recall = load(filenameTable3);

% set up a path for the main data
% PLV_data1 =  'PLV_VTC-IPL_first 2s delay_vs_bs_all_trials_200permut_2024-03.mat';
% PLV_data2 =  'PLV_VTC-IPL_last 2s delay_vs_bs_all_trials_200permut_2024-03.mat';
% PLV_data3 =  'PLV_VTC-IPL_0.5s recall_vs_bs_all_trials_200permut_2024-04.mat';
PLV_data1 =  'PLV_Hip-IPL_first 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
PLV_data2 =  'PLV_Hip-IPL_last 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
PLV_data3 =  'PLV_Hip-IPL_0.5s recall_vs_bs_all_trials_200permut_2024-04.mat';
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % set up pacienti_memact by selecting which patients to analyze

% create an aggregated table
numPatients = numel(first_delay.tablePLV_allSubj); % number of patients in the structures
aggregTable = []; % aggregated table from 3 periods

for i=1:numPatients
    allSignificantidxChan = [first_delay.tablePLV_allSubj(i).idxChan_in_Pairs; last_delay.tablePLV_allSubj(i).idxChan_in_Pairs; recall.tablePLV_allSubj(i).idxChan_in_Pairs];
    unique_idxChan_in_Pairs = unique(allSignificantidxChan, 'rows');
    aggregTable(i).patient = first_delay.tablePLV_allSubj(i).patient;
    aggregTable(i).total_n_chnPairs = first_delay.tablePLV_allSubj(i).total_n_chnPairs;
    aggregTable(i).n_signif_chnPairs = size(unique_idxChan_in_Pairs,1);
    aggregTable(i).idxChan_in_Pairs = unique_idxChan_in_Pairs;
end   
    
% save an aggregated table
[filepath, name] = fileparts(filenameTable1);
parts = split(name, '_');
filename = [strjoin(parts(1:2), '_') '_aggregated_significant_chnPairs'];
full_path = fullfile(filepath, filename);
save(full_path, 'aggregTable')

% Extract important info from the file name
ROI1 = regexp(PLV_data1, '(?<=PLV_)[^-]+', 'match', 'once');
ROI2 = regexp(PLV_data1, '(?<=-)[^-_]+', 'match', 'once');

%% plot significant difference for all signif pairs - individual points with real PLV difference between any period and bs

% Create a new figure
figure('Name','significant PLV difference between any period and baseline');

% Define a colormap for the scatter plot. This is a simple way to ensure each
% channel pair gets a unique color, but with many pairs, some colors may be similar.
% colors = lines(max([tablePLV_allSubj.n_signif_chnPairs]));  % 'lines' is a MATLAB colormap with distinct colors
colors = distinguishable_colors(max([aggregTable.n_signif_chnPairs]));

% Iterate over each subject
for p = 1:numel(pacienti)
    if ~pacienti(p).todo
        continue; % Skip if the subject is not marked 'todo'
    end
    
    % Match the patient by name
    subjIndex = find(strcmp({aggregTable.patient}, pacienti(p).folder));
    if isempty(subjIndex)
        continue;
    end
    
    filePath = [basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\'];
    %filePath = [basedir pacienti(p).folder '/' subfolder '/' PLV_folder '/'];

    if isfile([filePath PLV_data1])
        % Load the PLV data for all 3 periods
        first_delay_plv = load([filePath PLV_data1]);
        last_delay_plv = load([filePath PLV_data2]);
        recall_plv = load([filePath PLV_data3]);
        
        % Get indices of aggragated significant channel pairs for the current subject
        significantidxChan = aggregTable(subjIndex).idxChan_in_Pairs;
        isignPair = ismember(first_delay_plv.ROI_chanpairs, significantidxChan, 'rows');
        significantPairsIndicesForSubject = find(isignPair);
             
        % Create subplot for the current patient
        subplot(ceil(numPatients/2),2, subjIndex);
        hold on
        % Loop over each significant channel pair for the current patient
        for isignChPair = 1:numel(significantPairsIndicesForSubject)
%             signFreqIndices = find(dataFiltered(isignChPair, :) > 0); % significant freq bins indices
%             freqsToPlot = PLVCond2.freq(signFreqIndices);
%             magsToPlot = dataFiltered(isignChPair, signFreqIndices);
            signFreqIndices_first_delay = find(first_delay_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject(isignChPair), :)>0);
            signFreqIndices_last_delay = find(last_delay_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject(isignChPair), :)>0);
            signFreqIndices_recall = find(recall_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject(isignChPair), :)>0);
            
             PLV_diff_first_delay = first_delay_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject(isignChPair), signFreqIndices_first_delay);
             PLV_diff_last_delay = last_delay_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject(isignChPair), signFreqIndices_last_delay);
             PLV_diff_last_recall = recall_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject(isignChPair), signFreqIndices_recall);   
             
             magsToPlot = [PLV_diff_first_delay, PLV_diff_last_delay, PLV_diff_last_recall];            
             freqsToPlot = first_delay_plv.PLVCond2.freq([signFreqIndices_first_delay, signFreqIndices_last_delay, signFreqIndices_recall]);
                         
            % Choose color for the current channel pair
            currentColor = colors(mod(isignChPair, size(colors, 1)) + 1, :);            
            
            % Plot the dots for the current channel pair
            scatter(freqsToPlot, magsToPlot, 18, currentColor, 'filled');
        end
        
        % Customize the subplot
        title([pacienti(p).folder ', ' num2str(numel(significantPairsIndicesForSubject)) ' chn pairs'], 'Interpreter', 'None');
        xlabel('Frequency (Hz)');
        ylabel('PLV Difference');
        xlim([first_delay_plv.PLVCond2.freq(1) first_delay_plv.PLVCond2.freq(end)]);                
        hold off;
    end
end

% Get handles to all axes objects
allAxes = findobj(gcf, 'type', 'axes');

% Find the minimum and maximum y-axis limits across all subplots
minY = min(arrayfun(@(x) x.YLim(1), allAxes));
maxY = max(arrayfun(@(x) x.YLim(2), allAxes));

% Set the same y-axis limits for all subplots
arrayfun(@(x) set(x, 'YLim', [minY, maxY]), allAxes);

%% plot histogram with ratio of significant pairs in each freq bin

% Create a new figure
figure('Name','ratio of significant channel pairs across freq bins');

ratioPairsAllPatients = zeros(numPatients, numel(first_delay_plv.PLVCond2.freq)); % ratio of significant pairs in each freq bin in all patients, patients x freq bins
countsSigPRatioAll = zeros(numPatients, numel(first_delay_plv.PLVCond2.freq)); % number of significant pairs in each freq bin in all patients, patients x freq bins

% Iterate over each subject
for p = 1:numel(pacienti)
    if ~pacienti(p).todo
        continue; % Skip if the subject is not marked 'todo'
    end
    
    % Match the patient by name
    subjIndex = find(strcmp({aggregTable.patient}, pacienti(p).folder));
    if isempty(subjIndex)
        continue;
    end

    %filePath = [basedir pacienti(p).folder '/' subfolder '/' PLV_folder '/'];
    filePath = [basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\'];
    
    if isfile([filePath PLV_data1])
        % Load the PLV data for all 3 periods
        first_delay_plv = load([filePath PLV_data1]);
        last_delay_plv = load([filePath PLV_data2]);
        recall_plv = load([filePath PLV_data3]);
        
        % Get indices of aggragated significant channel pairs for the current subject
        significantidxChan = aggregTable(subjIndex).idxChan_in_Pairs;
        isignPair = ismember(first_delay_plv.ROI_chanpairs, significantidxChan, 'rows');
        significantPairsIndicesForSubject = find(isignPair);
        
        % Filter the data to include only the significant channel pairs in three periods
        dataFiltered1 = first_delay_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject, :);
        dataFiltered2 = last_delay_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject, :);
        dataFiltered3 = recall_plv.plv_signif_allPairs_clustcorr(significantPairsIndicesForSubject, :);
        dataFilteredAll = dataFiltered1 + dataFiltered2 + dataFiltered3;

        % Create subplot for the current patient
        subplot(ceil(numPatients/2),2, subjIndex);
        
        % Sum the number of positive values across rows for each frequency bin
        positiveCounts = sum(dataFilteredAll > 0, 1);
%         positiveCounts1 = sum(dataFiltered1 > 0, 1);
%         positiveCounts2 = sum(dataFiltered2 > 0, 1);
%         positiveCounts3 = sum(dataFiltered3 > 0, 1);
%         positiveCounts = positiveCounts1 + positiveCounts2 + positiveCounts3;
        countsSigPRatioAll(subjIndex,:) = positiveCounts;

        % normalize by total number of ch pairs
        ratioPairs = positiveCounts/aggregTable(subjIndex).total_n_chnPairs;
        ratioPairsAllPatients(subjIndex,:) = ratioPairs;
        
        % Plot histogram
        bar(first_delay_plv.PLVCond2.freq, ratioPairs, 'FaceColor', 'b');

        % Add labels and title
        xlabel('Frequency (Hz)');
        ylabel('Sig.P.Ratio');
        title([pacienti(p).folder ', ' num2str(size(dataFiltered1, 1)) '/' num2str(aggregTable(subjIndex).total_n_chnPairs) ' chn pairs'], 'Interpreter', 'None');
    end
end


% Get handles to all axes objects
allAxes = findobj(gcf, 'type', 'axes');

% Find the minimum and maximum y-axis limits across all subplots
minY = min(arrayfun(@(x) x.YLim(1), allAxes));
maxY = max(arrayfun(@(x) x.YLim(2), allAxes));

% Set the same y-axis limits for all subplots
arrayfun(@(x) set(x, 'YLim', [minY, maxY]), allAxes);

% % normalize by total number of ch pairs in all patients
% proportionsTotal = sum(ratioPairsAllPatients,1); % proportions for all freq bins - sum across all patients
% 
% % Plot histogram for mean of percentage of significant pairs across all patients 
% subplot(ceil(numPatients/2),2, subjIndex+1);
% bar(first_delay_plv.PLVCond2.freq, proportionsTotal, 'FaceColor', 'r');
% 
% chance_level = median(proportionsTotal);  % The median Sig.P.Ratio across all freq bins is used as the chance level for binomial testing similar to Park et al. 2019
% hold on
% line(xlim, [chance_level chance_level], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2); % plot median
% 
% % Add labels and title
% xlabel('Frequency (Hz)');
% ylabel('Sig.P.Ratio');
% title(['Sum across all ' num2str(numPatients) ' patients, ' num2str(sum([aggregTable.n_signif_chnPairs])) '/' num2str(sum([aggregTable.total_n_chnPairs])) ' chn pairs']);
% legend('', 'median across all freq bins');


%% binomial test 1 version - probably less correct than the second approach
% proportionsTotal = sum(ratioPairsAllPatients,1); % proportions for all freq bins - sum across all patients
% chnPairsTotal = sum([aggregTable.total_n_chnPairs]); % total number of channel pairs - sum across all patients
% chance_level = median(proportionsTotal);  % The median Sig.P.Ratio across all freq bins is used as the chance level for binomial testing similar to Park et al. 2019
% 
% x = round(proportionsTotal*chnPairsTotal); % observed number of 'successes'
% 
% % Perform binomial test for each proportion - freq bin
% p_values = zeros(size(proportionsTotal));
% for i = 1:numel(proportionsTotal)
%     p_values(i) = binocdf(x(i), chnPairsTotal , chance_level, 'upper' ); % one-sided test
% end
% 
% % FDR correction 
% [~, ~, adj_p_values] = fdr_bh(p_values, 0.05, 'dep');
% 
% isignificant_bins = find(adj_p_values<0.05);
% significant_bins = first_delay_plv.PLVCond2.freq(isignificant_bins);
% fprintf('Significant Frequency Bins in Hz:\n');
% disp(significant_bins);

%% binomial test 2 version
chnPairsTotal = sum([aggregTable.total_n_chnPairs]); % total number of channel pairs - sum across all patients
SigPRatioAll = sum(countsSigPRatioAll,1)/chnPairsTotal ; % proportions for all freq bins across all subejcts
chance_level = median(SigPRatioAll);  % The median Sig.P.Ratio across all freq bins is used as the chance level for binomial testing similar to Park et al. 2019
x = sum(countsSigPRatioAll,1); % observed number of 'successes'
% Perform binomial test for each proportion - freq bin
p_values = zeros(size(SigPRatioAll ));
for i = 1:numel(SigPRatioAll )
    p_values(i) = binocdf(x(i), chnPairsTotal, chance_level, 'upper' ); % one-sided test
end

isignificant_bins_uncorr = find(p_values<0.05); % indexes of significant freq bins before fdr corr

% FDR correction - do we need it here?
[~, ~, adj_p_values] = fdr_bh(p_values, 0.05, 'dep');

isignificant_bins = find(adj_p_values<0.05);
significant_bins = first_delay_plv.PLVCond2.freq(isignificant_bins);
fprintf('Significant Frequency Bins in Hz:\n');
disp(significant_bins);

% new histogram
figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]); % Maximize the figure window
barHandle = bar(first_delay_plv.PLVCond2.freq, SigPRatioAll, 'FaceColor', 'r');
hold on
linehandle = line(xlim, [chance_level chance_level], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2); % plot median

% Add labels and title
xlabel('Frequency (Hz)');
ylabel('Sig.P.Ratio');
title(['Sig.P.Ratio across all ' num2str(numPatients) ' patients, in total: ' num2str(sum([aggregTable.n_signif_chnPairs])) '/' num2str(sum([aggregTable.total_n_chnPairs])) ' chn pairs significant in ' ROI1 '-' ROI2]);
%legend('', 'median across all freq bins');

% Get the x and y values for the tops of the bars
x = barHandle(1).XData;
y = barHandle(1).YData;
ylim([0, max(y)* 1.1]);

% Loop over the significant bars to place asterisks (uncorrected)
hold on; % Keep the current bar
for i = 1:numel(isignificant_bins_uncorr)
    text(x(isignificant_bins_uncorr(i)), y(isignificant_bins_uncorr(i)) + (max(y) - min(y)) * 0.025, '*', 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 16);
end

% Loop over the significant bars to place asterisks (fdr-corrected)
hold on; % Keep the current bar
for i = 1:numel(isignificant_bins)
    text(x(isignificant_bins(i)), y(isignificant_bins(i)) + (max(y) - min(y)) * 0.05, '**', 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 16);
end

% Create an invisible plot for the legend entry
hInvisible = plot(nan, nan, 'w'); % 'w' makes it white/invisible
hInvisible1 = plot(nan, nan, 'w'); % 'w' makes it white/invisible

% Add the legend entry for the asterisks
legend([linehandle, hInvisible, hInvisible1], {'median across all freq bins', ' *  - p < 0.05 uncorr, binomial test', ' **  - p < 0.05 fdr-corr, binomial test'}, 'Location', 'best');
hold off;

% Define the image name
figureName = ['SigPairRatio_of_PLV_difference_between_anyPeriod_and_bs_in_' ROI1 '-' ROI2 '.jpg'];
filename2 = [outputPath '\' figureName];

% Save the figure as JPEG in high resolution
print(gcf, filename2, '-djpeg', '-r300'); 

