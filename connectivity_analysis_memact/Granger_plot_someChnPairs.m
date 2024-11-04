%%% plot some examples of Granger spectrum for individual channel pairs

% Load Granger data for both conditions (same and diff)
% same = load('D:\eeg\motol\pacienti\p2179801 Syk VT62\memact\Granger_permut_stat\Granger_VTC-IPL_last 1.9s delay_same_trials_2024-09.mat');
% diff = load('D:\eeg\motol\pacienti\p2179801 Syk VT62\memact\Granger_permut_stat\Granger_VTC-IPL_last 1.9s delay_diff_trials_2024-09.mat');

same = load('D:\eeg\motol\pacienti\p2179801 Syk VT62\memact\Granger_permut_stat\Granger_Hip-IPL_first 1.9s delay_same_trials_2024-09.mat');
diff = load('D:\eeg\motol\pacienti\p2179801 Syk VT62\memact\Granger_permut_stat\Granger_Hip-IPL_first 1.9s delay_diff_trials_2024-09.mat');

% same = load('d:\eeg\motol\pacienti\p1883612 And VT59\memact\Granger_permut_stat\Granger_Hip-IPL_first 1.9s delay_same_trials_2024-09.mat');
% diff = load('d:\eeg\motol\pacienti\p1883612 And VT59\memact\Granger_permut_stat\Granger_Hip-IPL_first 1.9s delay_diff_trials_2024-09.mat');


% Select channel pairs to plot
chanPairsToPlot = 1:size(same.ROI_chanpairs,1); % All channel pairs in the data
iPairsGranger = 1:2:2*size(same.ROI_chanpairs, 1); % As we have 2 rows of values for each channel pair (two directions)

% Define colors for the plot (Same: Red shades, Diff: Blue shades)
colors_same = [1 0.5 0.5; 1 0 0]; % Light red for Dir1, Dark red for Dir2
colors_diff = [0.5 0.5 1; 0 0 1]; % Light blue for Dir1, Dark blue for Dir2

% Number of channel pairs to plot
numPairsToPlot = 16; % Change this to however many you want to plot

for i = 1:numPairsToPlot
    
    % Get the channel pair index and the corresponding Granger indices
    ipair = same.ROI_chanpairs(chanPairsToPlot(i), :); % Real indexes of channels in the data
    ipairGranger = iPairsGranger(chanPairsToPlot(i)); % Index of channel pair in Granger spectrum

    % Create a new figure for each channel pair
    figure(i);
    hold on;
    
    % Plot Granger data for the "same" condition
    plot(same.grangerPeriod.freq, same.grangerPeriod.grangerspctrm(ipairGranger,:), 'Color', colors_same(1,:), 'LineWidth', 2, 'DisplayName', 'Same Condition (Dir1)');
    plot(same.grangerPeriod.freq, same.grangerPeriod.grangerspctrm(ipairGranger+1,:), 'Color', colors_same(2,:), 'LineWidth', 2, 'DisplayName', 'Same Condition (Dir2)');
    
    % Plot Granger data for the "diff" condition
    plot(diff.grangerPeriod.freq, diff.grangerPeriod.grangerspctrm(ipairGranger,:), 'Color', colors_diff(1,:), 'LineWidth', 2, 'DisplayName', 'Diff Condition (Dir1)');
    plot(diff.grangerPeriod.freq, diff.grangerPeriod.grangerspctrm(ipairGranger+1,:), 'Color', colors_diff(2,:), 'LineWidth', 2, 'DisplayName', 'Diff Condition (Dir2)');
    
    % Customize the plot
    ylim([0 0.1]); % Adjust as needed
    yticks(0:0.02:0.1); % Set y-ticks with step 0.02
    xlim([1 20]);
    xlabel('Frequency (Hz)');
    ylabel('Granger Causality');
    
    % Get channel labels for the legend
    chlabel1 = split(same.grangerPeriod.labelcmb{ipairGranger,1}, '[');
    chlabel2 = split(same.grangerPeriod.labelcmb{ipairGranger+1,1}, '[');
    
    % Add a title to each plot
    title(['Spectral Granger between ' chlabel1{1} ' (' same.dataPeriod.channelInfo(ipair(2)).ROI ') and ' chlabel2{1} ' (' same.dataPeriod.channelInfo(ipair(1)).ROI ')']);
    
    % Add a legend with meaningful labels
    legend(['Same from ' chlabel1{1} same.dataPeriod.channelInfo(ipair(2)).ROI ' to ' chlabel2{1} same.dataPeriod.channelInfo(ipair(1)).ROI], ['Same from ' chlabel2{1} same.dataPeriod.channelInfo(ipair(1)).ROI ' to ' chlabel1{1} same.dataPeriod.channelInfo(ipair(2)).ROI],...
        ['Diff from ' chlabel1{1} same.dataPeriod.channelInfo(ipair(2)).ROI ' to ' chlabel2{1} same.dataPeriod.channelInfo(ipair(1)).ROI], ['Diff from ' chlabel2{1} same.dataPeriod.channelInfo(ipair(1)).ROI ' to ' chlabel1{1} same.dataPeriod.channelInfo(ipair(2)).ROI], 'Location', 'Best', 'Interpreter', 'None');
    
    % Set the figure properties
    set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 20, 15]); % Adjust the figure size as needed
    
    % save the figure
    % saveas(gcf, ['Granger_ChannelPair_' num2str(i) '.png']);
end

