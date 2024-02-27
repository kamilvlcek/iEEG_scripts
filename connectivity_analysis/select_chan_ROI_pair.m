function [chan_labels, ROI_labels, ROI_chanpairs] = select_chan_ROI_pair(table, ROI1, ROI2, patient_id, significant, patient_path)
% table - full filename of the table in which to search channels in ROI
% ROI1 - e.g. 'VTC'
% ROI2 - e.g. 'IPL'
% patient_id - e.g. 'p1883612'
% significant - if 1, select only significant chan pairs that were obtained by function PLV_permutation_stat for all trials
% patient_path - path for the data with significant PLV

if(~exist('table','var')) || isempty(table)
    table = 'E:\work\PhD\MemoryActions\results\iEEG\delayed condition\January2024_9pat\Response2XLS_CM_Memact_CHilbert_8-13Hz_AnyResp_OTP_-0.5-5.9_refBipo_Ep2024-01_encod+del_CHMult_2024-01-25_13-01-46.xls';
end
if(~exist('significant','var')) || isempty(significant)
    significant = 0; % default all ROI1-ROI2 chan pairs
end

if significant == 0
    alphaIncreaseChan = readtable(table, 'ReadRowNames',true); % table with all channels from all patients showing alpha increase during the delay
    
    % find channels for this patient and ROIs in the table
    ichanROI = find(strcmp(alphaIncreaseChan.pacient, patient_id) & ...
        alphaIncreaseChan.any_vs_bs == 1 &...
        (strcmp(alphaIncreaseChan.newROI, ROI1) | strcmp(alphaIncreaseChan.newROI, ROI2))); % chan indexes with alpha increase to any condition
    chan_labels = strrep(alphaIncreaseChan.name(ichanROI), [patient_id ' '], ''); % find labels of these chan
    ROI_labels = alphaIncreaseChan.newROI(ichanROI); % store also their ROI labels
    
    % Channel pairs to run for connectivity
    ROI1_chan_idx = find(contains(ROI_labels, ROI1)); % channel indexes for ROI1
    ROI2_chan_idx = find(contains(ROI_labels, ROI2)); % channel indexes for ROI2
    
    % all possible pairs between 2 ROIs
    ROI_chanpairs = zeros(length(ROI1_chan_idx) * length(ROI2_chan_idx), 2);
    count = 1;
    for i = 1:length(ROI1_chan_idx)
        for j = 1:length(ROI2_chan_idx)
            ROI_chanpairs(count, :) = [ROI1_chan_idx(i) ROI2_chan_idx(j)];
            count = count + 1;
        end
    end
    
else
    
    % Specify the directory where the data are stored
    PLV_path = [patient_path '\PLV_permut_stat\'];
    partialName = ['PLV_' ROI1 '-' ROI2 '_last 2s'];
    PLV_data_file = dir(fullfile(PLV_path, [partialName '*' 'bs_all_trials_2024-02.mat'])); % find the file
    
    if exist([PLV_path PLV_data_file.name],'file') == 2
        load([PLV_path PLV_data_file.name]);
        
        % channels in the data
        chan_labels = dataCond1.label;
        ROI_labels = {dataCond1.channelInfo.ROI}';
        
        %     % find indexes of chan pairs with significant PLV difference, only positive delay > bs
        %     significant_chanPairs = sum(plv_signif_allPairs_clustcorr,2) > 0;
        %     ROI_chanpairs = ROI_chanpairs(significant_chanPairs, :);
        
        %     % find rows (ch pairs) with at least 2 positive values (2 freq bins) delay > bs (more strict selection of pairs)
        %     positive_values_count = sum(plv_signif_allPairs_clustcorr > 0, 2);
        %     significant_chanPairs = find(positive_values_count >= 2);
        
        % Find rows with at least 2 consecutive positive values
        rows_with_consecutive_positives = [];
        
        for i = 1:size(plv_signif_allPairs_clustcorr, 1)
            row = plv_signif_allPairs_clustcorr(i, :);
            % Use a sliding window to check for consecutive positive values
            for j = 1:(length(row)-1)
                if row(j) > 0 && row(j+1) > 0
                    rows_with_consecutive_positives = [rows_with_consecutive_positives, i];
                    break; % Break once you find the first occurrence
                end
            end
        end
        
        significant_chanPairs = rows_with_consecutive_positives';
        ROI_chanpairs = ROI_chanpairs(significant_chanPairs, :);
    else
        chan_labels = [];
        ROI_labels = [];
        ROI_chanpairs = [];
    end
end

