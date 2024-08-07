function [chan_labels, ROI_labels, ROI_chanpairs] = select_chan_ROI_pair(table, ROI1, ROI2, patient_id, significant, patient_path)
% table - full filename of the table in which to search channels in ROI
% ROI1 - e.g. 'VTC'
% ROI2 - e.g. 'IPL'
% patient_id - e.g. 'p1883612'
% significant - if 1, select only significant chan pairs that were obtained by function PLV_permutation_stat for all trials
% patient_path - path for the data with significant PLV

if(~exist('table','var')) || isempty(table)
    %     table = 'E:\work\PhD\MemoryActions\results\iEEG\delayed condition\January2024_9pat\Response2XLS_CM_Memact_CHilbert_8-13Hz_AnyResp_OTP_-0.5-5.9_refBipo_Ep2024-01_encod+del_CHMult_2024-01-25_13-01-46.xls';
    %     table = 'F:\Sofia\MemoryActions\results\iEEG\delayed condition\January2024_9pat\Response2XLS_CM_Memact_CHilbert_8-13Hz_AnyResp_OTP_-0.5-5.9_refBipo_Ep2024-01_encod+del_CHMult_2024-01-25_13-01-46.xls';
    table = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\StructFind PAC_memact_all_chan_9pat_ROI_v2.xlsx';
end
if(~exist('significant','var')) || isempty(significant)
    significant = 0; % default all ROI1-ROI2 chan pairs
end

if significant == 0
    %     alphaIncreaseChan = readtable(table, 'ReadRowNames',true); % table with all channels from all patients showing alpha increase during the delay
    %
    %     % find channels for this patient and ROIs in the table
    %     ichanROI = find(strcmp(alphaIncreaseChan.pacient, patient_id) & ...
    %         alphaIncreaseChan.any_vs_bs == 1 &...
    %         (strcmp(alphaIncreaseChan.newROI, ROI1) | strcmp(alphaIncreaseChan.newROI, ROI2))); % chan indexes with alpha increase to any condition
    %     chan_labels = strrep(alphaIncreaseChan.name(ichanROI), [patient_id ' '], ''); % find labels of these chan
    %     ROI_labels = alphaIncreaseChan.newROI(ichanROI); % store also their ROI labels
    
    % find channels for this patient and ROIs in the table of all implanted channels
    PAC = readtable(table);
    ichanROI = find(strcmp(PAC.pacient, patient_id) & ...
        (strcmp(PAC.brainlabel, ROI1) | strcmp(PAC.brainlabel, ROI2))); % chan indexes of 2 ROIs
    chan_labels = PAC.name(ichanROI); % find labels of these chan
    ROI_labels = PAC.brainlabel(ichanROI); % store also their ROI labels
    
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
    % Specify the directory where the PLV data are stored
    PLV_path = [patient_path '\PLV_permut_stat\'];
    partialName = ['PLV_' ROI1 '-' ROI2 '_last 1.9s'];
    PLV_data_file = dir(fullfile(PLV_path, [partialName '*' 'bs_all_trials_200permut_2024-04.mat'])); % find the file
    
    % find significant ch pairs from all periods in the aggregated table
    aggregTable_path = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\'; % specify the directory where the aggregated data are stored
    aggregTableName = ['PLV_' ROI1 '-' ROI2 '_aggregated'];
    aggregTable = dir(fullfile(aggregTable_path, [aggregTableName '*' '.mat'])); % find the file
    
    if exist([PLV_path PLV_data_file.name],'file') == 2
        load([PLV_path PLV_data_file.name]);
        
        % channels in the data
        chan_labels = dataCond1.label;
        ROI_labels = {dataCond1.channelInfo.ROI}';
        
        if exist([aggregTable_path aggregTable.name],'file') == 2
            load([aggregTable_path aggregTable.name]);
            
            % Match the patient by name
            subjIndex = strcmp({aggregTable.patient}, patient_id);
            
            % Get indices of aggregated significant channel pairs for the current patient
            ROI_chanpairs = aggregTable(subjIndex).idxChan_in_Pairs;                        
        end
    else
        chan_labels = [];
        ROI_labels = [];
        ROI_chanpairs = [];
    end
    
end

