% find all significant ch pairs with PLV across all patients
% and save a summary table with their numbers and indices in .mat file

%% set up a path for data
PLV_data =  'PLV_VTC-IPL_0.5s recall_vs_bs_all_trials_200permut_2024-04.mat';
% PLV_data =  'PLV_VTC-IPL_last 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
% PLV_data =  'PLV_VTC-IPL_first 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
% PLV_data =  'PLV_Hip-IPL_0.5s recall_vs_bs_all_trials_200permut_2024-04.mat';
% PLV_data =  'PLV_Hip-IPL_first 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
% PLV_data =  'PLV_Hip-IPL_last 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
PLV_folder = 'PLV_permut_stat';
direction = 1;  % 1 - significant ch pairs where cond 1 > cond 2; 2 - significant ch pairs where cond 2 > cond 1

setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % 

%%
tablePLV_allSubj = [];

ip = 1; % index of patients with PLV data

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if isfile([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data])
            % load the PLV data with statistics
            load([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data]);
            
            if direction == 1 % cond 1 > cond 2
                % find all rows with at least one positive value
                positives = plv_signif_allPairs_clustcorr > 0;
                % Sum across columns to find rows with at least one positive
                rowsWithPositives = sum(positives, 2) > 0;
                % rows with negative values
                negatives = plv_signif_allPairs_clustcorr < 0;
                % Sum across columns to ensure no negatives in the row
                rowsWithoutNegatives = sum(negatives, 2) == 0;
                % Combine conditions: rows with at least one positive and no negatives
                desiredRows = rowsWithPositives & rowsWithoutNegatives;
                
            elseif direction == 2 % cond 2 > cond 1
                % find all rows with at least one negative value
                negatives = plv_signif_allPairs_clustcorr < 0;
                % Sum across columns to find rows with at least one negative
                rowsWithNegatives = sum(negatives, 2) > 0;
                % rows with positive values
                positives = plv_signif_allPairs_clustcorr > 0;
                % Sum across columns to ensure no positives in the row
                rowsWithoutPositives = sum(positives, 2) == 0;
                % Combine conditions: rows with at least one positive and no negatives
                desiredRows = rowsWithNegatives & rowsWithoutPositives;
            end
            
            % Find indices of these rows
            significant_chanPairs = find(desiredRows);
            
            idxChan = ROI_chanpairs(significant_chanPairs, :); % channel indices in each pair            
            chn_labels = [{PLVCond1.label{idxChan(:,1)}}' {PLVCond1.label{idxChan(:,2)}}'];  % original labels of chan in pairs
            tablePLV_allSubj(ip).patient = pacienti(p).folder;
            
            % number of channels in each ROI
            unique_ROI = unique({dataCond1.channelInfo.ROI}, 'stable'); % names of ROIs
            ROI1_fieldname = [unique_ROI{1} '_n_channels'];
            nROI1 = sum(strcmp({dataCond1.channelInfo.ROI},unique_ROI{1}));
            ROI2_fieldname = [unique_ROI{2} '_n_channels'];
            nROI2 = sum(strcmp({dataCond1.channelInfo.ROI},unique_ROI{2}));
            tablePLV_allSubj(ip).(ROI1_fieldname) = nROI1;
            tablePLV_allSubj(ip).(ROI2_fieldname) = nROI2;
            
            tablePLV_allSubj(ip).total_n_chnPairs = size(plv_signif_allPairs_clustcorr,1);
            tablePLV_allSubj(ip).n_signif_chnPairs = numel(significant_chanPairs);
            tablePLV_allSubj(ip).idxChan_in_Pairs = idxChan;
            tablePLV_allSubj(ip).chn_labels_in_Pairs = chn_labels;
            ip = ip+1;
            
        end
    end
end

%% save
[~, name] = fileparts(PLV_data);
filename = [name '_summaryNewSignif.mat'];
% filepath = 'E:\work\PhD\MemoryActions\results\iEEG\connectivity\group data';
filepath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data';
full_path = fullfile(filepath, filename);
save(full_path, 'tablePLV_allSubj')

% filename2 = [name '_summary.xlsx']; % export to xls
% full_path2 = fullfile(filepath, filename2);
% writetable(struct2table(tablePLV_allSubj), full_path2)

