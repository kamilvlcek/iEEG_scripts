% find all significant ch pairs with PLV across all patients
% and save a summary table with their numbers and indices
 
%% set up a path for data
% PLV_data =  'PLV_VTC-IPL_last 2s delay_vs_bs_all_trials_2024-02.mat';
PLV_data =  'PLV_VTC-IPL_first 2s delay_vs_bs_all_trials_2024-03.mat';
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact();

%%
tablePLV_allSubj = [];

ip = 1; % index of patients with PLV data

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if isfile([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data])
            % load the PLV data with statistics
            load([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data]);
            
            %             % Find rows (ch pairs) with at least 2 positive values (2 freq bins)
            %             positive_values_count = sum(plv_signif_allPairs_clustcorr > 0, 2);
            %             significant_chanPairs = find(positive_values_count >= 2);
            
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
filename = [name '_summary.mat'];
% filepath = 'E:\work\PhD\MemoryActions\results\iEEG\connectivity\group data';
filepath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data';
full_path = fullfile(filepath, filename);
save(full_path, 'tablePLV_allSubj')

% filename2 = [name '_summary.xlsx']; % export to xls
% full_path2 = fullfile(filepath, filename2);
% writetable(struct2table(tablePLV_allSubj), full_path2)

