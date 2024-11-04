% find all significant ch pairs with Granger across all patients 
% and save a summary table with their numbers and indices

%% set up a path for data
% granger_data = 'Granger_VTC-IPL last 2s delay all_trials_2024-02.mat';
granger_data = 'Granger_VTC-IPL last 2s delay all_trials_2024-03.mat';
granger_folder = 'Granger_permut_stat';

setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact();

%%
tableGranger_allSubj = [];

ip = 1; % index of patients with PLV data

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if isfile([basedir pacienti(p).folder '\' subfolder '\' granger_folder '\' granger_data])
            % load the granger data with statistics
            load([basedir pacienti(p).folder '\' subfolder '\' granger_folder '\' granger_data]);
            
            % first find indexes of chan pairs with significant direction VTC -> IPL (negative difference)
%             significant_chanPairs = find(sum(Granger_signif_allPairs_clustcorr,2) < 0);
            significant_chanPairs = find(any(Granger_signif_allPairs_clustcorr < 0,2) &...
                ~any(Granger_signif_allPairs_clustcorr > 0, 2)); % at least one negative value and no pisitive values
            idxChan = ROI_chanpairs(significant_chanPairs, :); % channel indices in each pair            
            chn_labels = [{dataDelay.label{idxChan(:,1)}}' {dataDelay.label{idxChan(:,2)}}'];  % original labels of chan in pairs
            
            tableGranger_allSubj(ip).patient = pacienti(p).folder;
            tableGranger_allSubj(ip).total_n_chnPairs = size(Granger_signif_allPairs_clustcorr,1);
            tableGranger_allSubj(ip).n_signif_chnPairs = numel(significant_chanPairs);
            tableGranger_allSubj(ip).idxChan_in_Pairs = idxChan;
            tableGranger_allSubj(ip).chn_labels_in_Pairs = chn_labels;
            ip = ip+1;
            
        end
    end
end

%% save
[~, name] = fileparts(granger_data);
filename = [name '_summary.mat'];
% filepath = 'E:\work\PhD\MemoryActions\results\iEEG\connectivity\group data';
filepath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data';
full_path = fullfile(filepath, filename);
save(full_path, 'tableGranger_allSubj')

% filename2 = [name '_summary.xlsx']; % export to xls
% full_path2 = fullfile(filepath, filename2);
% writetable(struct2table(tableGranger_allSubj), full_path2) 

