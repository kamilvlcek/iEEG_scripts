%% script creates a nod file for BrainNet Viewer with coordinates of all channels in LOC and IPL (or any other ROI)
% The node file is defined as an ASCII text file with the suffix ‘node’. In the node file, there are 6 columns: 
% columns 1-3 represent node coordinates, column 4 represents node colors, column 5 represents node sizes, and the last column represents node labels. 
% Please note, a symbol ‘-‘ (no ‘’) in column 6 means no labels and blank characters in label would cause an error

% to be configured for each pair of ROIs you want to visualize
ROI1 = 'VTC';
ROI2 = 'IPL';

%% table with all channels from all patients showing an alpha increase during the delay
% alphaIncrease = 'E:\work\PhD\MemoryActions\results\iEEG\delayed condition\January2024_9pat\Response2XLS_CM_Memact_CHilbert_8-13Hz_AnyResp_OTP_-0.5-5.9_refBipo_Ep2024-01_encod+del_CHMult_2024-01-25_13-01-46.xls';
% alphaIncreaseChan = readtable(alphaIncrease, 'ReadRowNames',true);
% 
% % find all channels with an alpha increase in ROI1 and ROI2 for patients that have both of them
% % ROI_chan = find((strcmp(alphaIncreaseChan.mybrainlabel, ROI1) | strcmp(alphaIncreaseChan.mybrainlabel, ROI2))...
% %     & ~(strcmp(alphaIncreaseChan.pacient, 'p65639') | strcmp(alphaIncreaseChan.pacient, 'p1883612_VT66') )); % exclude 2 patients that have only one ROI
% ROI_chan = find(strcmp(alphaIncreaseChan.newROI, 'AnG') | strcmp(alphaIncreaseChan.newROI, 'SMG')...
%     | strcmp(alphaIncreaseChan.newROI, 'ITG') | strcmp(alphaIncreaseChan.newROI, 'VTC') | strcmp(alphaIncreaseChan.newROI, 'MTL')...
%     | strcmp(alphaIncreaseChan.newROI, 'precun'));

%% table with all implanted channels from all patients
chan_filepath = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\StructFind PAC_memact_all_chan_9pat_ROI_v2.xlsx';
allChan = readtable(chan_filepath);

% find all channels in ROI1 and ROI2 for patients that have both of them
ROI_chan = find((strcmp(allChan.brainlabel, ROI1) | strcmp(allChan.brainlabel, ROI2))...
    & ~(strcmp(allChan.pacient, 'p1239007 Sko VT61') | strcmp(allChan.pacient, 'p65639 Kam VT63') )); % exclude 2 patients that have only one ROI


%% create a nod table
nod_table = table('Size', [size(ROI_chan,1), 6], 'VariableTypes', {'double','double','double','double','double', 'string'}...
    , 'VariableNames', {'MNI_x','MNI_y','MNI_z','node_color','node_size', 'node_label'});

nod_table.MNI_x = abs(allChan.MNI_x(ROI_chan)); % convert all to the right hemisphere
nod_table.MNI_y = allChan.MNI_y(ROI_chan);
nod_table.MNI_z = allChan.MNI_z(ROI_chan);
% nod_table.node_label  = allChan.name(ROI_chan); % original names of channels
nod_table.node_label  = allChan.neurologyLabel(ROI_chan); % neurologyLabel
% nod_table.node_label  = alphaIncreaseChan.newROI(ROI_chan); % ROI labels

% % assign values to node_color based on node_label
% unique_labels = unique(nod_table.node_label);
% for i = 1:length(unique_labels)
%     nod_table.node_color(strcmp(nod_table.node_label, unique_labels{i})) = i;
% end

% assign values to node_color for each patient
unique_patients = unique(allChan.pacient(ROI_chan), 'stable');
patient_names = allChan.pacient(ROI_chan);
for i = 1:length(unique_patients)
    for j = 1:length(ROI_chan)        
        if strcmp(patient_names{j}, unique_patients{i})
            nod_table.node_color(j) = i;
        end
    end
end

% for i = 1:length(unique_patients)
%     for j = 1:length(ROI_chan)
%         patient_name = regexp(nod_table.node_label(j),' ','split');
%         if strcmp(patient_name{1}{1}, unique_patients{i})
%             nod_table.node_color(j) = i;
%         end
%     end
% end

% % use node_size to distingish patients
% [unique_patients,~,ic] = unique(alphaIncreaseChan.pacient(ROI_chan), 'stable');
% nod_table.node_size = ic/10;

% replace all values in node_size with 1
nod_table.node_size = ones(size(nod_table,1),1);

% Modify the 'node_label' column to delete spaces, (), -, otherwise it causes an error in BrainNet
% nod_table.node_label = strrep(nod_table.node_label, ' ', '_');
% nod_table.node_label = strrep(nod_table.node_label, '(', '');
% nod_table.node_label = strrep(nod_table.node_label, ')', '');
% nod_table.node_label = strrep(nod_table.node_label, '-', ':');
nod_table.node_label = regexprep(nod_table.node_label, '\((\w+).*\)', '$1');


% export the table to a text file
% filename = ['E:\work\PhD\MemoryActions\results\iEEG\connectivity\' ROI1 '-' ROI2 '_channels_alphaIncrease.node'];
% filename = ['E:\work\PhD\MemoryActions\results\iEEG\connectivity\newROIs_channels_alphaIncrease_9pat.node'];
filename = ['F:\Sofia\MemoryActions\results\iEEG\connectivity\visualization of channels\' ROI1 '-' ROI2 '_' num2str(size(nod_table,1)) '_channels_in_' num2str(numel(unique_patients)) ' patients.node'];
writetable(nod_table, filename, 'Delimiter', ' ', 'QuoteStrings', false, 'FileType', 'text', 'WriteVariableNames',0);

%% create a nod table for all channels in 3 ROIs from 9 patients
load('F:\Sofia\MemoryActions\results\iEEG\PAC_9pat_3ROIs_final.mat')

nod_table = table('Size', [size(PAC,2), 6], 'VariableTypes', {'double','double','double','double','double', 'string'}...
    , 'VariableNames', {'MNI_x','MNI_y','MNI_z','node_color','node_size', 'node_label'});

nod_table.MNI_x = [PAC.MNI_x]'; 
nod_table.MNI_y = [PAC.MNI_y]';
nod_table.MNI_z = [PAC.MNI_z]';
nod_table.node_label  = {PAC.neurologyLabel}'; % neurologyLabel

% assign values to node_color based on ROI
unique_labels = unique({PAC.brainlabel}');
for i = 1:length(unique_labels)
    nod_table.node_color(strcmp({PAC.brainlabel}', unique_labels{i})) = i;
end

% replace all values in node_size with 1
nod_table.node_size = ones(size(nod_table,1),1);

% Modify the 'node_label' column to delete spaces, (), -, otherwise it causes an error in BrainNet
nod_table.node_label = regexprep(nod_table.node_label, '\((\w+).*\)', '$1');

filename = ['F:\Sofia\MemoryActions\results\iEEG\connectivity\visualization of channels\all_369Chan_3ROIs_9patients.node'];
writetable(nod_table, filename, 'Delimiter', ' ', 'QuoteStrings', false, 'FileType', 'text', 'WriteVariableNames',0);
