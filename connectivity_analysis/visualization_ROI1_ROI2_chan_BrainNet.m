%% script creates a nod file for BrainNet Viewer with coordinates of all channels in LOC and IPL (or any other ROI)
% The node file is defined as an ASCII text file with the suffix node. In the node file, there are 6 columns: 
% columns 1-3 represent node coordinates, column 4 represents node colors, column 5 represents node sizes, and the last column represents node labels. 
% Please note, a symbol - (no ) in column 6 means no labels and blank characters in label would cause an error

% to be configured for each pair of ROIs you want to visualize
ROI1 = 'LOC';
ROI2 = 'pVTC';

% table with all channels from all patients showing an alpha increase during the delay
alphaIncrease = 'E:\work\PhD\MemoryActions\results\iEEG\delayed condition\January2024_9pat\Response2XLS_CM_Memact_CHilbert_8-13Hz_AnyResp_OTP_-0.5-5.9_refBipo_Ep2024-01_encod+del_CHMult_2024-01-25_13-01-46.xls';
alphaIncreaseChan = readtable(alphaIncrease, 'ReadRowNames',true);

% find all channels with an alpha increase in ROI1 and ROI2 for patients that have both of them
% ROI_chan = find((strcmp(alphaIncreaseChan.mybrainlabel, ROI1) | strcmp(alphaIncreaseChan.mybrainlabel, ROI2))...
%     & ~(strcmp(alphaIncreaseChan.pacient, 'p65639') | strcmp(alphaIncreaseChan.pacient, 'p1883612_VT66') )); % exclude 2 patients that have only one ROI
ROI_chan = find(strcmp(alphaIncreaseChan.newROI, 'AnG') | strcmp(alphaIncreaseChan.newROI, 'SMG')...
    | strcmp(alphaIncreaseChan.newROI, 'ITG') | strcmp(alphaIncreaseChan.newROI, 'VTC') | strcmp(alphaIncreaseChan.newROI, 'MTL')...
    | strcmp(alphaIncreaseChan.newROI, 'precun'));

% create a nod table
nod_table = table('Size', [size(ROI_chan,1), 6], 'VariableTypes', {'double','double','double','double','double', 'string'}...
    , 'VariableNames', {'MNI_x','MNI_y','MNI_z','node_color','node_size', 'node_label'});

nod_table.MNI_x = abs(alphaIncreaseChan.MNI_x(ROI_chan)); % convert all to the right hemisphere
nod_table.MNI_y = alphaIncreaseChan.MNI_y(ROI_chan);
nod_table.MNI_z = alphaIncreaseChan.MNI_z(ROI_chan);
% nod_table.node_label  = alphaIncreaseChan.name(ROI_chan); % original names of channels
nod_table.node_label  = alphaIncreaseChan.newROI(ROI_chan); % ROI labels

% assign values to node_color based on node_label
unique_labels = unique(nod_table.node_label);
for i = 1:length(unique_labels)
    nod_table.node_color(strcmp(nod_table.node_label, unique_labels{i})) = i;
end

% % assign values to node_color for each patient
% unique_patients = unique(alphaIncreaseChan.pacient(ROI_chan), 'stable');
% 
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

% % Modify the 'node_label' column to delete spaces, (), -, otherwise it causes an error in BrainNet
% nod_table.node_label = strrep(nod_table.node_label, ' ', '_');
% nod_table.node_label = strrep(nod_table.node_label, '(', '');
% nod_table.node_label = strrep(nod_table.node_label, ')', '');
% nod_table.node_label = strrep(nod_table.node_label, '-', ':');

% export the table to a text file
% filename = ['E:\work\PhD\MemoryActions\results\iEEG\connectivity\' ROI1 '-' ROI2 '_channels_alphaIncrease.node'];
filename = ['E:\work\PhD\MemoryActions\results\iEEG\connectivity\newROIs_channels_alphaIncrease_9pat.node'];

writetable(nod_table, filename, 'Delimiter', ' ', 'QuoteStrings', false, 'FileType', 'text', 'WriteVariableNames',0);

