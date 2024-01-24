% extract data for BrainNet viewer, only one condition

% first look at PN.plotISPC.PlotISPCbins.max_bins to see for what window we want to extract ispc
% data and find indexes of this window - itime and ifreq 
itime = 9;
ifreq = 1;
icondition = 1;
condition = ' del_same';
% time_interval = ' 0-0.5 sec';
time_interval = ' 3.2-3.6 sec';
freq_range = ' 4-6Hz';

% then find indexes of channel pairs with significant ispc in this window; del_same condition
window1_chpairs = PN.plotISPC.PlotISPCbins.bins_idx_chanPairs{icondition,itime,ifreq}; % {kats x time x freq} [ch1 ch2]

% then for these chan pairs and this window extract binned ispc data (average or max per bin)
binned_data = zeros(size(window1_chpairs,1),1);
for ipair = 1:length(binned_data)
    binned_data(ipair) = squeeze(PN.plotISPC.PlotISPCbins.ispc_cats_bins{icondition}(window1_chpairs(ipair,1), window1_chpairs(ipair,2), itime, ifreq)); % one ispc value per each chan pair
end

% for BrainNet Viewer connectivity data should be in a symmetrical matrix
% Initialize the matrix
n = max(window1_chpairs(:));
plv_matrix = zeros(n, n);

% Convert subscripts to linear indices
index1 = sub2ind([n, n], window1_chpairs(:, 1), window1_chpairs(:, 2));
index2 = sub2ind([n, n], window1_chpairs(:, 2), window1_chpairs(:, 1));

% Fill the matrix
plv_matrix(index1) = binned_data;
plv_matrix(index2) = binned_data;

% export it to ASCII text file
% Set the file name and path
filename = ['mean_plv' condition freq_range time_interval '.edge'];
filepath = 'F:\Sofia\MemoryActions\for data analysis\connectivity';
full_path = fullfile(filepath, filename);
% Write the matrix to the ASCII text file
dlmwrite(full_path, plv_matrix, 'delimiter', '\t');

%% create a nod file for BrainNet Viewer with coordinates of all channels
% The node file is defined as an ASCII text file with the suffix ‘node’. In the node file, there are 6 columns: 
% columns 1-3 represent node coordinates, column 4 represents node colors, column 5 represents node sizes, and the last column represents node labels. 
% Please note, a symbol ‘-‘ (no ‘’) in column 6 means no labels and blank characters in label would cause an error

nod_table = table('Size', [n, 6], 'VariableTypes', {'double','double','double','double','double', 'string'}...
    , 'VariableNames', {'MNI_x','MNI_y','MNI_z','node_color','node_size', 'node_label'});

nod_table.MNI_x = [PN.E.CH.H.channels(PN.plotISPC.channels(1:n)).MNI_x]';
nod_table.MNI_y = [PN.E.CH.H.channels(PN.plotISPC.channels(1:n)).MNI_y]';
nod_table.MNI_z = [PN.E.CH.H.channels(PN.plotISPC.channels(1:n)).MNI_z]';
nod_table.node_label  = {PN.E.CH.brainlabels(PN.plotISPC.channels(1:n)).label}'; % all brainlabels = ROIs

% assign values to node_color based on node_label
unique_labels = unique(nod_table.node_label);
for i = 1:length(unique_labels)
    nod_table.node_color(strcmp(nod_table.node_label, unique_labels{i})) = i;
end

% replace all values in node_size with 1
nod_table.node_size = ones(size(nod_table,1),1);

% export the table to a text file
filename2 = ['PT_channels4_' condition time_interval '_VT62.node'];
full_path2 = fullfile(filepath, filename2);
writetable(nod_table, full_path2, 'Delimiter', ' ', 'QuoteStrings', false, 'FileType', 'text', 'WriteVariableNames',0);

%% 20.11.2023
%%% extract PLV difference between 2 conditions for BrainNet viewer from PN.plotISPC.ispc_cats_clust_corr{3}
%%% takes 1 max plv value for each chan pair

time_interval = PN.plotISPC.stat_cats_time(1:2);
freq_range = PN.plotISPC.stat_cats_freq;
plv_diff = PN.plotISPC.ispc_cats_clust_corr{3}; % ch1 x ch2 x time x freq
plv_diff_matrix = NaN(size(plv_diff,1), size(plv_diff,2));

% for each ch pair leave only abs max value
for ichn1 = 1:size(plv_diff,1)
    for ichn2 = ichn1:size(plv_diff,2) % iterate only through the upper triangular part of the matrix
        if ~isnan(plv_diff(ichn1,ichn2, 1, 1))
            plv_diff_matrix(ichn1,ichn2) = max(max(abs(squeeze(plv_diff(ichn1,ichn2, :, :)))));
            plv_diff_matrix(ichn2,ichn1) = plv_diff_matrix(ichn1,ichn2); % symmetrical position
        end
    end
end

% Find rows and columns with NaN values
nanRows = all(isnan(plv_diff_matrix), 2);
nanCols = all(isnan(plv_diff_matrix), 1);

% Use logical indexing to keep only rows and columns without NaN
plv = plv_diff_matrix(~nanRows, ~nanCols);

% replace remaining nan by zeros, otherwise it'll produce an error in brainnet
plv(isnan(plv)) = 0;

% export it to ASCII text file
% Set the file name and path

filename = ['VT63_max_plv_diff_between_conditions_' strrep(num2str(freq_range),'  ','-') ' Hz_' strrep(num2str(time_interval),'         ','-') '.edge'];
filepath = 'F:\Sofia\MemoryActions\for data analysis\connectivity';
full_path = fullfile(filepath, filename);
% Write the matrix to the ASCII text file
dlmwrite(full_path, plv, 'delimiter', '\t');

%% create a nod file for BrainNet Viewer with coordinates of all channels
% The node file is defined as an ASCII text file with the suffix ‘node’. In the node file, there are 6 columns: 
% columns 1-3 represent node coordinates, column 4 represents node colors, column 5 represents node sizes, and the last column represents node labels. 
% Please note, a symbol ‘-‘ (no ‘’) in column 6 means no labels and blank characters in label would cause an error

nod_table = table('Size', [size(plv,1), 6], 'VariableTypes', {'double','double','double','double','double', 'string'}...
    , 'VariableNames', {'MNI_x','MNI_y','MNI_z','node_color','node_size', 'node_label'});

nod_table.MNI_x = [PN.E.CH.H.channels(PN.plotISPC.channels(~nanRows)).MNI_x]';
nod_table.MNI_y = [PN.E.CH.H.channels(PN.plotISPC.channels(~nanRows)).MNI_y]';
nod_table.MNI_z = [PN.E.CH.H.channels(PN.plotISPC.channels(~nanRows)).MNI_z]';
nod_table.node_label  = {PN.E.CH.brainlabels(PN.plotISPC.channels(~nanRows)).label}'; % all brainlabels = ROIs

% assign values to node_color based on node_label
unique_labels = unique(nod_table.node_label);
for i = 1:length(unique_labels)
    nod_table.node_color(strcmp(nod_table.node_label, unique_labels{i})) = i;
end

% replace all values in node_size with 1
nod_table.node_size = ones(size(nod_table,1),1);

% export the table to a text file
filename2 = ['VT63_PT_channels_plv_diff_' strrep(num2str(freq_range),'  ','-') ' Hz_' strrep(num2str(time_interval),'         ','-') '.node'];
full_path2 = fullfile(filepath, filename2);
writetable(nod_table, full_path2, 'Delimiter', ' ', 'QuoteStrings', false, 'FileType', 'text', 'WriteVariableNames',0);
