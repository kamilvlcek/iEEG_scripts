%% script to create a nod and edge file for BrainNet Viewer for each patient
% with chn pairs showing significant PLV in any period

% load an aggregated table with indices of all significant pairs for all three periods across patients
% filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_VTC-IPL_aggregated_significant_chnPairs.mat';
filenameTable = 'F:\Sofia\MemoryActions\results\iEEG\connectivity\group data\PLV_Hip-IPL_aggregated_significant_chnPairs.mat';
load(filenameTable);

% PLV_data =  'PLV_VTC-IPL_last 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
PLV_data =  'PLV_Hip-IPL_last 1.9s delay_vs_bs_all_trials_200permut_2024-04.mat';
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact

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
    
    if isfile([filePath PLV_data])
        load([filePath PLV_data]);
        
        %% first create edge file with significant connections
        
        % Get indices of aggregated significant channel pairs for the current subject
        significantidxChan = aggregTable(subjIndex).idxChan_in_Pairs; % n signif pairs x 2 (each pair with chan indices)
        
        % for BrainNet Viewer connectivity data should be in a symmetrical matrix
        plv_matrix = zeros(numel(PLVCond1.label), numel(PLVCond1.label)); % all channels x all channels
        
        % populate the matrix by 1 for each signif ch pair
        for i = 1:numel(PLVCond1.label)
            for j = 1:numel(PLVCond1.label)
                isignPair = ismember(significantidxChan, [i j], 'rows');
                if any(isignPair) > 0
                    plv_matrix(i, j) = 1; % indicate the significant chn pair
                    plv_matrix(j, i) = plv_matrix(i, j); % save also to the symmetrical position
                end
            end
        end
        
        % export it to ASCII text file
        [~, name] = fileparts(PLV_data);
        % extract roi names fron the name
        pattern = 'PLV_[A-Za-z]+-[A-Za-z]+';
        part_name = regexp(name, pattern, 'match', 'once');
        patient_filepath = [basedir pacienti(p).folder '\' subfolder '\'];
        filename = ['BrainNet_figures\' part_name '.edge'];
        full_path = fullfile(patient_filepath, filename);
        % Write the matrix to the ASCII text file
        dlmwrite(full_path, plv_matrix, 'delimiter', '\t');
        
        
        %% create a nod file with coordinates of all channels
        % The node file is defined as an ASCII text file with the suffix ‘node’. In the node file, there are 6 columns:
        % columns 1-3 represent node coordinates, column 4 represents node colors, column 5 represents node sizes, and the last column represents node labels.
        % Please note, a symbol ‘-‘ (no ‘’) in column 6 means no labels and blank characters in label would cause an error
        
        nod_table = table('Size', [numel(PLVCond1.label), 6], 'VariableTypes', {'double','double','double','double','double', 'string'}...
            , 'VariableNames', {'MNI_x','MNI_y','MNI_z','node_color','node_size', 'node_label'});
        
        nod_table.MNI_x = [dataCond1.channelInfo.MNI_x]';
        nod_table.MNI_y = [dataCond1.channelInfo.MNI_y]';
        nod_table.MNI_z = [dataCond1.channelInfo.MNI_z]';
        nod_table.node_label  = dataCond1.label; % original bipolar names of the channels
        
        % assign values to node_color based on ROI
        unique_ROI = unique({dataCond1.channelInfo.ROI}, 'stable');
        for i = 1:size(nod_table,1)
            for j = 1: numel(unique_ROI)
                if dataCond1.channelInfo(i).ROI == unique_ROI{j}
                    nod_table.node_color(i) = j;
                end
            end
        end
        
        % replace all values in node_size with 1
        nod_table.node_size = ones(size(nod_table,1),1);
        
        % export the table to a text file
        filename2 = ['BrainNet_figures\' part_name '.node'];
        full_path2 = fullfile(patient_filepath, filename2);
        writetable(nod_table, full_path2, 'Delimiter', ' ', 'QuoteStrings', false, 'FileType', 'text', 'WriteVariableNames',0);
    end
end

