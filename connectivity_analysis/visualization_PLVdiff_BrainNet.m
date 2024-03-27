%% script creates a nod and edge file for BrainNet Viewer for one patient
% with ch pairs showing stat difference between last 2 s of delay and bs
% in edge file, it saves the average significant diff between PLV delay and PLV bs

% % load file with PLV data
% % patient_filename = 'E:\work\eeg\motol\pacienti\p2092846 Gru VT68\memact\PLV_permut_stat\PLV_VTC-IPL_last 1.5s delay_vs_bs_all_trials_2024-02.mat';
% patient_filename = 'E:\work\eeg\motol\pacienti\p2092846 Gru VT68\memact\PLV_permut_stat\PLV_VTC-IPL_last 2s delay_vs_bs_all_trials_2024-02.mat';
% load(patient_filename)

% PLV_data =  'PLV_Hip-IPL_last 2s delay_vs_bs_all_trials_200permut_2024-03.mat';
PLV_data =  'PLV_VTC-IPL_last 2s delay_vs_bs_all_trials_200permut_2024-03.mat';
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact, some are marked not to analyze as they don't have ROIs we want

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if isfile([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data])
            % load the PLV data with statistics for this patient
            load([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data]);
            
            %% first create edge file with plv difference
            % % find indexes of chan pairs with significant PLV difference, only positive delay > bs (any positive value, even one)
            % significant_chanPairs = sum(plv_signif_allPairs_clustcorr,2) > 0;
            % ROI_chanpairs_signif = ROI_chanpairs(significant_chanPairs, :);
            
            % % Find rows (ch pairs) with at least 2 positive values (2 freq bins)
            % positive_values_count = sum(plv_signif_allPairs_clustcorr > 0, 2);
            % significant_chanPairs = find(positive_values_count >= 2);
            
            % Find rows (ch pairs) with at least 2 consecutive positive values
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
            
            ROI_chanpairs_signif = ROI_chanpairs(significant_chanPairs, :);
            
            % significant PLV diff
            PLVdiff_chpairs_signif = plv_signif_allPairs_clustcorr(significant_chanPairs, :);
            
            % for BrainNet Viewer connectivity data should be in a symmetrical matrix
            % Initialize the matrix - all channels x all channels
            plv_diff_matrix = zeros(numel(PLVCond1.label), numel(PLVCond1.label));
            
            % populate the matrix by average significant plv diff for each signif ch pair
            for i = 1:numel(PLVCond1.label)
                for j = 1:numel(PLVCond1.label)
                    isignPair = ismember(ROI_chanpairs_signif, [i j], 'rows');
                    if any(isignPair) > 0
                        % find significant freq bins for this ch pair
                        ifreq_sign = PLVdiff_chpairs_signif(isignPair, :)>0;
                        plv_diff_matrix(i, j) = mean(PLVdiff_chpairs_signif(isignPair, ifreq_sign)); % average PLV diff
                        plv_diff_matrix(j, i) = plv_diff_matrix(i, j); % save also to the symmetrical position
                        
                    end
                end
            end
            
            % export it to ASCII text file
            % Set the file name and path
%             [patient_filepath, name] = fileparts(patient_filename);
            [~, name] = fileparts(PLV_data);
            patient_filepath = [basedir pacienti(p).folder '\' subfolder '\'];
            filename = ['BrainNet_figures\' name '.edge'];
            full_path = fullfile(patient_filepath, filename);
            % Write the matrix to the ASCII text file
            dlmwrite(full_path, plv_diff_matrix, 'delimiter', '\t');
            
            
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
            filename2 = ['BrainNet_figures\' name '.node'];
            full_path2 = fullfile(patient_filepath, filename2);
            writetable(nod_table, full_path2, 'Delimiter', ' ', 'QuoteStrings', false, 'FileType', 'text', 'WriteVariableNames',0);
        end
    end
end
