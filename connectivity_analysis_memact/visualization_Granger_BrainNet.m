%%% script creates an edge file with Granger data for BrainNet Viewer for one patient
% with ch pairs showing stat difference between 2 directions during last 2s of delay 
% in edge file, it saves the average significant diff between 2 directions(the net information flow)
% the connections are directional: negative value - direction VTC -> IPL; positive - IPL -> VTC
% the node file is the same as for PLV

%% load file with Granger data
patient_filename = 'E:\work\eeg\motol\pacienti\p2092846 Gru VT68\memact\Granger_permut_stat\Granger_VTC-IPL last 2s delay all_trials_2024-02.mat';
load(patient_filename)

%% create edge file with Granger difference

% significant pairs with significant direction VTC -> IPL
significant_chanPairs_direction1 = find(any(Granger_signif_allPairs_clustcorr < 0,2) &...
                ~any(Granger_signif_allPairs_clustcorr > 0, 2)); % at least one negative value and no positive values
idxChan_direction1 = ROI_chanpairs(significant_chanPairs_direction1, :); % channel indices in each pair            
       
% significant pairs with significant direction IPL -> VTC
significant_chanPairs_direction2 = find(any(Granger_signif_allPairs_clustcorr > 0,2) &...
                ~any(Granger_signif_allPairs_clustcorr < 0, 2)); % at least one positive value and no negative values
idxChan_direction2 = ROI_chanpairs(significant_chanPairs_direction2, :); % channel indices in each pair            

% significant net information flow
Granger_chpairs_signif1 = Granger_signif_allPairs_clustcorr(significant_chanPairs_direction1, :);
Granger_chpairs_signif2 = Granger_signif_allPairs_clustcorr(significant_chanPairs_direction2, :);

% for BrainNet Viewer connectivity data should be in a symmetrical matrix
% Initialize the matrix - all channels x all channels
Granger_diff_matrix = zeros(numel(dataDelay.label), numel(dataDelay.label));
edge_color_matrix = zeros(numel(dataDelay.label), numel(dataDelay.label)); % matrix for color - to set a different color for each direction

% populate the matrix by average significant Granger diff for each signif ch pair
for i = 1:numel(dataDelay.label)
    for j = 1:numel(dataDelay.label)
        isignPair1 = ismember(idxChan_direction1, [i j], 'rows'); % one direction
        if any(isignPair1) > 0
            % find significant freq bins for this ch pair
            ifreq_sign1 = Granger_chpairs_signif1(isignPair1, :)<0; % negative in direction VTC -> IPL
            Granger_diff_matrix(i, j) = mean(Granger_chpairs_signif1(isignPair1, ifreq_sign1)); % average Granger diff        
            Granger_diff_matrix(j, i) = Granger_diff_matrix(i, j); % same for symmetrical position
            edge_color_matrix(i, j) = 1;
            edge_color_matrix(j, i) = edge_color_matrix(i, j);
        end
        isignPair2 = ismember(idxChan_direction2, [i j], 'rows'); % another direction
        if any(isignPair2) > 0
            % find significant freq bins for this ch pair
            ifreq_sign2 = Granger_chpairs_signif2(isignPair2, :)>0; % positive in direction IPL -> VTC
            Granger_diff_matrix(j, i) = mean(Granger_chpairs_signif2(isignPair2, ifreq_sign2)); % average Granger diff 
            Granger_diff_matrix(i, j) = Granger_diff_matrix(j, i); % same for symmetrical position
            edge_color_matrix(i, j) = 2;
            edge_color_matrix(j, i) = edge_color_matrix(i, j);
        end
    end
end

% export it to ASCII text file
% Set the file name and path
[patient_filepath, name] = fileparts(patient_filename);
filename = ['BrainNet_figures\' name '.edge'];
full_path = fullfile(fileparts(patient_filepath), filename);
% Write the matrix to the ASCII text file
dlmwrite(full_path, Granger_diff_matrix, 'delimiter', '\t');

% save also color matrix
filename2 = ['BrainNet_figures\' name '_edge_color.txt'];
full_path = fullfile(fileparts(patient_filepath), filename2);
dlmwrite(full_path, edge_color_matrix, 'delimiter', '\t')
