% script for defining ROIs in PAC structure in memact test
% CB=CBrainPlot();
% PAC = CB.StructFindLoad ('d:\eeg\motol\iEEG_scripts\logs\StructFind PAC_memact_()_(TOTZ-LgG-LPHT-PHG-FG-FFG-FuG-ITG-MTG-STG-TP-Hi-Ent-PCun-precun-SMG-AnG-IPL).xlsx',1);
% PAC = CB.StructFindLoad('d:\eeg\motol\iEEG_scripts\logs\StructFind PAC_memact_()_().xlsx',1);

%9.02.2024
% table with all channels from all patients showing an alpha increase during the delay
% table = 'E:\work\PhD\MemoryActions\results\iEEG\delayed condition\January2024_9pat\Response2XLS_CM_Memact_CHilbert_8-13Hz_AnyResp_OTP_-0.5-5.9_refBipo_Ep2024-01_encod+del_CHMult_2024-01-25_13-01-46.xls';

% table with all channels from 9 patients
table = 'F:\Sofia\MemoryActions\raw data\electrodes localization\StructFind PAC_memact_all_chan_9pat.xlsx';
% tableChan = readtable(table, 'ReadRowNames',true);
tableChan = readtable(table);
PAC = table2struct(tableChan); % I want to change ROIs in this table, they will be in a field newROI
%%

for i = 1:length(PAC)
    % Extract MNI coordinates for the current channel
    MNI_x = PAC(i).MNI_x;
    MNI_y = PAC(i).MNI_y;
    MNI_z = PAC(i).MNI_z;
    
%     % name ROIs(brainlabels) based on neurologyLabel and MNI coordinates % 13.12.2023
%     if contains(PAC(i).neurologyLabel, 'Ang', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'SMG', 'IgnoreCase',true) % parietal ROIs, from 11.12.2023
%         PAC(i).brainlabel = 'IPL';
%     elseif contains(PAC(i).neurologyLabel, 'PCun', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel,'precun', 'IgnoreCase',true)
%         PAC(i).brainlabel = 'precun';
%     elseif contains(PAC(i).neurologyLabel, 'SPL', 'IgnoreCase',true)
%         PAC(i).brainlabel = 'SPL';
%     elseif (abs(MNI_x) >= 29 && abs(MNI_x) <= 59) && (MNI_y >= -85 && MNI_y <= -45) && (MNI_z >= -20 && MNI_z <= 8) % 13.12.2023
% %     elseif (abs(MNI_x) >= 25 && abs(MNI_x) <= 68) && (MNI_y >= -90 && MNI_y <= -37) && (MNI_z >= -25 && MNI_z <= 10)
%         PAC(i).brainlabel = 'LOC';  % temporal ROIs
%     elseif contains(PAC(i).neurologyLabel, 'Hi', 'IgnoreCase',true) 
%         PAC(i).brainlabel = 'Hip';
%     elseif (contains(PAC(i).neurologyLabel, 'FuG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'PHG', 'IgnoreCase',true)...
%             || contains(PAC(i).neurologyLabel, 'LgG', 'IgnoreCase',true)) && MNI_y > -35 
%         PAC(i).brainlabel = 'aVTC';
%     elseif (contains(PAC(i).neurologyLabel, 'MTG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'STG', 'IgnoreCase',true)...
%             || contains(PAC(i).neurologyLabel, 'TP', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'ITG', 'IgnoreCase',true)) && MNI_y > -20 
%         PAC(i).brainlabel = 'aLTC';
%     elseif (contains(PAC(i).neurologyLabel, 'MTG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'STG', 'IgnoreCase',true)...
%             || contains(PAC(i).neurologyLabel, 'ITG', 'IgnoreCase',true)) && MNI_y <= -20 
%         PAC(i).brainlabel = 'pLTC';
%     elseif (contains(PAC(i).neurologyLabel, 'FuG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'PHG', 'IgnoreCase',true) ...
%             || contains(PAC(i).neurologyLabel, 'LgG', 'IgnoreCase',true)) && MNI_y <= -35  % 13.12.2023
%         PAC(i).brainlabel = 'pVTC';    
%     end
    
    % new ROIs(brainlabels) based on neurologyLabel only (anatomical, not functional) % 12.02.2024
    if (contains(PAC(i).neurologyLabel, 'Ang', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'SMG', 'IgnoreCase',true))...
            && (~contains(PAC(i).neurologyLabel, 'MOG', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'POP', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'INS', 'IgnoreCase',true)...
            && ~contains(PAC(i).neurologyLabel, 'STG', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'MTG', 'IgnoreCase',true) ...
            && ~contains(PAC(i).neurologyLabel, 'SPL', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'Het', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'bone', 'IgnoreCase',true))
        PAC(i).brainlabel = 'IPL';
    elseif (contains(PAC(i).neurologyLabel, 'PCun', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel,'precun', 'IgnoreCase',true)) && ~contains(PAC(i).neurologyLabel, 'Het', 'IgnoreCase',true)
        PAC(i).brainlabel = 'precun';
    elseif contains(PAC(i).neurologyLabel, 'SPL', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'Het', 'IgnoreCase',true)
        PAC(i).brainlabel = 'SPL';
%     elseif (contains(PAC(i).neurologyLabel, 'Hi', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'PHG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'Ent', 'IgnoreCase',true))
%         PAC(i).brainlabel = 'MTL';
    elseif contains(PAC(i).neurologyLabel, 'Hi', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'Het', 'IgnoreCase',true)
        PAC(i).brainlabel = 'Hip';
    elseif (contains(PAC(i).neurologyLabel, 'ITG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'FuG', 'IgnoreCase',true)...
            || contains(PAC(i).neurologyLabel, 'LgG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'MOG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'PHG', 'IgnoreCase',true))...
            && (MNI_y < 0 && MNI_y > -82) && (MNI_z < 15)...
            && (~contains(PAC(i).neurologyLabel, 'Amg', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'Ent', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'Het', 'IgnoreCase',true)...
            && ~contains(PAC(i).neurologyLabel, 'V1', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'AnG', 'IgnoreCase',true) && ~contains(PAC(i).neurologyLabel, 'bone', 'IgnoreCase',true))
        
        PAC(i).brainlabel = 'VTC';
    end

end

% % Convert the structure to a table
PAC_table = struct2table(PAC);
% 
% % Write the table to an Excel file
% writetable(PAC_table, 'd:\eeg\motol\iEEG_scripts\logs\StructFind PAC_memact_ROI_7pat.xlsx', 'Sheet', 1, 'Range', 'A1');
writetable(PAC_table, 'F:\Sofia\MemoryActions\results\iEEG\connectivity\StructFind PAC_memact_all_chan_9pat_ROI2.xlsx', 'Sheet', 1, 'Range', 'A1');


