% script for defining ROIs in PAC structure in memact test
% CB=CBrainPlot();
% PAC = CB.StructFindLoad ('d:\eeg\motol\iEEG_scripts\logs\StructFind PAC_memact_()_(TOTZ-LgG-LPHT-PHG-FG-FFG-FuG-ITG-MTG-STG-TP-Hi-Ent-PCun-precun-SMG-AnG-IPL).xlsx',1);
% PAC = CB.StructFindLoad('d:\eeg\motol\iEEG_scripts\logs\StructFind PAC_memact_()_().xlsx',1);

for i = 1:length(PAC)
    % Extract MNI coordinates for the current channel
    MNI_x = PAC(i).MNI_x;
    MNI_y = PAC(i).MNI_y;
    MNI_z = PAC(i).MNI_z;
    
    % name ROIs(brainlabels) based on neurologyLabel and MNI coordinates
%     if contains(PAC(i).neurologyLabel, 'Ang', 'IgnoreCase',true) % parietal ROIs from 1.11.2023
%         PAC(i).brainlabel = 'AnG';
%     elseif contains(PAC(i).neurologyLabel, 'SMG', 'IgnoreCase',true) 
%         PAC(i).brainlabel = 'SMG';
    if contains(PAC(i).neurologyLabel, 'Ang', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'SMG', 'IgnoreCase',true) % parietal ROIs, from 11.12.2023
        PAC(i).brainlabel = 'IPL';
    elseif contains(PAC(i).neurologyLabel, 'PCun', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel,'precun', 'IgnoreCase',true)
        PAC(i).brainlabel = 'precun';
    elseif contains(PAC(i).neurologyLabel, 'SPL', 'IgnoreCase',true)
        PAC(i).brainlabel = 'SPL';
    elseif (abs(MNI_x) >= 29 && abs(MNI_x) <= 59) && (MNI_y >= -85 && MNI_y <= -45) && (MNI_z >= -20 && MNI_z <= 8) % 13.12.2023
%     elseif (abs(MNI_x) >= 25 && abs(MNI_x) <= 68) && (MNI_y >= -90 && MNI_y <= -37) && (MNI_z >= -25 && MNI_z <= 10)
        PAC(i).brainlabel = 'LOC';  % temporal ROIs
    elseif contains(PAC(i).neurologyLabel, 'Hi', 'IgnoreCase',true) 
        PAC(i).brainlabel = 'Hip';
    elseif (contains(PAC(i).neurologyLabel, 'FuG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'PHG', 'IgnoreCase',true)...
            || contains(PAC(i).neurologyLabel, 'LgG', 'IgnoreCase',true)) && MNI_y > -35 
        PAC(i).brainlabel = 'aVTC';
    elseif (contains(PAC(i).neurologyLabel, 'MTG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'STG', 'IgnoreCase',true)...
            || contains(PAC(i).neurologyLabel, 'TP', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'ITG', 'IgnoreCase',true)) && MNI_y > -20 
        PAC(i).brainlabel = 'aLTC';
    elseif (contains(PAC(i).neurologyLabel, 'MTG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'STG', 'IgnoreCase',true)...
            || contains(PAC(i).neurologyLabel, 'ITG', 'IgnoreCase',true)) && MNI_y <= -20 
        PAC(i).brainlabel = 'pLTC';
    elseif (contains(PAC(i).neurologyLabel, 'FuG', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'PHG', 'IgnoreCase',true) ...
            || contains(PAC(i).neurologyLabel, 'LgG', 'IgnoreCase',true)) && MNI_y <= -35  % 13.12.2023
        PAC(i).brainlabel = 'pVTC';    
    end
    
%     if contains(PAC(i).neurologyLabel, 'Ang', 'IgnoreCase',true) || contains(PAC(i).neurologyLabel, 'SMG', 'IgnoreCase',true) % parietal ROIs, from 6.11.2023
%         PAC(i).brainlabel = 'IPL';
%     elseif (abs(MNI_x) >= 25 && abs(MNI_x) <= 68) && (MNI_y >= -90 && MNI_y <= -37) && (MNI_z >= -25 && MNI_z <= 10)
%         PAC(i).brainlabel = 'LOC';  % temporal ROIs
%     end
end

% % Convert the structure to a table
% PAC_table = struct2table(PAC);
% 
% % Write the table to an Excel file
% writetable(PAC_table, 'd:\eeg\motol\iEEG_scripts\logs\StructFind PAC_memact_ROI_7pat.xlsx', 'Sheet', 1, 'Range', 'A1');
