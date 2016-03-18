function []=anonymize(adresar,id)
% ANONYMIZE prepise header.pacientID hodnotou id
% napriklad anonymize('D:\eeg\motol\pacienti\VT8\','p79');
% vse ulozi znova pod jmenem s pismenem X na konci


files = dir(fullfile(adresar, '*.mat'));

for f = 1:numel(files)
    disp(['*** nacitam ' files(f).name ' ***']);
    load([ adresar files(f).name]);
    newfilename = strrep(files(f).name, '.mat', '_X.mat');
    ulozit = false;
    if exist('header','var')
        if isfield(header,'patientID')
            if strcmp(header.patientID,id)
                disp(['puvodni id totozne s novym: ' header.patientID]);
            else
                disp(['puvodni id jine nez nove: ' header.patientID]);
                header.patientID = id;
                ulozit = true;
            end
        else
            disp('prazdne id ');
            header.patientID = id;
            ulozit = true;
        end        
        if isfield(header,'recordID')
            header = rmfield(header, 'recordID'); %pokud tam je
            disp('recordID smazano ');
            ulozit = true;
        end
        if ulozit
            disp(['ukladam ' newfilename]);
            save([ adresar newfilename], '-regexp', '^(?!(adresar|id|files|f|newfilename)$).','-v7.3');    
            disp(['ulozeno s header: ' newfilename]);
        else
            disp(['zadna zmena, neulozeno: ' newfilename]);
        end        
    elseif exist('H','var')
        if isfield(H,'patientID')
            if strcmp(H.patientID,id)
                disp(['puvodni id totozne s novym: ' H.patientID]);
            else
                disp(['puvodni id jine nez nove: ' H.patientID]);
                H.patientID = id;
                ulozit = true;
            end
        else
            disp('prazdne id ');
            H.patientID = id;
            ulozit = true;
        end
        H.patientID = id;
        if isfield(H,'recordID')
            H = rmfield(H, 'recordID'); %pokud tam je
            disp('recordID smazano ');
            ulozit = true;
        end
        if ulozit
            disp(['ukladam ' newfilename]);
            save([ adresar newfilename], '-regexp', '^(?!(adresar|id|files|f|newfilename)$).','-v7.3')
            disp(['ulozeno s H: ' newfilename]);
        else
            disp(['zadna zmena, neulozeno: ' newfilename]);
        end
    else
        disp(['zadny header, neulozeno: ' files(f).name]);
    end
    clearvars -except adresar id files f; %smazu vsechny promenne
end
disp('hotovo');