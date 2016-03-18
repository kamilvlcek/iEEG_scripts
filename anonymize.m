function []=anonymize(adresar,id)
% ANONYMIZE prepise header.pacientID hodnotou id
% napriklad anonymize('D:\eeg\motol\pacienti\VT8\','p79');
% vse ulozi znova pod jmenem s pismenem X na konci


files = dir(fullfile(adresar, '*.mat'));

for f = 1:numel(files)
    disp(['nacitam ' files(f).name]);
    load([ adresar files(f).name]);
    newfilename = strrep(files(f).name, '.mat', '_X.mat');
    if exist('header','var')
        if isfield(header,'patientID')
            disp(['puvodni id: ' header.patientID]);
        else
            disp('prazdne id ');
        end
        header.patientID = id;
        header = rmfield(header, 'recordID'); %pokud tam je
        disp(['ukladam ' newfilename]);
        save([ adresar newfilename], '-regexp', '^(?!(adresar|id|files|f|newfilename)$).','-v7.3');    
        disp(['ulozeno s header: ' newfilename]);
    elseif exist('H','var')
        if isfield(H,'patientID')
            disp(['puvodni id: ' H.patientID]);
        else
            disp('prazdne id ');
        end
        H.patientID = id;
        disp(['ukladam ' newfilename]);
        save([ adresar newfilename], '-regexp', '^(?!(adresar|id|files|f|newfilename)$).','-v7.3')
        disp(['ulozeno s H: ' newfilename]);
    else
        disp(['zadny header, neulozeno: ' files(f).name]);
    end
    clearvars -except adresar id files f; %smazu vsechny promenne
end
disp('hotovo');