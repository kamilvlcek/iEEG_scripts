function [pacienti] = pacientFolderSelect(pacienti,setup)
%PACIENTFOLDERSELECT selects the existing folder for each pacient
%   if pacienti().folder is cell with multiple possible folders, it selects the correct one with test-name subfolder 
%   it also sets the field 'found' indicating if the pacient folder and the test-name (e.g. aedist) subfolder exist
  
for p = 1:numel(pacienti)
    if iscell(pacienti(p).folder)
        found = false;
        for ifolder = 1:numel(pacienti(p).folder)
            ff = fullfile(setup.basedir, pacienti(p).folder{ifolder},setup.subfolder);
            if exist(ff,'dir')==7 %the folder exists
                pacienti(p).folder = pacienti(p).folder{ifolder};
                found = true;
                break; %we do not want to go to next ifolder
            end
        end
        if ~found %we have to select one folder from the cell array
            pacienti(p).folder = pacienti(p).folder{1};
            pacienti(p).found = false;
        else
            pacienti(p).found = true;
        end
    else
        ff = fullfile(setup.basedir, pacienti(p).folder,setup.subfolder);
        pacienti(p).found = exist(ff,'dir')==7;        
    end
    
end
end

