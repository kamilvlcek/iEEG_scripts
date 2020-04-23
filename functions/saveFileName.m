function [filename] = saveFileName(oldName, suffix)
%SAVEFILENAME Normalize the save file's filename by removing old suffixes and adding a new one
%TODO: Enclose the function in an utility class
    filename = oldName;
    % Remove the old save suffixes
    filename=strrep(filename,'_CHilb','');
    filename=strrep(filename,'_CiEEG','');
    filename=strrep(filename,'_CHMult','');
    [pathstr, fname, ext] = CiEEGData.matextension(filename);         
    filename = fullfile(pathstr, [fname '_' suffix ext]);
end