function [ dirs ] = adresar_najdi( pacientid )
%ADRESAR vrati adresar z kodu pacienta
%   hleda podadresare, ktere vyhovuji v eeg\motol\pacienti na D a W

dirs = {};
dirindex = 1;
basedirs = {'d:\eeg\motol\pacienti\','w:\eeg\motol\pacienti\'};
for d = 1:length(basedirs)
    files = dir(fullfile(basedirs{d}));
    for f = 1:numel(files)
        dirname = [ basedirs{d} files(f).name '\'];
        if ~isempty(strfind(dirname,pacientid)) && isdir(dirname)            
            dirs{dirindex} = dirname; %#ok<AGROW>            
            dirindex = dirindex +1;
        end
    end
end

dirs = dirs{1};

end

