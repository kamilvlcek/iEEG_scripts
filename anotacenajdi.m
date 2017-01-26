function []=anotacenajdi(adresar,filter)
% ANOTACENAJDI vypise anotace ze vsech souboru mat 
% napriklad anotacenajdi('D:\eeg\motol\pacienti\Pluhackova\');
% filter je retezec, ktery musi obsahovat jmeno souboru

if ~exist('filter','var'), filter = []; end %pokud neni definovan filter, vytvorim jako prazdny
files = dir(fullfile(adresar, '*.mat'));

for f = 1:numel(files)
    filename = [ adresar files(f).name];
    if isempty(filter) || ~isempty(strfind(filename,filter))
        %zpracovavam jen soubory vyhovujici filtru, pokud nejaky existuje
        vars = whos('-file',filename);
        if ismember('header', {vars.name})
            load(filename,'header','tabs');
            disp(['*** ' files(f).name]);
            anotace(header,tabs);    
            clear header d tabs fs;
        elseif ismember('evts', {vars.name})
            load(filename,'evts','tabs');
            disp(['*** ' files(f).name]);
            anotace([],tabs,evts);
            clear header d tabs fs;
        else
           disp(['*** ' files(f).name '- no header']);
        end
    end
end