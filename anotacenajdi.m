function []=anotacenajdi(adresar)
% ANOTACENAJDI vypise anotace ze vsech souboru mat 
% napriklad anotacenajdi('D:\eeg\motol\pacienti\Pluhackova\');

files = dir(fullfile(adresar, '*.mat'));

for f = 1:numel(files)
    load([ adresar files(f).name],'header','tabs');
    disp(['*** ' files(f).name]);
    anotace(header,tabs);
    clear header d tabs fs;
end