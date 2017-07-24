function [ E ] = pacient_load( nick,test,filename )
%PACIENT_LOAD nacte soubor pacienta, bud ze zdrojvych data (pokud neni udano filename) nebo uz ulozeny soubor
%   zdrojova data rovnou zpracuje, rozepochuje atd 
if strcmp(test,'aedist')
    pacienti = pacienti_aedist(); %nactu celou strukturu pacientu
    subfolder = 'Aedist';
end
nalezen = false;
for p = 1:numel(pacienti)
    if strfind(pacienti(p).folder,nick)
        nalezen = true;
        break; %nasel jsem pacienta
    end        
end
if ~nalezen
    error(['pacient nenalezen: ' nick]);    
end
disp(['loading pacient ' pacienti(p).folder ]);
basedir = 'D:\eeg\motol\pacienti\';

if exist('filename','var')
    fullfilename = [basedir pacienti(p).folder '\' subfolder '\' filename];
    if exist(fullfilename,'file')==2
        if strfind(filename,'CHilbert') 
            E = CHilbert(fullfilename);
        else
            E = CiEEGData(fullfilename);
        end      
    else
        E = [];
        disp(['soubor neexistuje: ' fullfilename]);
    end
else
    %EEG data
    load([ basedir pacienti(p).folder '\' pacienti(p).data]);
    if ~exist('mults','var'), mults = [];  end
    E = CiEEGData(d,tabs,fs,mults);

    %header
    load([ basedir pacienti(p).folder '\' pacienti(p).header]);
    E.GetHHeader(H);

    %epievents
    filename = [ basedir pacienti(p).folder '\' subfolder '\' pacienti(p).epievents];
    if exist(filename,'file')~=2
        disp(['no epievents:' pacienti(p).epievents]);
        return; 
    end
    load(filename);
    E.GetEpiEvents(DE);

    %vyradim kanaly
    if ~isfield(pacienti(p),'rjch'), return; end
    E.RejectChannels(pacienti(p).rjch);

    %vytvorim epochy
    filename = [ basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy];
    if exist(filename,'file')~=2
        disp(['no psychopy data:' pacienti(p).psychopy]);
        return; 
    end
    load(filename);
    E.ExtractEpochs(aedist,[-0.2 1.2],[-0.5 -0.2]);

    %vyradim epochy
    E.RjEpochsEpi([],0);
    E.RjEpochsEpi(30);

    E.PlotElectrode();
end
end

