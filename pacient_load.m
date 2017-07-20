function [ E ] = pacient_load( nick,test )
%PACIENT_LOAD Summary of this function goes here
%   Detailed explanation goes here
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

