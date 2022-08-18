function [memact] = memact_data(pacientid,RT_corr, U1,U2,tabs,eegfile)
%%%% creates and returns the behavioral data structure of the psychopy 
%%%% that contains all trials and their timestamps in MemoryActions test 
% pacientid - id pacienta, e.g. p85
% RT_corr - if to use behavioral RT of reaching the correct object (=1) or RT of start moving joystick (0)
% U1 and U2 - data from the function udalosti2() - timestamps of sync pulses of stimuli and responses (U2 = RT_corr)
% tabs - tabs value from data from Motol, possibly truncated with datatrim()
% eegfile - the original file name with the EEG data

if ~exist('RT_cor','var')  || isempty(RT_cor)
    RT_corr = 0; % default: RT of start moving joystick, if 1 - RT of hiting the correct object 
end

dir = 'd:\eeg\motol\PsychoPydata\MemoryActions\';
% load mat file with all behav data
load([dir pacientid '_MemoryActions.mat']);

% get only table with RT, accuracy and timing, without coordinates of joystick 
dataS = MemoryActions.Gdata;

% rearrange data according to Kamil's format (like in aedist_data.m)
data(:,1) = dataS(:,7);  % instead of soubor - delay
data(:,2) = dataS(:,10); % answer_button (to the question in diff conditions)
data(:,3) = dataS(:,4); % spravne = accuracy
if RT_corr == 0
    data(:,4) = dataS(:,6); % rt of start
else
    data(:,4) = dataS(:,5); % rt of correct
end
data(:,5) = dataS(:,3); % opakovani = block
data(:,6) = dataS(:,2); % feedback = zpetnavazba
data(:,7) = dataS(:,1); % condition = kategorie
data(:,8) = dataS(:,8); % t_encod_del = exact duration of encoding phase in delayed trials

if size(U1,1) ~= size(data,1) || size(U2,1) ~= size(data,1)
    disp(['data:' num2str(size(data,1)) ' U1:' num2str(size(U1,1)) ' U2:' num2str(size(U2,1))]);
    error('different lengths of data and events, cannot be processed!');    
end

data(:,9)=U1(:,2); % timestamp of sync pulses of stimuli 
data(:,10)=U2(:,2); % timestamp of sync pulses of responses

% column names in the table
sloupce = {};
sloupce.delay=1;   % instead of soubor - time of delay in delayed conditions
sloupce.klavesa=2; % answer_button response (to the question in diff conditions)
sloupce.spravne=3;
sloupce.rt = 4;
sloupce.opakovani=5;
sloupce.zpetnavazba=6;
sloupce.kategorie=7;
sloupce.t_encod_del=8;
sloupce.ts_podnet=9;
sloupce.ts_odpoved=10;

% text code for button response and test conditions
klavesa = cell(3,2);
klavesa(1,:)={'incorrect' 0};
klavesa(2,:)={'correct' 1};
klavesa(3,:)={'none' NaN};

podminka = cell(4,2);
podminka(1,:)={'immed_same' 0};
podminka(2,:)={'immed_diff' 1};
podminka(3,:)={'del_same' 2};
podminka(4,:)={'del_diff' 3};

memact = struct('data',data,'sloupce',sloupce);
memact.strings.klavesa = klavesa;
memact.strings.podminka = podminka; 

% timestamps of the start and end of the data from the test
memact.interval = [tabs(1) tabs(end)];
disp(['MemActions data od ' datestr(tabs(1),'dd-mmm-yyyy HH:MM:SS.FFF') ' do ' datestr(tabs(end),'dd-mmm-yyyy HH:MM:SS.FFF')]);

memact.eegfile = eegfile;
memact.pacientid = pacientid; 
end