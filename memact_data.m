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
data(:,1) = dataS(:,8);  % instead of soubor - delay
imissed = logical(dataS(:,7)); % find indexes of missed trials in original data
data(imissed,2) = -1;  % in column klavesa replace missed response by -1 (to be compatible with other tests, e.g. aedist), other cells = 0
data(:,3) = dataS(:,4); % spravne = accuracy for joystick response
data(data(:,3)==-1,3) = 0; % in column spravne - replace incorrect responses (-1) by 0 (to be compatible with other tests)
data(imissed,3) = 0; % add missed responses also in column spravne as 0  

if RT_corr == 0
    data(:,4) = dataS(:,6); % rt of start
else
    data(:,4) = dataS(:,5); % rt of correct
end
data(:,5) = dataS(:,3); % opakovani = block
data(:,6) = dataS(:,2); % feedback = zpetnavazba
data(:,7) = dataS(:,1); % condition = kategorie 0,1 = immed; 2,3=delayed; 0,2=same(=circle); 1,3=different(=square,triangle)
data(:,8) = dataS(:,9); % t_encod_del = exact duration of encoding phase in delayed trials

if size(U1,1) ~= size(data,1) || size(U2,1) ~= size(data,1)
    disp(['data:' num2str(size(data,1)) ' U1:' num2str(size(U1,1)) ' U2:' num2str(size(U2,1))]);
    error('different lengths of data and events, cannot be processed!');    
end

data(:,9)=U1(:,2); % timestamp of sync pulses of stimuli 
data(:,10)=U2(:,2); % timestamp of sync pulses of responses
data(:,11) = dataS(:,11); % answer_button_corr (accuracy of response to the question in diff conditions), 0 means both incorrect and missed; if answer_button_rt == NaN and answer_button_corr == 0, it's a missed response 
data(:,12) = dataS(:,12); % answer_button_rt (rt of response to the question in diff conditions)

% column names in the table
sloupce = {};
sloupce.delay=1;   % instead of soubor - time of delay in delayed conditions (variable: 3.9 - 4.1 s)
sloupce.klavesa=2; % -1 - missed reponses for the joystick response, all other cases = 0
sloupce.spravne=3; % accuracy for the joystick response: 1 - correct, 0 - both missed and incorrect 
sloupce.rt = 4; % RT of reaching the correct object or RT of start moving joystick according to the parameter RT_corr
sloupce.opakovani=5; % number of block
sloupce.zpetnavazba=6; % if feedback was given (only in the training)
sloupce.kategorie=7; % condition: 'immed_same' - 0; 'immed_diff' - 1; 'del_same' - 2; 'del_diff' - 3
sloupce.t_encod_del=8; % exact duration of encoding phase in delayed trials (it's always 2 sec, but here the precise duration is given with numbers after the decimal point to determine the precise start of the whole delayed epoch)
sloupce.ts_podnet=9; % timestamp of sync pulses of stimuli (in delayed trials - after delay, when a green cross and a word 'ted' are presented)  
sloupce.ts_odpoved=10; % timestamp of sync pulses of responses, when subject hits the correct object 
sloupce.answer_button_corr=11; % accuracy of response to the question in diff conditions: 1 - correct, 0 - both missed and incorrect
sloupce.answer_button_rt=12; % rt of response to the question in diff conditions

% text code for klavesa column and test conditions
klavesa = cell(2,2);
klavesa(1,:)={'None' -1};
klavesa(2,:)={'Responses' 0};

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