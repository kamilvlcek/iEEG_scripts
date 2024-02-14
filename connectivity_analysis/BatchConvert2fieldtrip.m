function BatchConvert2fieldtrip(typeEpochs, bipolarRef)
% convert raw preprocessed ieeg data to the fieldtrip data format for all patients

%%%   typeEpochs = to use one of four different epochs with different epoch time and baseline and alignment to stimulus:
% 0 - immediate epochs
% 1 - epochs before delay (or encoding + delay, can be changed in setup_memact)
% 2 - epochs after delay
% 3 - epochs within delay

%%% bipolarRef = if to change the original reference to bipolar
% 1 - bipolar
% 0 - original
setup = setup_memact(typeEpochs); % setup for the specific epoch type
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if(exist([basedir pacienti(p).folder '\' pacienti(p).data],'file')~=2)
            if(exist([basedir pacienti(p).folder '\' subfolder '\'  pacienti(p).data],'file')~=2)
                datafolder = ['\' subfolder];
            end
        else
            datafolder = '' ;
        end
        % load all files
        load([basedir pacienti(p).folder datafolder '\' pacienti(p).data]); % raw ieeg data
        load([basedir pacienti(p).folder '\' pacienti(p).header]); % header
        load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy]); % psychopy data
        load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).rjepoch]); % rejected epochs
        load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).epievents]); % epievents
        
        % some preprocessing, including the epochs extraction
        E = CiEEGData (d, tabs, fs); % create ieeg data object
        E.GetHHeader(H); % read header
        E.GetEpiEvents(DE); % read epi events
        E.RejectChannels(pacienti(p).rjch); % reject bad channels
        E.Filter({0.5 [50 100 150] },[],[],0); % notch filter for [50 100 150]+-.5Hz
        if bipolarRef
            E.ChangeReference('b'); % change the reference to bipolar
            ref = 'refBipo';
        else
            ref = 'refOrig';
        end
        
        E.ExtractEpochs(memact,setup_memact(typeEpochs)); % extract delayed epochs (encoding+delay) or immediate
        E.RejectEpochs(rjepoch(setup.index).RjEpoch, rjepoch(setup.index).RjEpochCh); % save rejected epochs to the object
        
        % create the fieldtrip data structure
        data = {};
        data.fsample = E.fs;
        data.trial = cell(1, size(E.d,3)); % cell 1 x ntrials
        
        % convert 3D matrix E.d (time x channels x epochs) to 2D (channels x time) in cell of epochs
        for epoch = 1:size(E.d,3)
            data.trial{epoch} = squeeze(E.d(:, :, epoch))';
        end
        
        data.time = cell(1, size(E.d,3)); % cell of all trials, each trial - 1 x time points
        for triali = 1 : size(E.d,3)
            data.time{triali} = E.epochtime(1): 1/E.fs : E.epochtime(2);
        end
        
        data.label = {E.CH.H.channels(:).name}'; % original names of all channels
        
        data.channelInfo = E.CH.H.channels; % info about channels
        
        data.RjEpochChannel = E.RjEpochCh; % channel x trials (epochs),
        % epochs in each channel containing interictal epileptiform discharges, which were identified by a spike detector
        % marked by 1 should be rejected
        
        % create a separate table with info about trials
        TrialNumber = 1:size(E.PsyData.P.data,1);
        Condition = E.PsyData.P.data(:,7);
        Correct = E.PsyData.P.data(:,3);
        ResponseTime = E.PsyData.P.data(:,4);
        Block = E.PsyData.P.data(:,5);
        Training = E.PsyData.P.data(:,6);
        DelayLength = E.PsyData.P.data(:,1);
        AnswerButtonCorrect = E.PsyData.P.data(:,11);
        AnswerButtonRT = E.PsyData.P.data(:,12);
        Trials2Reject = zeros(size(E.PsyData.P.data,1),1); % epochs that should be rejected globally (over all channels) with too many epi events (from rjepoch struct)
        Trials2Reject(E.RjEpoch) = 1;
        
        TrialInformationTable = table(TrialNumber', Condition, Correct, ResponseTime, Block, Training, DelayLength, AnswerButtonCorrect, AnswerButtonRT, Trials2Reject,...
            'VariableNames',...
            {'TrialNumber','Condition', 'Correct', 'ResponseTime', 'Block', 'Training', 'DelayLength', 'AnswerButtonCorrect', 'AnswerButtonRT','Trials2Reject'} );
        
        % save data
        outfilename = [ basedir pacienti(p).folder '\' subfolder '\' setup.prefix '_' ref ' ' sprintf('%.1f-%.1f',setup.epochtime(1:2)) ' ' setup.suffix ' fieldtrip_' datestr(now,'YYYY-mm') '.mat'];
        save(outfilename, 'data', 'TrialInformationTable');
        clear E d tabs fs mults header RjEpoch RjEpochCh memact H data TrialInformationTable; 
    end
end
end