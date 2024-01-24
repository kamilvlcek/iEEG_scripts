%% which type of epoch I want to analyze
immediate = 1; 
beforedelay = 0;  
afterdelay = 0; 
withindelay = 1;
baseline4delay = 1;

%% an analysis of immediate epochs
 if immediate
    disp(' ++++ ANALYSIS 1 - immediate epochs ++++');
    pacienti = {'VT62'}; 
    cfg = struct('hybernovat',0,'suffix','Ep2024-01'); 
    cfg.pacienti = pacienti; %which patients to analyse - has to be added after struct creation
    %cfg.overwrite=1; %overwrite old output files?
    cfg.freqepochs=1; %save all frequencies from all epochs?
    %cfg.normalization='z'; %use different normalization than orig (orig =divide by mean)
    %cfg.epochfilter = {7,[0 1]}; %which epochs to extract and save - {column, [values]}; %stored already in setup
    cfg.typeEpochs = 0; % immediate
    cfg.normalizeEpochs = 1;
    filenames = BatchHilbert('memact',cfg);
end
%% an analysis of epochs before delay
if beforedelay
    disp(' ++++ ANALYSIS 2 - epochs before delay ++++');
%     pacienti = {'VT67', 'VT68'};
    cfg = struct('hybernovat',0,'suffix','Ep2024-01'); 
%     cfg.pacienti = pacienti;
    cfg.typeEpochs = 1; % before delay
%     cfg.normalizeEpochs = 0; % we won't normalize parts of delayed epochs and compute statistics in BatchHilbert before connecting two parts together     
    cfg.normalizeEpochs = 1;
    filenames = BatchHilbert('memact',cfg);
end
%% an analysis of epochs after delay
if afterdelay
    disp(' ++++ ANALYSIS 3 - epochs after delay ++++');
%     pacienti = {'VT66'};
    cfg = struct('hybernovat',0,'suffix','Ep2024-01');
%     cfg.pacienti = pacienti;
    cfg.typeEpochs = 2; % after delay
    cfg.normalizeEpochs = 0;    
    filenames = BatchHilbert('memact',cfg);
end
%% an analysis of epochs during the whole delay
if withindelay
    disp(' ++++ ANALYSIS 4 - epochs within delay ++++');
    pacienti = {'VT62'};
    cfg = struct('hybernovat',0,'suffix','Ep2024-01'); 
    cfg.pacienti = pacienti;
    cfg.freqepochs=1;
    cfg.typeEpochs = 3; % within delay
    cfg.normalizeEpochs = 0; % doesn't compute statistics in BatchHilbert, we should compute it only after appending with the baseline
    filenames = BatchHilbert('memact',cfg);
end
%% extraction of baseline activity for analysis of epochs during the whole delay
if baseline4delay
    disp(' ++++ ANALYSIS 4 - extraction of baseline activity for epochs within delay ++++');
    pacienti = {'VT62'};
    cfg = struct('hybernovat',0,'suffix','Ep2024-01'); 
    cfg.pacienti = pacienti;
    cfg.freqepochs=1;
    cfg.typeEpochs = 4; % baseline activity before encoding phase
    cfg.normalizeEpochs = 0;   
    filenames = BatchHilbert('memact',cfg);
end
