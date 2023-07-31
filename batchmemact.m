%% which type of epoch I want to analyze
immediate = 0; 
beforedelay = 1;  
afterdelay = 1; 
withindelay = 0;

%% an analysis of immediate epochs
 if immediate
    disp(' ++++ ANALYSIS 1 - immediate epochs ++++');
    pacienti = {'VT66'}; 
    cfg = struct('hybernovat',0,'suffix','Ep2023-07'); 
    cfg.pacienti = pacienti; %which patients to analyse - has to be added after struct creation
    %cfg.overwrite=1; %overwrite old output files?
    %cfg.freqepochs=1; %save all frequencies from all epochs?
    %cfg.normalization='z'; %use different normalization than orig (orig =divide by mean)
    %cfg.epochfilter = {7,[0 1]}; %which epochs to extract and save - {column, [values]}; %stored already in setup
    cfg.typeEpochs = 0; % immediate
    cfg.normalizeEpochs = 1;
    filenames = BatchHilbert('memact',cfg);
end
%% an analysis of epochs before delay
if beforedelay
    disp(' ++++ ANALYSIS 2 - epochs before delay ++++');
%     pacienti = {'VT66'};
    cfg = struct('hybernovat',0,'suffix','Ep2023-07'); 
%     cfg.pacienti = pacienti;
    cfg.typeEpochs = 1; % before delay
    cfg.normalizeEpochs = 0; % we won't normalize parts of delayed epochs before connecting them together     
    filenames = BatchHilbert('memact',cfg);
end
%% an analysis of epochs after delay
if afterdelay
    disp(' ++++ ANALYSIS 3 - epochs after delay ++++');
%     pacienti = {'VT66'};
    cfg = struct('hybernovat',0,'suffix','Ep2023-07');
%     cfg.pacienti = pacienti;
    cfg.typeEpochs = 2; % after delay
    cfg.normalizeEpochs = 0;    
    filenames = BatchHilbert('memact',cfg);
end
%% an analysis of epochs during the whole delay
if withindelay
    disp(' ++++ ANALYSIS 4 - epochs within delay ++++');
    cfg = struct('hybernovat',0,'suffix','Ep2023-07'); 
    cfg.typeEpochs = 3; % within delay
    cfg.normalizeEpochs = 0;    
    filenames = BatchHilbert('memact',cfg);
end