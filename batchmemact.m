%% which analyses I want to run
podnet = 1; %stimulus-based align of epochs
odpovedi = 0;  %response align
podilcasuodpovedi =0; %response time percentage

%% a normal analysis with a stimulus-based align of epochs
 if podnet
    disp(' ++++ ANALYZA 1 - RAZENI PODLE PODNETU ++++');
    pacienti = {'VT59','VT61'}; 
    cfg = struct('hybernovat',0,'suffix','Ep2023-06'); %'Ep2018-04', 'Ep2018-06'
    cfg.pacienti = pacienti; %which patients to analyse - has to be added after struct creation
    %cfg.overwrite=1; %overwrite old output files?
    %cfg.freqepochs=1; %save all frequencies from all epochs?
    %cfg.normalization='z'; %use different normalization than orig (orig =divide by mean)
    cfg.epochfilter = {7,[0 1]}; %which epochs to extract and save - {column, [values]}; only immediate epochs
    filenames = BatchHilbert('memact',cfg);
end
%% an analysis with aligning according to the responses
if odpovedi
    disp(' ++++ ANALYZA 2 - RAZENI PODLE ODPOVEDI ++++');
    cfg = struct('hybernovat',0,'srovnejresp',1);
    filenames = BatchHilbert('aedist',cfg);
end
%% finally, the analysis with response time percentage
% hibernating at the end
if podilcasuodpovedi
    disp(' ++++ ANALYZA 3 - PODIL CASU ODPOVEDI ++++');
    cfg = struct('hybernovat',1,'podilcasuodpovedi',1);
    filenames = BatchHilbert('aedist',cfg);
end