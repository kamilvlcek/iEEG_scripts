%% ktere analyzy chci spustit
podnet = 1;
odpovedi = 0;
podilcasuodpovedi =0;

%% nejdriv normalni analyzu s razenim podle podnetu
 if podnet
    disp(' ++++ ANALYZA 1 - RAZENI PODLE PODNETU ++++');
    pacienti = {'p073'}; 
    cfg = struct('hybernovat',0,'suffix','Ep2018-07'); %'Ep2018-04', 'Ep2018-06'
    cfg.pacienti = pacienti; %kdyz to tam vlozim rovnou, tak se mi udela struct array
    %cfg.overwrite=1; %vyjimecne
    %cfg.freqepochs=1; %vyjimecne
    filenames = BatchHilbert('aedist',cfg);
end
%% potom analyza s razenim podle odpovedi
if odpovedi
    disp(' ++++ ANALYZA 2 - RAZENI PODLE ODPOVEDI ++++');
    cfg = struct('hybernovat',0,'srovnejresp',1);
    filenames = BatchHilbert('aedist',cfg);
end
%% nakonec analyza s podilem casu odpovedi
% uz budu na konci hybernovat
if podilcasuodpovedi
    disp(' ++++ ANALYZA 3 - PODIL CASU ODPOVEDI ++++');
    cfg = struct('hybernovat',1,'podilcasuodpovedi',1);
    filenames = BatchHilbert('aedist',cfg);
end