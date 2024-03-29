%% ktere analyzy chci spustit
podnet = 1;
odpovedi = 0;
podilcasuodpovedi =0;

%% nejdriv normalni analyzu s razenim podle podnetu
if podnet
    disp(' ++++ ANALYZA 1 - RAZENI PODLE PODNETU ++++');
    pacienti = {'p181'}; 
    cfg = struct('hybernovat',0,'suffix','Ep2023-06'); %zatim necham puvodni priponu aby stejne soubory
    cfg.pacienti = pacienti; %kdyz to tam vlozim rovnou, tak se mi udela struct array
    cfg.overwrite=1; %vyjimecne
    cfg.decimatefactor = 1; %for single trial analysis. We do not want to reduce the sampling frequency
    filenames = BatchHilbert('menrot',cfg);    
end
%% potom analyza s razenim podle odpovedi
if odpovedi
    disp(' ++++ ANALYZA 2 - RAZENI PODLE ODPOVEDI ++++');
    cfg = struct('hybernovat',0,'srovnejresp',1,'suffix','Ep2018-07'); %,'suffix','Ep2018-01'
    cfg.overwrite=0; %vyjimecne
    filenames = BatchHilbert('menrot',cfg);
end
%% nakonec analyza s podilem casu odpovedi
% uz budu na konci hybernovat
if podilcasuodpovedi
    disp(' ++++ ANALYZA 3 - PODIL CASU ODPOVEDI ++++');
    cfg = struct('hybernovat',0,'podilcasuodpovedi',1); %,'suffix','Ep2018-01'
    filenames = BatchHilbert('menrot',cfg);
end