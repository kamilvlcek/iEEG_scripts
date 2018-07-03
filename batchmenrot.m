%% nejdriv normalni analyzu s razenim podle podnetu
disp(' ++++ ANALYZA 1 - RAZENI PODLE PODNETU ++++');
%pacienti = {'p082'}; 
cfg = struct('hybernovat',0,'suffix','Ep2018-01');
%cfg.pacienti = pacienti; %kdyz to tam vlozim rovnou, tak se mi udela struct array
%cfg.overwrite=1; %vyjimecne
BatchHilbert('menrot',cfg);
return; %nic dalsiho zatim nechci 
%% potom analyza s razenim podle odpovedi
disp(' ++++ ANALYZA 2 - RAZENI PODLE ODPOVEDI ++++');
cfg = struct('hybernovat',0,'srovnejresp',1); %,'suffix','Ep2018-01'
BatchHilbert('menrot',cfg);

%% nakonec analyza s podilem casu odpovedi
% uz budu na konci hybernovat
disp(' ++++ ANALYZA 3 - PODIL CASU ODPOVEDI ++++');
cfg = struct('hybernovat',1,'podilcasuodpovedi',1); %,'suffix','Ep2018-01'
BatchHilbert('menrot',cfg);