%% nejdriv normalni analyzu s razenim podle podnetu
disp(' ++++ ANALYZA 1 - RAZENI PODLE PODNETU ++++');
pacienti = {'p160'}; 
cfg = struct('hybernovat',0,'suffix','Ep2018-04');
cfg.pacienti = pacienti; %kdyz to tam vlozim rovnou, tak se mi udela struct array
BatchHilbert('aedist',cfg);
return; %nic dalsiho zatim nechci 
%% potom analyza s razenim podle odpovedi
disp(' ++++ ANALYZA 2 - RAZENI PODLE ODPOVEDI ++++');
cfg = struct('hybernovat',0,'srovnejresp',1);
BatchHilbert('aedist',cfg);

%% nakonec analyza s podilem casu odpovedi
% uz budu na konci hybernovat
disp(' ++++ ANALYZA 3 - PODIL CASU ODPOVEDI ++++');
cfg = struct('hybernovat',1,'podilcasuodpovedi',1);
BatchHilbert('aedist',cfg);