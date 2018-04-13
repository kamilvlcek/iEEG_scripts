%% nejdriv normalni analyzu s razenim podle podnetu
disp(' ++++ ANALYZA 1 - RAZENI PODLE PODNETU ++++');
cfg = struct('hybernovat',0);
BatchHilbert('aedist',cfg);

%% potom analyza s razenim podle odpovedi
disp(' ++++ ANALYZA 2 - RAZENI PODLE ODPOVEDI ++++');
cfg = struct('hybernovat',0,'srovnejresp',1);
BatchHilbert('aedist',cfg);

%% nakonec analyza s podilem casu odpovedi
% uz budu na konci hybernovat
disp(' ++++ ANALYZA 3 - PODIL CASU ODPOVEDI ++++');
cfg = struct('hybernovat',1,'podilcasuodpovedi',1);
BatchHilbert('aedist',cfg);