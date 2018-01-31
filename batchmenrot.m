%% nejdriv normalni analyzu s razenim podle podnetu
cfg = struct('hybernovat','0');
BatchHilbert('menrot',cfg);

%% potom analyza s razenim podle odpovedi
cfg = struct('hybernovat','0','srovnejresp',1);
BatchHilbert('menrot',cfg);

%% nakonec analyza s podilem casu odpovedi
% uz budu na konci hybernovat
cfg = struct('hybernovat','1','podilcasuodpovedi',1);
BatchHilbert('menrot',cfg);