%% nejdriv normalni analyzu s razenim podle podnetu
disp(' ++++ ANALYZA 1 - RAZENI PODLE PODNETU ++++');
%cfg = struct('hybernovat',0); 
pacienti = {'p073','p097','p132','p136'}; 
cfg = struct('hybernovat',0,'freqepochs',1); %frekvence vsech epoch
cfg.pacienti = pacienti; %kdyz to tam vlozim rovnou, tak se mi udela struct array
BatchHilbert('ppa',cfg);

%dalsi analyza uz nema smysl
%% potom analyza s razenim podle odpovedi
% disp(' ++++ ANALYZA 2 - RAZENI PODLE ODPOVEDI ++++');
% cfg = struct('hybernovat',0,'srovnejresp',1);
% BatchHilbert('menrot',cfg);

%% nakonec analyza s podilem casu odpovedi
% uz budu na konci hybernovat
% disp(' ++++ ANALYZA 3 - PODIL CASU ODPOVEDI ++++');
% cfg = struct('hybernovat',1,'podilcasuodpovedi',1);
% BatchHilbert('menrot',cfg);