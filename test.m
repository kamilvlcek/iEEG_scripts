%% MORLET WAVELETS PACIENT VT10

clear all;

pacienti = pacienti_aedist();

% pacient_cislo = 'VT10';
% pacient = pacienti(4);
% 
% % hippocamus
% % channels = [58 59 60 61 66 67];
% channels = [58];
% frequencies = [4:8];
% %frequencies = logspace(log10(2),log10(150),40);

pacient_cislo = 'VT8';
pacient = pacienti(2);

%channels = [88 89 90 91 92];
channels = [56]
frequencies = [1:0.5:4];

load(['D:\DATA\', pacient_cislo, '\aedist\', pacient.data]);

E = CMorlet(d,tabs,fs);
%E = CHilbert(d,tabs,fs)
decimatefactor = 8 ;

load(['D:\DATA\', pacient_cislo, '\', pacient.header]);

E.GetHHeader(H);
E.RejectChannels(pacient.rjch);
E.ChangeReference('b');

E.PasmoFrekvence(frequencies,channels,0, decimatefactor);

load(['D:\DATA\', pacient_cislo, '\aedist\', pacient.psychopy]);
E.ExtractEpochs(aedist,[-0.2 1.2],[-0.5 -0.2]);

load(['D:\DATA\', pacient_cislo, '\aedist\', pacient.rjepoch]);
E.RejectEpochs(RjEpoch, RjEpochCh);

E.ResponseSearch(0.1,[0 1 2]);

E.plotEpochs.channels = channels;

