%% TEST SCRIPT

clear all;
pacienti = pacienti_aedist();
location = 'D:\DATA\';

%%%%% p073 Pech VT6
%     pacient_cislo = 1;
%%%%% p079 Plu VT8
    pacient_cislo = 2;
    channels = [50 51];
%%%%% p082 Vov VT9
%     pacient_cislo = 3;
%%%%% p083 Kol VT10
%     pacient_cislo = 4;
%     channels = [58 59 60 61 66 67];
%%%%% p095 Hav VT11
%     pacient_cislo = 5;
%     channels = [88 89 90 91 92];
%%%%% p096 Gro VT12
%     pacient_cislo = 6;
%     channels = [20 21 22 23 24 31 32 33 66 67 68 69 77 78 79];

pacient = pacienti(pacient_cislo);

%frequencies = logspace(log10(2),log10(150),40);

frequencies = [1:0.5:4];

load([location, pacient.folder, '\aedist\', pacient.data]);

E = CMorlet(d,tabs,fs);
%E = CHilbert(d,tabs,fs)
decimatefactor = 8 ;

load([location, pacient.folder, '\', pacient.header]);

E.GetHHeader(H);
E.RejectChannels(pacient.rjch);
E.ChangeReference('b');

E.PasmoFrekvence(frequencies,channels,0, decimatefactor);

load([location, pacient.folder, '\aedist\', pacient.psychopy]);
E.ExtractEpochs(aedist,[-0.2 1.2],[-0.5 -0.2]);

load([location, pacient.folder, '\aedist\', pacient.rjepoch]);
E.RejectEpochs(RjEpoch, RjEpochCh);

E.ResponseSearch(0.1,[0 1 2]);

E.plotEpochs.channels = channels;


%% VYPISANIE NAZVOV BIP. CHANNELS

for i=1:length(E.CH.H.channels)
    display([num2str(i), ' channels = ', E.CH.H.channels(i).name])
end

