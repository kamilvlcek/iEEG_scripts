pacients = {  'p079'}; %'p079' 'p083' 'p095' 'p097'
freq = 10:10:150; freqname = '10-150';
basedir = 'd:\eeg\motol\pacienti\';
if sum(ismember(pacients,'p079'))>0
    disp('p079 Plu VT8');
    load([basedir 'p079 Plu VT8\VT8_2015-04-09_09-46_001_concat_X_aedist.mat']);
    load([basedir 'p079 Plu VT8\P79_header.mat']);
    load([basedir 'p079 Plu VT8\p79_aedist.mat']);
    E = CHilbert(d,tabs,fs,[],header);
    clear d;
    E.GetHHeader(H);
    E.RejectChannels([47 64 114]);
    E.ChangeReference('b');
    E.PasmoFrekvence(freq);
    E.ExtractEpochs(aedist,[-.2 1.2]);
    load([basedir 'p079 Plu VT8\aedist RjEpoch.mat']);
    E.RejectEpochs(RjEpoch);
    E.Save([ basedir 'p079 Plu VT8\aedist CHilbert' freqname ' Ep']);
    disp('p079 OK');
    clear E d tabs fs mults header RjEpoch aedist H ans;
end

if sum(ismember(pacients,'p083'))>0
    %freq = 1:2:9; freqname = '1-9';
    disp('p083 Kol VT10');
    load([basedir 'p083 Kol VT10\VT10_2015-05-19_10-00_001_X_aedist.mat']);
    load([basedir 'p083 Kol VT10\P83_header.mat']);
    load([basedir 'p083 Kol VT10\p83_aedist.mat']);
    E = CHilbert(d,tabs,fs,[],header);
    clear d;
    E.GetHHeader(H);
    E.RejectChannels([47]);
    E.ChangeReference('b');
    E.PasmoFrekvence(freq);
    E.ExtractEpochs(aedist,[-.2 1.2]);
    load([basedir 'p083 Kol VT10\aedist RjEpoch.mat']);
    E.RejectEpochs(RjEpoch);
    E.Save([ basedir 'p083 Kol VT10\aedist CHilbert' freqname ' Ep']);
    disp('p083 OK');
    clear E d tabs fs mults header RjEpoch aedist H ans;
end


if sum(ismember(pacients,'p095'))>0
    disp('p095 Hav VT11');
    load([basedir 'p095 Hav VT11\VT11_2015-12-15_aedist.mat']);
    load([basedir 'p095 Hav VT11\P95_header.mat']);
    load([basedir 'p095 Hav VT11\p95_aedist.mat']);
    E = CHilbert(d,tabs,fs,[],header);
    clear d;
    E.GetHHeader(H);
    E.RejectChannels([47]);
    E.ChangeReference('b');
    E.PasmoFrekvence(freq);
    E.ExtractEpochs(aedist,[-.2 1.2]);
    load([basedir 'p095 Hav VT11\aedist rjepoch.mat']);
    E.RejectEpochs(RjEpoch);
    E.Save([basedir 'p095 Hav VT11\aedist CHilbert' freqname ' Ep']);
    disp('p095 OK');
    clear E d tabs fs mults header RjEpoch aedist H ans;
end

if sum(ismember(pacients,'p097'))>0
    disp('p097 Nov VT13');
    load([basedir 'p097 Nov VT13\VT13_2016-02-11_09-20_001 aedist.mat']);
    load([basedir 'p097 Nov VT13\P97_header.mat']);
    load([basedir 'p097 Nov VT13\p97_aedist.mat']);
    E = CHilbert(d,tabs,fs,mults);
    clear d;
    E.GetHHeader(H);
    E.RejectChannels([47]);
    E.ChangeReference('b');
    E.PasmoFrekvence(freq);
    E.ExtractEpochs(aedist,[-.2 1.2]);
    load([basedir 'p097 Nov VT13\aedist rj epoch.mat']);
    E.RejectEpochs(RjEpoch);
    E.Save([basedir 'p097 Nov VT13\aedist CHilbert' freqname ' Ep']);
    disp('p097 OK');
    clear E d tabs fs mults header RjEpoch aedist H ans;
end