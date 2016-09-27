%15.9.2016 - AlloEgo zarovnani podle odpovedi
basedir = 'd:\eeg\motol\pacienti\';
timewindow = [-1 1]; %zarovnani podle odpovedi/podnetu
baseline = [-1 0.8]; %baseline
suffix = ' Ep';
prefix = 'AlloEgo';

frekvence = struct;
f=1;
frekvence(f).todo = 0;
frekvence(f).freq = 4:2:150;
frekvence(f).freqname = '4-150'; % all range
f=2;
frekvence(f).todo = 1;
frekvence(f).freq = 50:5:150;
frekvence(f).freqname = '50-150'; % slow gamma
f=3;
frekvence(f).todo = 1;
frekvence(f).freq = 7:2:15;
frekvence(f).freqname = '7-15'; % alpha
f=4;
frekvence(f).todo = 0;
frekvence(f).freq = 50:10:90;
frekvence(f).freqname = '50-90'; %gamma 2
f=5;
frekvence(f).todo = 1;
frekvence(f).freq = 30:5:50;
frekvence(f).freqname = '30-50'; % gamma 1
f=6;
frekvence(f).todo = 1;
frekvence(f).freq = 15:3:31;
frekvence(f).freqname = '15-31'; % gamma 1
f=7;
frekvence(f).todo = 1;
frekvence(f).freq = 4:1:8;
frekvence(f).freqname = '4-8'; % theta

pacienti = struct;
p = 1;
pacienti(p).todo = 0;
pacienti(p).folder = 'p079 Plu VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_aedist.mat';
pacienti(p).header = 'P79_header.mat';
pacienti(p).aedist = 'p79_aedist.mat';
pacienti(p).rjepoch = 'aedist RjEpoch Resp.mat';
pacienti(p).rjch = [47 64 114]; %#ok<NBR%#ok<MSNU> AK> 
pacienti(p).frekvence = struct;

p=2;
pacienti(p).todo = 0;
pacienti(p).folder = 'p083 Kol VT10';
pacienti(p).data = 'VT10_2015-05-19_10-00_001_X_aedist.mat';
pacienti(p).header = 'P83_header.mat';
pacienti(p).aedist = 'p83_aedist.mat';
pacienti(p).rjepoch = 'aedist RjEpoch.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=3;
pacienti(p).todo = 0;
pacienti(p).folder = 'p095 Hav VT11';
pacienti(p).data = 'VT11_2015-12-15_aedist.mat';
pacienti(p).header = 'P95_header.mat';
pacienti(p).aedist = 'p95_aedist.mat';
pacienti(p).rjepoch = 'aedist rjepoch Resp.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=4;
pacienti(p).todo = 0;
pacienti(p).folder = 'p097 Nov VT13';
pacienti(p).data = 'VT13_2016-02-11_09-20_001 aedist.mat';
pacienti(p).header = 'P97_header.mat';
pacienti(p).aedist = 'p97_aedist.mat';
pacienti(p).rjepoch = 'aedist rj epoch.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=5;
pacienti(p).todo = 0;
pacienti(p).folder = 'p073 Pech VT6';
pacienti(p).data = 'VT6_INV Test Vlcek1_X_aedist.mat';
pacienti(p).header = 'p73_header_kamil.mat';
pacienti(p).aedist = 'p73_aedist.mat';
pacienti(p).rjepoch = 'aedist_RjEpoch.mat';
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=6;
pacienti(p).todo = 1;
pacienti(p).folder = 'p082 Vov VT9';
pacienti(p).data = 'VT9_2015-04-23_11-00_001_X_500hz_concat BVA.mat';
pacienti(p).header = 'P82_64ch_header.mat';
pacienti(p).aedist = 'p82_alloego.mat';
pacienti(p).rjepoch = 'alloego RjEpoch.mat';
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=7;
pacienti(p).todo = 1;
pacienti(p).folder = 'p083 Kol VT10';
pacienti(p).data = 'VT10_2015-06-02_13-39_001_X_500hz_concat BVA.mat';
pacienti(p).header = 'P83_64ch_header.mat';
pacienti(p).aedist = 'p83_alloego.mat';
pacienti(p).rjepoch = 'alloego RjEpoch.mat';
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

%nejdriv overim, jestli existuje vsechno co potrebuju nacist
for p = 1:numel(pacienti)
    if pacienti(p).todo 
        assert(exist([basedir pacienti(p).folder '\' pacienti(p).data],'file')==2,      ['Data neexistuji: ' pacienti(p).folder]);
        assert(exist([basedir pacienti(p).folder '\' pacienti(p).header],'file')==2,    ['Header neexistuje: ' pacienti(p).folder]);
        assert(exist([basedir pacienti(p).folder '\' pacienti(p).aedist],'file')==2,    ['Aedist neexistuje: ' pacienti(p).folder]);
        assert(exist([basedir pacienti(p).folder '\' pacienti(p).rjepoch],'file')==2,   ['rjepoch neexistuje: ' pacienti(p).folder]);
    end
end
disp('vsechny soubory ok');
 
for f=1:numel(frekvence)        
    if frekvence(f).todo
        disp(frekvence(f).freqname); 
        for p = 1:numel(pacienti)            
            if pacienti(p).todo
                disp(pacienti(p).folder);
                load([basedir pacienti(p).folder '\' pacienti(p).data]);
                load([basedir pacienti(p).folder '\' pacienti(p).header]);
                load([basedir pacienti(p).folder '\' pacienti(p).aedist]);
                load([basedir pacienti(p).folder '\' pacienti(p).rjepoch]);
                if ~exist('mults','var'),  mults = []; end
                if ~exist('header','var'), header = []; end
                E = CHilbert(d,tabs,fs,mults,header);
                clear d;
                E.GetHHeader(H);
                E.RejectChannels(pacienti(p).rjch);
                E.ChangeReference('b');
                E.PasmoFrekvence(frekvence(f).freq);
                disp('extracting epochs ...');
                E.ExtractEpochs(alloego,timewindow,baseline);        
                E.RejectEpochs(RjEpoch);
                E.ResponseSearch(0.1,[0]); %statistika s klouzavym oknem 100ms
                disp('saving data ...');
                E.Save([ basedir pacienti(p).folder '\' prefix ' CHilbert ' frekvence(f).freqname ' ' suffix]);
                disp([ pacienti(p).folder ' OK']);
                clear E d tabs fs mults header RjEpoch aedist H ans; 
            end
        end
    end
end

system('shutdown -s')
