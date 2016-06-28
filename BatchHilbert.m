basedir = 'd:\eeg\motol\pacienti\';
timewindow = [-.2 1.2];

frekvence = struct;
f=1;
frekvence(f).todo = 1;
frekvence(f).freq = 4:1:8;
frekvence(f).freqname = '4-8'; % theta
f=2;
frekvence(f).todo = 1;
frekvence(f).freq = 50:5:90;
frekvence(f).freqname = '50-90'; % slow gamma
f=3;
frekvence(f).todo = 0;
frekvence(f).freq = 7:2:15;
frekvence(f).freqname = '7-15'; % alpha
f=4;
frekvence(f).todo = 0;
frekvence(f).freq = 50:10:120;
frekvence(f).freqname = '50-120'; %gamma 2
f=5;
frekvence(f).todo = 0;
frekvence(f).freq = 30:5:50;
frekvence(f).freqname = '30-50'; % gamma 1

pacienti = struct;
p = 1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p079 Plu VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_aedist.mat';
pacienti(p).header = 'P79_header.mat';
pacienti(p).aedist = 'p79_aedist.mat';
pacienti(p).rjepoch = 'aedist RjEpoch.mat';
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
pacienti(p).todo = 1;
pacienti(p).folder = 'p095 Hav VT11';
pacienti(p).data = 'VT11_2015-12-15_aedist.mat';
pacienti(p).header = 'P95_header.mat';
pacienti(p).aedist = 'p95_aedist.mat';
pacienti(p).rjepoch = 'aedist rjepoch.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=4;
pacienti(p).todo = 1;
pacienti(p).folder = 'p097 Nov VT13';
pacienti(p).data = 'VT13_2016-02-11_09-20_001 aedist.mat';
pacienti(p).header = 'P97_header.mat';
pacienti(p).aedist = 'p97_aedist.mat';
pacienti(p).rjepoch = 'aedist rj epoch.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>


for f=1:numel(frekvence)    
    disp(frekvence(f).freqname); 
    if frekvence(f).todo
        for p = 1:numel(pacienti)
            disp(pacienti(p).folder);
            if pacienti(p).todo
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
                E.ExtractEpochs(aedist,timewindow);        
                E.RejectEpochs(RjEpoch);
                E.ResponseSearch(0.1,[0 1 2]); %statistika s klouzavym oknem 100ms
                E.Save([ basedir pacienti(p).folder '\aedist CHilbert' frekvence(f).freqname ' Ep']);
                disp([ pacienti(p).folder 'OK']);
                clear E d tabs fs mults header RjEpoch aedist H ans; 
            end
        end
    end
end

