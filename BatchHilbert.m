%15.9.2016 - AlloEgo zarovnani podle odpovedi
%25.5.2017 - Pridani reference a ERP
hybernovat = 0; %jestli chci po konci skriptu pocitac uspat - ma prednost
vypnout = 0; %jestli chci po konci skriptu pocitac vypnout (a nechci ho hybernovat) 

basedir = 'd:\eeg\motol\pacienti\';
timewindow =  [-0.3 0.8];  % hranice epochy [-0.3 0.8] PPA, zarovnani podle odpovedi/podnetu [-1 1];
baseline = []; %baseline [-1 0.8]
suffix = 'Ep';
prefix = 'PPA'; %musi byt bud AlloEgo, PPA, AEdist
stat_kats = [2 3 1];  % PPA [2 3 1] Face, Object, Scene ; AlloEgo [0 1 2] Control, Ego, Allo; 
stat_opak = {}; %{[1 2],[4 5]}; %PPA opakovani 12 vs 45
subfolder = 'PPA'; %podadresar, specificky pro test, muze byt prazdne pokud se nepouzivaji podadresare

frekvence = struct;
f=1;
frekvence(f).todo = 1;
frekvence(f).freq = [];
frekvence(f).freqname = 'ERP'; % slow gamma
f=2;
frekvence(f).todo = 1;
frekvence(f).freq = 50:5:150;
frekvence(f).freqname = '50-150'; % slow gamma
f=3;
frekvence(f).todo = 0;
frekvence(f).freq = 7:2:15;
frekvence(f).freqname = '7-15'; % alpha
f=4;
frekvence(f).todo = 0;
frekvence(f).freq = 50:10:90;
frekvence(f).freqname = '50-90'; %gamma 2
f=5;
frekvence(f).todo = 0;
frekvence(f).freq = 30:5:50;
frekvence(f).freqname = '30-50'; % gamma 1
f=6;
frekvence(f).todo = 0;
frekvence(f).freq = 15:3:31;
frekvence(f).freqname = '15-31'; % gamma 1
f=7;
frekvence(f).todo = 0;
frekvence(f).freq = 4:1:8;
frekvence(f).freqname = '4-8'; % theta
f=8;
frekvence(f).todo = 0;
frekvence(f).freq = 4:2:150;
frekvence(f).freqname = '4-150'; % all range

refence = struct;
r=1;
reference(r).todo = 1;
reference(r).name = 'refOrig';
reference(r).char = '';
r=2;
reference(r).todo = 1;
reference(r).name = 'refEle';
reference(r).char = 'e';
r=3;
reference(r).todo = 1;
reference(r).name = 'refHead';
reference(r).char = 'h';
r=4;
reference(r).todo = 1;
reference(r).name = 'refBipo';
reference(r).char = 'b';

pacienti = struct;
p = 1;
pacienti(p).todo = 0;
pacienti(p).folder = 'p079 Plu VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_ppa.mat'; %'VT8_2015-04-09_09-46_001_concat_X_aedist.mat';
pacienti(p).header = 'P79_header.mat';
pacienti(p).psychopy = 'p79_ppa.mat'; %'p79_aedist.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p79.mat'; %'aedist RjEpoch Resp.mat';
pacienti(p).rjch = [47 64 114]; %#ok<NBR%#ok<MSNU> AK> 
pacienti(p).frekvence = struct;

p=2;
pacienti(p).todo = 0;
pacienti(p).folder = 'p083 Kol VT10';
pacienti(p).data = 'VT10_2015-05-19_10-00_001_X_aedist.mat';
pacienti(p).header = 'P83_header.mat';
pacienti(p).psychopy = 'p83_aedist.mat';
pacienti(p).rjepoch = 'aedist RjEpoch.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=3;
pacienti(p).todo = 0;
pacienti(p).folder = 'p095 Hav VT11';
pacienti(p).data = 'VT11_2015-12-15_aedist.mat';
pacienti(p).header = 'P95_header.mat';
pacienti(p).psychopy = 'p95_aedist.mat';
pacienti(p).rjepoch = 'aedist rjepoch Resp.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=4;
pacienti(p).todo = 0;
pacienti(p).folder = 'p097 Nov VT13';
pacienti(p).data = 'VT13_2016-02-11_09-20_001 aedist.mat';
pacienti(p).header = 'P97_header.mat';
pacienti(p).psychopy = 'p97_aedist.mat';
pacienti(p).rjepoch = 'aedist rj epoch.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=5;
pacienti(p).todo = 1;
pacienti(p).folder = 'p073 Pech VT6';
pacienti(p).data = 'VT6_INV Test Vlcek1_X_ppa.mat'; %'VT6_INV Test Vlcek1_X_aedist.mat';
pacienti(p).header = 'p73_header.mat'; % p73_header.mat je ten nejnovejsi od Jirky - 26.5.2017 'p73_header_kamil.mat';
pacienti(p).psychopy = 'p73_ppa.mat'; %'p73_aedist.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat'; %'aedist_RjEpoch.mat';
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=6;
pacienti(p).todo = 0;
pacienti(p).folder = 'p082 Vov VT9';
pacienti(p).data = 'VT9_2015-04-23_11-00_001_X_500hz_concat BVA.mat';
pacienti(p).header = 'P82_64ch_header.mat';
pacienti(p).psychopy = 'p82_alloego.mat';
pacienti(p).rjepoch = 'alloego RjEpoch.mat';
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=7;
pacienti(p).todo = 0;
pacienti(p).folder = 'p083 Kol VT10';
pacienti(p).data = 'VT10_2015-06-02_13-39_001_X_500hz_concat BVA.mat';
pacienti(p).header = 'P83_64ch_header.mat';
pacienti(p).psychopy = 'p83_alloego.mat';
pacienti(p).rjepoch = 'alloego RjEpoch.mat';
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

%nejdriv overim, jestli existuje vsechno co potrebuju nacist
for p = 1:numel(pacienti)
    if pacienti(p).todo 
        assert(exist([basedir pacienti(p).folder '\' pacienti(p).data],'file')==2,      ['Data neexistuji: ' pacienti(p).folder '\' pacienti(p).data]);
        assert(exist([basedir pacienti(p).folder '\' pacienti(p).header],'file')==2,    ['Header neexistuje: ' pacienti(p).folder '\' pacienti(p).header]);
        assert(exist([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy],'file')==2,  ['Psychopy soubor neexistuje: ' pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy]);
        assert(exist([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).rjepoch],'file')==2,   ['rjepoch neexistuje: ' pacienti(p).folder '\' subfolder '\' pacienti(p).rjepoch]);
    end
end
disp('vsechny soubory ok');
 
for f=1:numel(frekvence)        
    if frekvence(f).todo
        disp([ '*****' frekvence(f).freqname '*****' ]); 
        if numel(frekvence(f).freq) == 0
            ERP = 1; 
        else
            ERP = 0; 
        end
        for p = 1:numel(pacienti)            
            if pacienti(p).todo
                disp( [ ' ---- ' pacienti(p).folder ' ---- ']);
                for r = 1:numel(reference)
                    if reference(r).todo
                        disp( [ ' ..... ' reference(r).name ' ..... ']);
                        load([basedir pacienti(p).folder '\' pacienti(p).data]);
                        load([basedir pacienti(p).folder '\' pacienti(p).header]);
                        load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy]);
                        if strcmp(prefix,'PPA')
                            psychopy = ppa; clear ppa;
                        elseif strcmp(prefix ,'AEdist')
                            psychopy = aedist; clear aedistl
                        else
                            error(['neznamy typ testu ' prefix']);
                        end
                        load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).rjepoch]);
                        if ~exist('mults','var'),  mults = []; end
                        if ~exist('header','var'), header = []; end
                        if ERP
                            E = CiEEGData(d,tabs,fs,mults,header);
                            E.GetHHeader(H);
                            E.Filter([0 60],[],[],0); %odfiltruju vsechno nad 60Hz, nekreslim obrazek
                            
                            E.Decimate(4); % ze 512 Hz na 128Hz. To staci na 60Hz signal
                            classname = 'CiEEG';
                        else
                            E = CHilbert(d,tabs,fs,mults,header);
                            E.GetHHeader(H);
                            classname = 'CHilbert';
                        end
                        clear d;                        
                        E.RejectChannels(pacienti(p).rjch);
                        if numel(reference(r).char)>0 %pokud se ma zmenit reference
                            E.ChangeReference(reference(r).char);
                        end
                        if ~ERP 
                            E.PasmoFrekvence(frekvence(f).freq);
                        end
                        disp('extracting epochs ...');
                        E.ExtractEpochs(psychopy,timewindow,baseline);        
                        E.RejectEpochs(RjEpoch);
                        E.ResponseSearch(0.1,stat_kats, stat_opak); %statistika s klouzavym oknem 100ms
                        disp('saving data ...');
                        
                        E.Save([ basedir pacienti(p).folder '\' subfolder '\' prefix ' ' classname ' ' frekvence(f).freqname ' ' reference(r).name ' ' suffix]);
                        disp([ pacienti(p).folder ' OK']);
                        clear E d tabs fs mults header RjEpoch psychopy H ans; 
                    end
                end
            end
        end
    end
end

if hybernovat
    system('shutdown -h') %#ok<UNRCH>
elseif vypnout
    system('shutdown -s') %#ok<UNRCH>
end
