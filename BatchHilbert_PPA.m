%15.9.2016 - AlloEgo zarovnani podle odpovedi
%25.5.2017 - Pridani reference a ERP
hybernovat = 1; %jestli chci po konci skriptu pocitac uspat - ma prednost
vypnout = 0; %jestli chci po konci skriptu pocitac vypnout (a nechci ho hybernovat) 
pouzetest = 0; %jestli chci jen otestovat pritomnost vsech souboru 

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
frekvence(f).todo = 0;
frekvence(f).freq = [];
frekvence(f).freqname = 'ERP'; % slow gamma
f=2;
frekvence(f).todo = 0;
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
frekvence(f).todo = 1;
frekvence(f).freq = 4:2:150;
frekvence(f).freqname = '4-150'; % all range

refence = struct;
r=1;
reference(r).todo = 1;
reference(r).name = 'refOrig';
reference(r).char = '';
r=2;
reference(r).todo = 0;
reference(r).name = 'refEle';
reference(r).char = 'e';
r=3;
reference(r).todo = 0;
reference(r).name = 'refHead';
reference(r).char = 'h';
r=4;
reference(r).todo = 0;
reference(r).name = 'refBipo';
reference(r).char = 'b';

pacienti = struct;
p = 1;
pacienti(p).todo = 0;
pacienti(p).folder = 'p073 Pech VT6'; %'p073_VT6';
pacienti(p).data = 'VT6_INV Test Vlcek1_X_ppa.mat';
pacienti(p).header = 'p73_header_kamil.mat';
pacienti(p).psychopy = 'p73_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p = 2;
pacienti(p).todo = 0;
pacienti(p).folder = 'p079 Plu VT8'; %'p79_VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_ppa.mat';
pacienti(p).header = 'P79_header.mat';
pacienti(p).psychopy = 'p79_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p79.mat'; 
pacienti(p).rjch = [47 64 68]; %#ok<NBR%#ok<MSNU> AK> puvodne 47 64 114
pacienti(p).frekvence = struct;

p=3;
pacienti(p).todo = 0;
pacienti(p).folder = 'p082 Vov VT9'; %'VT09';
pacienti(p).data = 'VT9_2015-04-21_09-46_001_concat_X_ppa.mat';
pacienti(p).header = 'P82_header.mat';
pacienti(p).psychopy = 'p82_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = [47 68 126]; %#ok<NB%#ok<MSNU> RAK>

p=4;
pacienti(p).todo = 0;
pacienti(p).folder = 'p083 Kol VT10'; %'VT10';
pacienti(p).data = 'VT10_2015-05-19_10-00_001_X_ppa.mat';
pacienti(p).header = 'P83_header.mat';
pacienti(p).psychopy = 'p83_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p83.mat';
pacienti(p).rjch = [47 64]; %  64 ?

p=5;
pacienti(p).todo = 0;
pacienti(p).folder = 'p095 Hav VT11'; %'VT11';
pacienti(p).data = 'VT11_2015-12-15_ppa.mat';
pacienti(p).header = 'P95_header.mat';
pacienti(p).psychopy = 'p95_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p95.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK> 64?

p=6;
pacienti(p).todo = 0;
pacienti(p).folder = 'p096 Gro VT12'; %'VT12';
pacienti(p).data = 'VT12_2016-01-26_09-16_001_concat_ppa.mat';
pacienti(p).header = 'P96_header.mat';
pacienti(p).psychopy = 'p96_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p96.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK> 64?

p=7;
pacienti(p).todo = 0;
pacienti(p).folder = 'p097 Nov VT13'; %'VT13';
pacienti(p).data = 'VT13_2016-02-11_09-20_001_concat_ppa.mat';
pacienti(p).header = 'P97_header.mat';
pacienti(p).psychopy = 'p97_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p97.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=8;
pacienti(p).todo = 0;
pacienti(p).folder = 'p110 Sou VT14'; %'VT14';
pacienti(p).data = 'P110_2016-06-08_15-56_001_concat_ppa.mat';
pacienti(p).header = 'p110_header.mat';
pacienti(p).psychopy = 'p110_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = [47 15]; %#ok<NB%#ok<MSNU> RAK>

p=9;% moc zachvatu?
pacienti(p).todo = 0;
pacienti(p).folder = 'p126 Sve VT15'; %'VT15';
pacienti(p).data = 'VT15_2016-09-06_09-18_001_concat_ppa.mat';
pacienti(p).header = 'P126_header.mat';
pacienti(p).psychopy = 'p126_ppa.mat';
pacienti(p).rjepoch = ''; %muze byt prazne, pak se nevyrazuji zadne epochy
pacienti(p).rjch = [47 50]; %#ok<NB%#ok<MSNU> RAK>

p=10;
pacienti(p).todo = 0;
pacienti(p).folder = 'p119 Buc VT16'; %'VT16';
pacienti(p).data = 'VT16_2016-10-10_17-11_001_concat_ppa.mat';
pacienti(p).header = 'p119_header_kamil.mat';%chybi header Hammer
pacienti(p).psychopy = 'p119_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = [47 57 64]; %#ok<NB%#ok<MSNU> RAK>

p= 11;
pacienti(p).todo = 1;
pacienti(p).folder = 'p130 Koc VT17'; %'VT17';
pacienti(p).data = 'VT17_2016-10-24_15-41_001_concat_ppa.mat';
pacienti(p).header = 'p130_header.mat';
pacienti(p).psychopy = 'p130_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = [47]; %#ok<NB%#ok<MSNU> RAK>

p= 12;
pacienti(p).todo = 1;
pacienti(p).folder = 'p132 Pol VT18'; %'VT18';
pacienti(p).data = 'VT18_2016-12-07_16-20_001_ppa.mat';
pacienti(p).header = 'P132_header_kamil.mat';% lisi se
pacienti(p).psychopy = 'p132_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = []; %#ok<NB%#ok<MSNU> RAK>

p= 13;
pacienti(p).todo = 1;
pacienti(p).folder = 'p129 Kuch VT19'; %'VT19';
pacienti(p).data = 'VT19_2017-01-16_09-38_001_concat_ppa.mat';
pacienti(p).header = 'p129_header.mat';
pacienti(p).psychopy = 'p129_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = [47];  %POZOR na elektrody 34 a 35 od epochy 166-169 + elektroda 67 epocha 182 

p= 14;
pacienti(p).todo = 1;
pacienti(p).folder = 'p136 Men VT20'; %'VT20';
pacienti(p).data = 'VT20_2017-01-30_09-40_001_concat_ppa.mat';
pacienti(p).header = 'p136_header.mat'; %'P83_64ch_header.mat'
pacienti(p).psychopy = 'p136_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = [47]; %#ok<NBRAK>

p= 15;
pacienti(p).todo = 1;
pacienti(p).folder = 'p138 Ven VT21'; %'VT21';
pacienti(p).data = 'VT21_2017-02-28_09-35_001_500hz_concat_ppa.mat';
pacienti(p).header = 'p138_header_kamil.mat'; %'P83_64ch_header.mat'
pacienti(p).psychopy = 'p138_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).rjch = []; 

[fileID,message] = fopen(['logs\BatchHilbert_PPA_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.log'],'wt'); %soubor na logovani prubehu
%nejdriv overim, jestli existuje vsechno co potrebuju nacist
chybasoubor = false;
for p = 1:numel(pacienti)
    if pacienti(p).todo 
        if(exist([basedir pacienti(p).folder '\' pacienti(p).data],'file')~=2)
            msg = ['Data neexistuji: ' pacienti(p).folder '\\' pacienti(p).data];
            disp(msg); fprintf(fileID,[msg '\n']);
            chybasoubor = true; 
        else
            fprintf(fileID,[ 'OK: ' pacienti(p).folder '\\' pacienti(p).data  '\n']);
        end;
        if(exist([basedir pacienti(p).folder '\' pacienti(p).header],'file')~=2)
            msg = ['Header neexistuje: ' pacienti(p).folder '\\' pacienti(p).header];
            disp(msg); fprintf(fileID,[msg '\n']);
            chybasoubor = true;
        else
            fprintf(fileID,[ 'OK: ' pacienti(p).folder '\\' pacienti(p).header  '\n']);
        end;
        if(exist([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy],'file')~=2)
            msg = ['Psychopy soubor neexistuje: ' pacienti(p).folder '\\' subfolder '\\' pacienti(p).psychopy]; 
            disp(msg); fprintf(fileID,[msg '\n']);
            chybasoubor = true;  
        else
            fprintf(fileID,[ 'OK: ' pacienti(p).folder '\\' pacienti(p).psychopy  '\n']);
        end;
        if ~isempty(pacienti(p).rjepoch)  %muze byt prazne, pak se nevyrazuji zadne epochy     
            if(exist([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).rjepoch],'file')~=2)
                msg = ['rjepoch neexistuje: ' pacienti(p).folder '\\' subfolder '\\' pacienti(p).rjepoch]; 
                disp(msg); fprintf(fileID,[msg '\n']);
                chybasoubor = true;                 
            else
                fprintf(fileID,[ 'OK: ' pacienti(p).folder '\\' pacienti(p).rjepoch  '\n']);
            end;            
        end
    end
end
if chybasoubor 
   fclose(fileID);
   error('nenalezeny nektere soubory');
else
    msg  = 'vsechny soubory ok';
    disp(msg);  fprintf(fileID,[msg '\n']); 
    if pouzetest
        disp('pouze test existence souboru');
        fclose(fileID);
        return; 
    end; %#ok<UNRCH>
end
if (vypnout),    disp('system se po dokonceni vypne'); end %#ok<UNRCH>
if (hybernovat), disp('system se po dokonceni uspi'); end

for f=1:numel(frekvence)        
    if frekvence(f).todo
        msg  =[ '*****' frekvence(f).freqname '*****' ]; 
        disp(msg);  fprintf(fileID,[msg '\n']); 
        if numel(frekvence(f).freq) == 0
            ERP = 1; 
        else
            ERP = 0; 
        end
        for p = 1:numel(pacienti)            
            if pacienti(p).todo
                msg  = [ ' ---- ' pacienti(p).folder ' ---- '];
                disp(msg);  fprintf(fileID,[msg '\n']); 
                for r = 1:numel(reference)
                    if reference(r).todo
                        msg  = [ ' ..... ' reference(r).name ' ..... ' datestr(now)]; %datum a cas
                        disp(msg);  fprintf(fileID,[msg '\n']); 
                        
                        try %i kdy bude nejaka chyba pri vyhodnoceni, chci pokracovat dalsimi soubory
                            load([basedir pacienti(p).folder '\' pacienti(p).data]);
                            load([basedir pacienti(p).folder '\' pacienti(p).header]);
                            load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy]);
                            if strcmp(prefix,'PPA')
                                psychopy = ppa; clear ppa;
                            elseif strcmp(prefix ,'AEdist')
                                psychopy = aedist; clear aedistl
                            else
                                msg = ['neznamy typ testu ' prefix'];
                                fprintf(fileID,[msg '\n']);
                                error(msg);
                            end
                            if ~isempty(pacienti(p).rjepoch)                         
                                load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).rjepoch]);
                            end
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
                            if ~isempty(pacienti(p).rjepoch) %muze byt prazne, pak se nevyrazuji zadne epochy
                                E.RejectEpochs(RjEpoch);
                            end
                            E.ResponseSearch(0.1,stat_kats, stat_opak); %statistika s klouzavym oknem 100ms
                            disp('saving data ...');
                            
                            outfilename = [ basedir pacienti(p).folder '\' subfolder '\' prefix ' ' classname ' ' frekvence(f).freqname ' ' reference(r).name ' ' suffix];
                            E.Save(outfilename);                            
                            disp([ pacienti(p).folder ' OK']); 
                            fprintf(fileID,[ 'OK: ' strrep(outfilename,'\','\\') ' - ' datestr(now) '\n']);
                            clear E d tabs fs mults header RjEpoch psychopy H ans; 
                        catch exception 
                            msg = ['Exception message: ' exception.message];
                            disp(msg);  fprintf(fileID,[msg '\n']);  %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal
                        end    
                    end
                end
            end
        end
    end
end
fclose(fileID);
if hybernovat
    system('shutdown -h') %#ok<UNRCH>
elseif vypnout            %#ok<UNRCH>
    system('shutdown -s') %#ok<UNRCH>
end

