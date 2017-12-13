%15.9.2016 - AlloEgo zarovnani podle odpovedi
%25.5.2017 - Pridani reference a ERP
hybernovat = 1; %jestli chci po konci skriptu pocitac uspat - ma prednost
vypnout = 0; %jestli chci po konci skriptu pocitac vypnout (a nechci ho hybernovat) 
pouzetest = 0; %jestli chci jen otestovat pritomnost vsech souboru 
overwrite = 0; %jestil se maji prepsat puvodni data, nebo ohlasit chyba a pokracovat v dalsim souboru 

setup = setup_aedist( 0 ); %nacte nastaveni testu - 11.12.2017
basedir = setup.basedir;
epochtime = setup.epochtime;
baseline = setup.baseline;
suffix = setup.suffix;
prefix = setup.prefix;
stat_kats = setup.stat_kats;
stat_opak = setup.stat_opak;
subfolder = setup.subfolder;

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
frekvence(f).todo = 1;
frekvence(f).freq = 50:5:120;
frekvence(f).freqname = '50-120'; %gamma 2
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
frekvence(f).freq = 2:2:150;
frekvence(f).freqname = '2-150'; % all range
frekvence(f).prekryv = 0.5; % 50% prekryv sousednich frekvencnich pasem 

refence = struct;
r=1;
reference(r).todo = 0;
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
reference(r).todo = 1;
reference(r).name = 'refBipo';
reference(r).char = 'b';

pacienti = pacienti_aedist(); %nactu celou strukturu pacientu
logfilename = ['logs\BatchHilbert_AEdist_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.log'];
[fileID,message] = fopen(logfilename,'wt'); %soubor na logovani prubehu
assert(fileID>=0,['nemohu otevrit soubor pro zapis: ' logfilename ]);

%nejdriv overim, jestli existuje vsechno co potrebuju nacist
chybasoubor = false;
for p = 1:numel(pacienti)
    if pacienti(p).todo 
        if(exist([basedir pacienti(p).folder '\' pacienti(p).data],'file')~=2)
            if(exist([basedir pacienti(p).folder '\' subfolder '\'  pacienti(p).data],'file')~=2)
               msg = ['Data neexistuji: ' pacienti(p).folder '\\' pacienti(p).data];
               disp(msg); fprintf(fileID,[msg '\n']);
               chybasoubor = true; 
            else
               datafolder = ['\' subfolder];
            end
        else
            datafolder = '' ;
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
    end; 
end
if (vypnout),    disp('system se po dokonceni vypne'); end 
if (hybernovat), disp('system se po dokonceni uspi'); end
clear E d tabs fs mults header RjEpoch psychopy H ans; %vymazu, kdyby tam byl nejaky zbytek z predchozich pacientu

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
                            if ERP
                                classname = 'CiEEG';
                                suffixclass = '.mat';
                            else
                                classname = 'CHilbert';
                                suffixclass = '_CHilb.mat';
                            end
                            outfilename = [ basedir pacienti(p).folder '\' subfolder '\' prefix ' ' classname ' ' frekvence(f).freqname ' ' reference(r).name ' ' suffix];
                            if exist([outfilename suffixclass],'file')==2 && overwrite == 0                                
                                disp([ outfilename ' NEULOZENO, preskoceno']); 
                                fprintf(fileID,[ 'NEULOZENO,preskoceno: ' strrep(outfilename,'\','\\') ' - ' datestr(now) '\n']); 
                                continue; %dalsi polozka ve for cyklu     
                            elseif overwrite == 0
                                disp(['soubor zatim neexistuje - zpracovavam: ' outfilename suffixclass]); 
                            end
                            load([basedir pacienti(p).folder datafolder '\' pacienti(p).data]);
                            load([basedir pacienti(p).folder '\' pacienti(p).header]);
                            load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy]);
                            if strcmp(prefix,'PPA')
                                psychopy = ppa; clear ppa;
                            elseif strcmp(prefix ,'AEdist')
                                psychopy = aedist; clear aedistl
                            elseif strcmp(prefix ,'Menrot')
                                psychopy = menrot; clear menrot;
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
                            else
                                E = CHilbert(d,tabs,fs,mults,header);
                                E.GetHHeader(H);                                
                            end
                            clear d;                        
                            E.RejectChannels(pacienti(p).rjch);
                            if numel(reference(r).char)>0 %pokud se ma zmenit reference
                                E.ChangeReference(reference(r).char);
                            end
                            if ~ERP
                                if isfield(frekvence(f),'prekryv') && ~isempty(frekvence(f).prekryv)
                                    prekryv = frekvence(f).prekryv;
                                else
                                    prekryv = 0;  %defaultne je nulovy prekryv pasem                                    
                                end
                                E.PasmoFrekvence(frekvence(f).freq,[],prekryv); %1
                            end
                            disp('extracting epochs ...');
                            E.ExtractEpochs(psychopy,epochtime,baseline);        
                            if exist('RjEpoch','var') %muze byt prazne, pak se nevyrazuji zadne epochy
                                E.RejectEpochs(RjEpoch); %globalne vyrazene epochy
                            end
                            if exist('RjEpochCh','var')
                                E.RejectEpochs(0,RjEpochCh); %epochy pro kazdy kanal zvlast
                            end
                            %E.ResampleEpochs(); % 27.11.2017 %resampluju na -1 1s podle casu odpovedi
                            %E.Decimate(8); %ze 512 na 64hz, protoze jsem predtim v PasmoFrekvence nedecimoval
                            E.ResponseSearch(0.1,stat_kats, stat_opak); %statistika s klouzavym oknem 100ms
                            disp('saving data ...');
                                                        
                            E.Save(outfilename);                            
                            disp([ pacienti(p).folder ' OK']); 
                            fprintf(fileID,[ 'OK: ' strrep(outfilename,'\','\\') ' - ' datestr(now) '\n']);
                            
                            clear E d tabs fs mults header RjEpoch psychopy H ans; 
                        catch exception 
                            errorMessage = sprintf('** Error in function %s() at line %d.\nError Message:\n%s', ...
                                exception.stack(1).name, exception.stack(1).line, exception.message);                            
                            disp(errorMessage);  fprintf(fileID,[errorMessage '\n']);  %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal                            
                            clear E d tabs fs mults header RjEpoch psychopy H ans; 
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