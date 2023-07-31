function filenames = BatchHilbert(testname,cfg)
%BATCHHILBERT prevedeno do funkce pro opakovane volani ve skriptu - 30.1.2018
% cfg je konfigurace 
%15.9.2016 - AlloEgo zarovnani podle odpovedi
%25.5.2017 - Pridani reference a ERP

%default configuration - values in cfg.
default = struct;
default.decimatefactor = 8;
default.extractepochs = 1;

if ~exist('cfg','var'), cfg = struct; end %pokud zadnou strukturu neuvedu, pouzivaji se defaultni nastaveni
if ~isfield(cfg,'hybernovat'), cfg.hybernovat = 0; end %jestli chci po konci skriptu pocitac uspat - ma prednost
if ~isfield(cfg,'vypnout'), cfg.vypnout = 0; end %jestli chci po konci skriptu pocitac vypnout (a nechci ho hybernovat) 
if ~isfield(cfg,'pouzetest'), cfg.pouzetest = 0; end %jestli chci jen otestovat pritomnost vsech souboru 
if ~isfield(cfg,'overwrite'), cfg.overwrite = 0;  end %jestil se maji prepsat puvodni data, nebo ohlasit chyba a pokracovat v dalsim souboru 
if ~isfield(cfg,'podilcasuodpovedi'), cfg.podilcasuodpovedi = 0; end  %jestli se maji epochy resamplovat na podil casu mezi podnetem a odpovedi
if ~isfield(cfg,'freqepochs'), cfg.freqepochs = 0; end %jestli se maji uklada frekvencni data od vsech epoch - velka data!
if ~isfield(cfg,'extractepochs'), cfg.extractepochs = default.extractepochs; end %muzu uklada nezepochovana data
if ~isfield(cfg,'typeEpochs'), cfg.typeEpochs = 0; end %jestli se maji epochy zarovnava podle odpovedi
if ~isfield(cfg,'suffix'), cfg.suffix = ['Ep' datestr(now,'YYYY-mm')]; end %defaultne automaticka pripona rok-mesic
if ~isfield(cfg,'pacienti'), cfg.pacienti = {}; end %muzu analyzovat jen vyber pacientu
if ~isfield(cfg,'normalization'), cfg.normalization = 'orig'; end %type of normalization after hilbert transform
if ~isfield(cfg,'statmethod'), cfg.statmethod = struct('test','wilcox','chn',1,'fdr',1); end %for explanation see BatchStat
if ~isfield(cfg,'epochfilter'), cfg.epochfilter = []; end % I can select only some epochs and include them in the resulting file, e.g. {8,[1 2]}; - epochs with values 1 or 2 in the 8th column of obj.P.data
if ~isfield(cfg,'decimatefactor'), cfg.decimatefactor = default.decimatefactor; end %decimate factor to use after hilbert tranform. 8 should be enough for trial-average analysis
if ~isfield(cfg,'normalizeEpochs'), cfg.normalizeEpochs = 1; end % if to normalize epochs by substracting the baseline - 20230/ because of memact

[ pacienti, setup,frekvence,reference  ] = pacienti_setup_load( testname,cfg.typeEpochs ); %11.1.2018 - 0 = zarovnani podle podnetu, 1=zarovnani podle odpovedi; 2023 - in memact, epoch type 0-3: immediate, epochs before delay, after delay, and within delay
if numel(cfg.pacienti)>0
    pacienti = filterpac(pacienti,cfg.pacienti);
end
basedir = setup.basedir;
epochtime = setup.epochtime;
if cfg.normalizeEpochs == 0
    baseline = [0 0]; % if we don't want to normalize epochs (in memact test, before connecting two parts of epochs together) 
else
    baseline = setup.baseline;
end
if isempty(cfg.epochfilter) && isfield(setup, 'filter') && ~isempty(setup.filter)
    cfg.epochfilter = setup.filter; % if we didn't specify filter in batchmemact.m, we'll use it from the test setup
end
suffix = cfg.suffix;  % napriklad 'Ep2018-01' + Resp pokud serazeno podle odpovedi
if cfg.podilcasuodpovedi == 1, suffix = [suffix 'PCO']; end %pokud, pridam jeste na konec priponu
% if cfg.typeEpochs > 0,  suffix = [suffix 'RES']; end %pokud zarovnavam podle odpovedi, pridavam priponu
if cfg.typeEpochs > 0, suffix = [suffix setup.suffix]; end  % 2023/07 Sifiia : to distinguish epoch types, eg. delayed in memact
if cfg.freqepochs == 1, suffix = [suffix ' FE']; end %pokud, pridam jeste na konec priponu
if cfg.extractepochs == 0, suffix = [suffix ' noEp']; end %pokud, pridam jeste na konec priponu
if cfg.decimatefactor ~= default.decimatefactor, suffix = [suffix ' D' num2str(cfg.decimatefactor)]; end %pokud, pridam jeste na konec priponu

prefix = setup.prefix;
stat_kats = setup.stat_kats;
stat_opak = setup.stat_opak;
subfolder = setup.subfolder;

logfilename = ['logs\BatchHilbert_' setup.prefix '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.log'];
[fileID,~] = fopen(logfilename,'wt'); %soubor na logovani prubehu
assert(fileID>=0,['nemohu otevrit soubor pro zapis: ' logfilename ]);
setuptext = setup2text(setup,cfg);
fprintf(fileID,setuptext); %ulozi setup do log souboru
filenames = cell(0,0); %tam si budu ukladat jmena vystupnich souboru, abych je mohl primo pouzit v BatchExtracts

%OVERENI - nejdriv overim, jestli existuje vsechno co potrebuju nacist 
chybasoubor = false;
for p = 1:numel(pacienti)
    if pacienti(p).todo        
        if(exist([basedir pacienti(p).folder '\' pacienti(p).data],'file')~=2)
            if(exist([basedir pacienti(p).folder '\' subfolder '\'  pacienti(p).data],'file')~=2)
                msg = ['Data neexistuji: ' pacienti(p).folder '\' pacienti(p).data];
                disp(msg); fprintf(fileID,[msg '\n']);
                chybasoubor = true; 
            else
                datafolder = ['\' subfolder];
            end
        else
            datafolder = '' ;
            fprintf(fileID,[ 'OK: ' pacienti(p).folder '\\' pacienti(p).data  '\n']);
        end
        if(exist([basedir pacienti(p).folder '\' pacienti(p).header],'file')~=2)
            msg = ['Header neexistuje: ' pacienti(p).folder '\\' pacienti(p).header];
            disp(msg); fprintf(fileID,[msg '\n']);
            chybasoubor = true;
        else
            fprintf(fileID,[ 'OK: ' pacienti(p).folder '\\' pacienti(p).header  '\n']);
        end
        if(exist([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy],'file')~=2)
            msg = ['Psychopy soubor neexistuje: ' pacienti(p).folder '\\' subfolder '\\' pacienti(p).psychopy]; 
            disp(msg); fprintf(fileID,[msg '\n']);
            chybasoubor = true;  
        else
            fprintf(fileID,[ 'OK: ' pacienti(p).folder '\\' pacienti(p).psychopy  '\n']);
        end
        if ~isempty(pacienti(p).rjepoch)  %muze byt prazne, pak se nevyrazuji zadne epochy     
            if(exist([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).rjepoch],'file')~=2)
                msg = ['rjepoch neexistuje: ' pacienti(p).folder '\\' subfolder '\\' pacienti(p).rjepoch]; 
                disp(msg); fprintf(fileID,[msg '\n']);
                chybasoubor = true;                 
            else
                fprintf(fileID,[ 'OK: ' pacienti(p).folder '\\' pacienti(p).rjepoch  '\n']);
            end         
        end
    end
end
if chybasoubor 
   fclose(fileID);
   error('nenalezeny nektere soubory');
else
    msg  = 'vsechny soubory ok';
    disp(msg);  fprintf(fileID,[msg '\n']); 
    if cfg.pouzetest
        disp('pouze test existence souboru');
        fclose(fileID);
        return; 
    end 
end

if (cfg.vypnout),    disp('system se po dokonceni vypne'); end 
if (cfg.hybernovat), disp('system se po dokonceni uspi'); end
clear E d tabs fs mults header RjEpoch psychopy H ans; %vymazu, kdyby tam byl nejaky zbytek z predchozich pacientu

pocetcyklu = sum([frekvence.todo]) * sum([reference.todo]) * sum([pacienti.todo]);
filestodo = pocetcyklu; %pocet souboru k vyhodnoceni, kvuli odhadu casu
souborystats = zeros(1,3); %statistika souboru - vynechane, ulozene, chybne
tablelog = cell(pocetcyklu+1,7); %frekvence, soubor, reference, status, chyba -  z toho bude vystupni xls tabulka s prehledem vysledku
tablelog(1,:) = {'freq','folder','cisloPac','ref','result','file','datetime'}; %hlavicky xls tabulky
cyklus = 1; fileno = 1;
batchtimer = tic;
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
                            elseif isfield(frekvence(f),'classname') && strcmp(frekvence(f).classname,'Morlet')
                                classname = 'CMorlet';
                                suffixclass = '_CHilb.mat';
                            else
                                classname = 'CHilbert';
                                suffixclass = '_CHilb.mat'; 
                            end
                            outfilename = [ basedir pacienti(p).folder '\' subfolder '\' prefix ' ' classname ' ' frekvence(f).freqname ' ' sprintf('%.1f-%.1f',epochtime(1:2)) ' ' reference(r).name ' ' suffix];
                            if exist([outfilename suffixclass],'file')==2 && cfg.overwrite == 0                                
                                disp([ outfilename ' NEULOZENO, preskoceno']); 
                                fprintf(fileID,[ 'NEULOZENO,preskoceno: ' strrep(outfilename,'\','\\') ' - ' datestr(now) '\n']); 
                                souborystats(1) = souborystats(1) + 1; %dalsi preskoceny soubor
                                tablelog(cyklus+1,:) = {['''' frekvence(f).freqname], pacienti(p).folder, num2str(p), reference(r).name, 'preskoceno','',datestr(now) };
                                cyklus = cyklus + 1;
                                filestodo = filestodo -1; %preskocene soubor nepocitam do celkoveho poctu
                                continue; %dalsi polozka ve for cyklu     
                            elseif cfg.overwrite == 0
                                disp(['soubor zatim neexistuje - zpracovavam: ' outfilename suffixclass]); 
                            end
                            load([basedir pacienti(p).folder datafolder '\' pacienti(p).data]); %#ok<LOAD>
                            load([basedir pacienti(p).folder '\' pacienti(p).header]); %#ok<LOAD>
                            load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).psychopy]); %#ok<LOAD>
                            if strcmp(prefix,'PPA')
                                psychopy = ppa; clear ppa;
                            elseif strcmp(prefix ,'AEdist')
                                psychopy = aedist; clear aedist;
                            elseif strcmp(prefix ,'Menrot')
                                psychopy = menrot; clear menrot;
                            elseif strcmp(prefix ,'MemAct')
                                psychopy = memact; clear memact;
                            else
                                msg = ['BatchHilbert: unknown test name ' prefix];
                                fprintf(fileID,[msg '\n']);
                                error(msg);
                            end
                            if ~isempty(pacienti(p).rjepoch)                         
                                load([basedir pacienti(p).folder '\' subfolder '\' pacienti(p).rjepoch]); %#ok<LOAD>
                            end
                            if ~exist('mults','var'),  mults = []; end
                            if ~exist('header','var'), header = []; end
                            if ERP
                                E = CiEEGData(d,tabs,fs,mults,header);                              
                            elseif strcmp(classname,'CMorlet')
                                E = CMorlet(d,tabs,fs,mults,header);                                
                            else
                                E = CHilbert(d,tabs,fs,mults,header);                                                                
                            end
                            E.GetHHeader(H);
                            if E.fs == 500 %pokud se jedna o wifi data s frekvenci 500Hz
                                E.Resample(512); %prevzorkuju na 512 Hz
                            elseif E.fs ~= 512
                                error(['soubor ma nekompatibilni vzorkovaci frekvenci ' num2str(obj.fs)]);
                            end
                            if ERP
                                E.Filter([0 60],[],[],0); %odfiltruju vsechno nad 60Hz, nekreslim obrazek
                                E.Decimate(4); % ze 512 Hz na 128Hz. To staci na 60Hz signal
                            else
                                E.Filter({0.5 [50 100 150] },[],[],0); %notch filter for [50 100 150]+-.5Hz
                            end
                            clear d;                        
                            E.RejectChannels(pacienti(p).rjch);
                            epieventfile = [basedir pacienti(p).folder '\' subfolder '\' pacienti(p).epievents];
                            if exist(epieventfile,'file')==2 %pokud existuji, nactu epieventy
                                 load(epieventfile); %#ok<LOAD>
                                 E.GetEpiEvents(DE); 
                            else
                                disp(['epievent soubor neexistuje: ' epieventfile]);
                            end
                            if numel(reference(r).char)>0 %pokud se ma zmenit reference
                                E.ChangeReference(reference(r).char);
                            end
                            if ~ERP
                                if isfield(frekvence(f),'prekryv') && ~isempty(frekvence(f).prekryv)
                                    prekryv = frekvence(f).prekryv;
                                else
                                    prekryv = 0;  %defaultne je nulovy prekryv pasem                                    
                                end                                
                                E.PasmoFrekvence(frekvence(f).freq,[],prekryv,cfg.decimatefactor); % nemuzu  najit puvod: iff(cfg.podilcasuodpovedi,2,[]),
                                    %pokud podilcasu, zdecimuju zatim jen malo, cele se mi ale nevejde do pameti
                                E.Normalize(cfg.normalization); %normalize the frequency bands
                            end                            
                            if cfg.extractepochs 
                                disp('extracting epochs ...');
                                if ERP
                                    E.ExtractEpochs(psychopy,epochtime,baseline,cfg.epochfilter);                                 
                                else
                                    E.ExtractEpochs(psychopy,epochtime,baseline,cfg.freqepochs,cfg.epochfilter);   
                                end
                                if exist('rjepoch','var') && isstruct(rjepoch) % 2023 Sofiia: in case of memact test, rjepoch is struct containing RjEpoch and RjEpochCh for each epoch type
                                    E.RejectEpochs(rjepoch(setup.index).RjEpoch);
                                    E.RejectEpochs(0,rjepoch(setup.index).RjEpochCh);
                                else % for other tests, use the old method
                                    if exist('RjEpoch','var') %muze byt prazne, pak se nevyrazuji zadne epochy
                                        E.RejectEpochs(RjEpoch); %globalne vyrazene epochy
                                    end
                                    if exist('RjEpochCh','var')
                                        E.RejectEpochs(0,RjEpochCh); %epochy pro kazdy kanal zvlast
                                    end
                                end
                                if cfg.podilcasuodpovedi == 1                            
                                    E.ResampleEpochs(); % 27.11.2017 %resampluju na -1 1s podle casu odpovedi
                                    E.Decimate(4); %ze 256 na 64hz, protoze jsem predtim v PasmoFrekvence decimoval jen 2x
                                end
                                %vypocet statistiky
                                E.ResponseSearchMulti(0.1,stat_kats);
%                                 if iscell(stat_kats) && numel(stat_kats)>1 && iscelldeep(stat_kats)  %pokud mam nekolik ruznych kontrastu na spocitani
%                                     for WpA = 1:numel(stat_kats)
%                                         E.SetStatActive(WpA);
%                                         disp(['pocitam kontrast ' num2str(WpA) ': ' cell2str(stat_kats{WpA}) ]);
%                                         E.ResponseSearch(0.1,stat_kats{WpA},stat_opak);
%                                     end
%                                 else  %jen jeden kontrast
%                                     E.ResponseSearch(0.1,stat_kats, stat_opak); %statistika s klouzavym oknem 100ms
%                                 end  
                            end
                            disp('saving data ...');                            
                            E.Save(outfilename);                            
                            disp([ pacienti(p).folder ' OK']); 
                            fprintf(fileID,[ 'OK: ' strrep(outfilename,'\','\\') ' - ' datestr(now) '\n']);
                            souborystats(2) = souborystats(2) + 1; %dalsi ulozeny soubor
                            tablelog(cyklus+1,:) = {['''' frekvence(f).freqname], pacienti(p).folder, num2str(p), reference(r).name, 'saved', outfilename,datestr(now) };                                
                            if isempty(find(~cellfun('isempty',strfind(filenames,basename(outfilename))))), filenames{end+1,1} = basename(outfilename); end %#ok<AGROW,EFIND>
                            clear E d tabs fs mults header RjEpoch RjEpochCh psychopy H ans; 
                        catch exception 
                            errorMessage = exceptionLog(exception);                         
                            disp(errorMessage);  fprintf(fileID,[errorMessage '\n']); %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal                            
                            souborystats(3) = souborystats(3) + 1; %dalsi chybny soubor
                            tablelog(cyklus+1,:) = {['''' frekvence(f).freqname], pacienti(p).folder, num2str(p), reference(r).name, 'error', exception.message , datestr(now)}; 
                            clear E d tabs fs mults header RjEpoch psychopy H ans; 
                        end    
                        cas = toc(batchtimer);
                        odhadcelehocasu = filestodo/fileno * cas;
                        fprintf(' %i/%i : cas zatim: %.1f min, zbyvajici cas %.1f min\n',cyklus,pocetcyklu,cas/60,(odhadcelehocasu - cas)/60); %vypisu v kolikatem jsem cyklu a kolik zbyva sekund do konce
                        cyklus = cyklus + 1; %cyklus pres vsechny i preskocene soubory
                        fileno = fileno + 1; %cyklus pres vsechny nepreskocene soubory
                        xlswrite([logfilename '.xls'],tablelog); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
                    end                    
                end
            end
        end
    end
end


stroutput = sprintf('preskoceno: %i, ulozeno: %i, chybnych: %i\n',souborystats(1),souborystats(2),souborystats(3));
fprintf(stroutput);  fprintf(fileID,stroutput); 
fclose(fileID);
if cfg.hybernovat
    system('shutdown -h') 
elseif cfg.vypnout            
    system('shutdown -s')
end

end  %function
function pacienti= filterpac(pacienti,filter)
    pacremove = []; %seznam pacientu k vyrazeni
    for p = 1 : numel(pacienti)
        nalezen = false;
        for f = 1:numel(filter)
            if strfind(pacienti(p).folder,filter{f})
                nalezen = true; %pacient je uveden ve filtru
                break;
            end
        end
        if ~nalezen
            pacremove = [pacremove p]; %#ok<AGROW> %for cyklus porad bere z puvodniho array, takze ho nemuzu zmensovat v ramci for cyklu
        end
    end
    pacienti(pacremove) = [];
end