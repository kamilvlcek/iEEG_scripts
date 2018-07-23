classdef CHilbertMulti < CHilbert 
    %CHILBERTMULTI zpracovava data z nekolika CHilbert vyexportovanych data
    %   nacita data vyprodukovana pomoci CHilbert.ExtractData
    
    properties
        orig; %struktura s puvodnimi daty 
        filenames; %cell array se jmeny souboru
        blokyprehazej; %matrix s originalnimi cisly bloku a cisly bloku k prehazeni z noveho souboru. Meni se pro kazdy pridany soubor
           %plni se v GetEpochData
        filesimported = 0;
        subjNames; %jmena subjektu z headeru
        label; %ze kterych extraktu byl vytvoren
        mfilename; %jmeno souboru CHilbertMulti
    end
    
    methods (Access = public)  
        function obj = CHilbertMulti(filename)              
            if exist('filename','var')
                obj.Load(filename);
            end
        end
        
        function FILES = TestExtract(obj,filenames)
            FILES = cell(size(filenames,1),2);
             for fileno = 1:size(filenames,1)
                filename = filenames{fileno,1}; %cell array, zatim to musi byt plna cesta                
                if exist(filename,'file') 
                    disp(obj.basename(filename)); %zobrazim jmeno souboru s pouze koncem 
                    obj.filenames{fileno,1} = filename;
                    clear d;
                    load(filename,'d','P','fs'); %nacte vsechny promenne
                    test = ~P.data(:,P.sloupce.zpetnavazba); %index testovych epoch
                    d = d(:,:,test);
                    %disp(['  velikost d (samples x channels x epochs):' num2str(size(d))]); 
                    FILES(fileno,:) = {obj.basename(filename), [ 'd: ' num2str(size(d),'%i ') ', fs: ' num2str(fs)]};
                else
                    disp(['soubor neexistuje ' filename]);
                    FILES(fileno,:) = {obj.basename(filename), 'soubor neexistuje'};
                end
             end
        end
        function obj = Clear(obj)
            %smaze data objektu, kvuli volani z ImportExtract, jinak je mozna jednodussi znova objekt vytvorit
            obj.d = [];
            obj.tabs = [];
            obj.Hf = [];
            obj.Hfmean  =[];
            obj.HFreq = [];
            obj.CH = [];
            obj.PsyData = {};
            obj.epochData = {};
            obj.els = [];
            obj.RjEpochCh = [];
            obj.mults = {};
            obj.epochtime = {};
            obj.baseline = {};
            
            obj.filenames = {};
            obj.filesimported = 0;         
            obj.orig = {};
            obj.blokyprehazej = {};
            obj.subjNames = {};
            obj.label = {};
            obj.mfilename = {};            
            
            obj.plotF = {};
            obj.plotRCh = {};
            obj.plotEp = {};
            obj.reference = {};
            obj.yrange = {};
            obj.DatumCas = {};
            disp('data objektu smazana');
        end
        function obj = ImportExtract(obj,filenames,label)
            if numel(obj.filenames)>0
                if obj.filesimported == 0 %nejaky soubor naimportovan castecne kvuli chybe - musim smazat
                    obj.Clear();
                else
                    disp(['pridavam data k existujicim souborum: ' num2str(obj.filesimported)]);
                end
            end
            if exist('label','var'), obj.label = label; end             
            for fileno = 1:size(filenames,1)
                filename = filenames{fileno,1}; %cell array, zatim to musi byt plna cesta                
                if exist(filename,'file') 
                    disp(obj.basename(filename)); %zobrazim jmeno souboru s pouze koncem 
                    obj.filenames{fileno,1} = filename;
                    load(filename); %nacte vsechny promenne
                    if ~exist('baseline','var'), baseline = [epochtime(1) 0]; end %fake baseline, pokud nebyla ulozena
                    test = ~P.data(:,P.sloupce.zpetnavazba); %index testovych epoch
                    
                    %identita jednotlivych epoch - epochData 
                    %predpokladam epochovana data, pro jina jsem ani nezkousel
                    obj.GetEpochData(epochData,fileno,test); %soucasne naplni obj.blokyprehazej, ktere muzu pouzivat dal                    
                    [d,tabs,RjEpochCh,P]=obj.PrehazejEpochy(d,tabs,RjEpochCh,P,test); %vyradi treningove trialy a prehazi epochy pokud je to treba
                    
                    %d hlavni data
                    obj.GetD(d); %ulozi nova data d, pripadne prehazi epochy
                    
                    %tabs a tabs_orig
                    obj.GetTabs(tabs,tabs_orig,fileno);                    
                    
                    %fs - vzorkovaci prekvence - musi byt pro vsechny soubory stejna
                    if ~isempty(obj.fs), assert(obj.fs == fs,'fs musi byt stejne');  else, obj.fs = fs; end
                                        
                    obj.GetEpochTime(epochtime,baseline);
                    
                    %vyrazene epochy x kanaly
                    obj.RjEpochCh = cat(1,obj.RjEpochCh,RjEpochCh); %spojim pres kanaly, pocet epoch musi by stejny                    
                    
                    %Hammer header
                    obj.GetHeader(H,fileno);
                    
                    %PsychoPy data
                    obj.GetPsyData(P,fileno);                                                           
                    
                    %frekvencni data
                    obj.GetHfreq(Hf,Hfmean,HFreq);
                    obj.GetRef(reference); 
                    
                    %obj.GetStat(Wp); %mam funkci, ale statistiku zatim neimportuju
                    %jen kopie - zatim nezpracovavam                                                           
                    obj.orig(fileno).DatumCas = DatumCas;
                    obj.orig(fileno).filename = filename;     
                    obj.orig(fileno).blokyprehazej = obj.blokyprehazej; %zaloha kvuli zpetne referenci
                    obj.filesimported = obj.filesimported +1;
                end
                
            end
            disp(['nacteno souboru: ' num2str(obj.filesimported)]);
        end
        function obj = SetLabel(obj,label)
            %kdyz chci pridat label jeste po tom, co jsem udelal ImportExtract
            obj.label = label;
        end
        function [d,tabs,RjEpochCh,P]= PrehazejEpochy(obj,d,tabs,RjEpochCh,P,test) 
            %vyradi treningove epochy ze vsech dat 
            %pokud je potreba, prehazi epochy podle blokyprehazej
            d = d(:,:,test); %jen testove epochy, trening vymazu
            if ~isempty(obj.d) %pokud uz existuji data v obj.d - nejedna se o prvni soubor
                assert(size(obj.d,1)==size(d,1),['pocet vzorku musi byt stejny: d:' num2str(size(d,1)) ' x predchozi:' num2str(size(obj.d,1))]);
                assert(size(obj.d,3)==size(d,3),['pocet epoch musi byt stejny: d:' num2str(size(d,3)) ' x predchozi:' num2str(size(obj.d,3))]);
                assert(size(obj.d,3)==size(obj.epochData,1), 'pocet epoch v epochData a d musi byt stejny');
            end              
            
            tabs = tabs(:,test); %velikost tabs se odviji od velikosti d, takze nemusim overovat            
            RjEpochCh = RjEpochCh(:,test);
            Psy = CPsyData(P);%vytvorim objekt CPsyData 
            Psy.RemoveTraining(); %odstranim treninkove trialy
            P = Psy.P; %zase vytahnu strukturu P ven
            
            if ~isempty(obj.blokyprehazej) %pokud je treba prehazet epochy - jine poradi bloku nez v predchozim souboru
                %tam budu ukladat kopie dat s prehazenymi epochami
                d2=zeros(size(d)); 
                tabs2 = zeros(size(tabs));
                RjEpochCh2 = zeros(size(RjEpochCh));
                Pdata2 = zeros(size(P.data));
                epocha2=1; %nove cislo epochy 1-n
                for blok2 = 1:size(obj.blokyprehazej,1) %blok 16ti epoch (u AEdist) - nova cisla bloku 1-n
                    blok1 = obj.blokyprehazej(blok2,2); %puvodni cislo bloku, ktere chci zmenit na nove                
                    for epocha1 = obj.blokyprehazej(blok1,3) :  obj.blokyprehazej(blok1,4) %jednotlive epochy - puvodni cisla epoch                 
                       
                       %vlozim epochu na novou pozici
                       d2(:,:,epocha2)= d(:,:,epocha1); %time x channels x epochs
                       tabs2(:,epocha2) = tabs(:,epocha1); %time x epochs
                       RjEpochCh2(:,epocha2) = RjEpochCh(:,epocha1); %channels x epochs
                       Pdata2(epocha2,:) = P.data(epocha1,:);
                       
                       epocha2 = epocha2 + 1;
                    end
                end
                d = d2; %clear d2;
                tabs = tabs2; %clear tabs2;
                RjEpochCh = RjEpochCh2;
                P.data = Pdata2; %clear Pdata2;
            end
            
        end
        function GetD(obj,d)
            %ulozi nova eeg data k predchazejicim - prida je do spolecneho pole obj.d                     
            obj.d = cat(2,obj.d,d); %spojim pres channels - jen data z testovych epoch            
            if isempty(obj.els)
                obj.els = size(d,2);
            else
                obj.els = [obj.els obj.els(end)+size(d,2)]; %konce elektrod, zde konce pacientu
            end
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
        end
        function GetHfreq(obj,Hf,Hfmean,HFreq)
            %spoji frekvencni data z predchozich a noveho souboru
            if isempty(obj.Hf)
                obj.Hf = Hf; %jen prvni soubor
            else
                assert(isequal(obj.Hf,Hf),'frekvencni pasma musi byt stejna');
            end
            if isempty(obj.Hfmean)
                obj.Hfmean = Hfmean; %jen prvni soubor
            end
            if isempty(obj.HFreq)
                obj.HFreq = HFreq;%time x channel x freq (x kategorie)
            else
                obj.HFreq = cat(2,obj.HFreq,HFreq);
            end
        end
        function GetRef(obj,reference)
            if isempty(obj.reference)
                obj.reference = reference; %jen prvni soubor
            else
                assert(isequal(obj.reference,reference),'reference musi byt stejna');
            end
        end
        function obj = GetTabs(obj,tabs) %,tabs_orig,f
            %spoji tabs z predchozich a noveho souboru
            if isempty(obj.tabs) %budu pouzivat tabs z prvniho souboru
                for ep = 1 : size(tabs,2) %pres vsechny epochy
                    if ep==1
                        tabs(:,1) = tabs(:,1) - tabs(1,1); %aby to zacinalo od 0
                    else
                        tabs(:,ep) = (tabs(:,ep) - tabs(1,ep)) + tabs(end,ep-1) + (tabs(end,ep-1)-tabs(end-1,ep-1));  %odectu prvni a pak prictu posledni z predchozi epochy a rozdil oproti predposlenimu
                    end
                end
                obj.tabs = tabs;                        
            end
            %ale pro kazdy soubor ulozim originalni tabs a tabs_orig   - zrusil jsem, soubory jsou hrozne velky, asi to nebudu potreboval, kdyz tak vygeneruju znova - 23.7.2018
            %obj.orig(f).tabs = tabs;                     
            %obj.orig(f).tabs_orig = tabs_orig; 
        end
        function obj = GetPsyData(obj,P,fileno)
            %pokud prvni soubor, ulozi P data
            %pokud dalsi soubor v poradi, ulozi P data z noveho souboru do zalohy v obj.ori            
            if isempty(P.pacientid) %v nekterych starych datech to neni uvedeno
                P.pacientid =  obj.subjNames{fileno}; %nahradim udajem z headeru, ktery jsem ziskal v GetHeader
            end
            obj.orig(fileno).P = P; %psychopy data  
            if isempty(obj.PsyData)
                obj.PsyData = CPsyDataMulti(P); %vytvorim objekt CPsyDataMulti - z prvniho souboru                
            else
                obj.PsyData.GetPsyData(P); %vlozim nova data pro tento subjekt a aktivuju je
                %obj.PsyData.P.pacientid = [ obj.PsyData.P.pacientid ',' P.pacientid]; %spojim pacient ID od vice pacientu
                %ostatni psydata zatim nepouzivam, reakcni uspesnost a rychlost bude u kazdeno jina
            end            
            obj.PsyData.RemoveTraining(); %odstranim treninkove trialy
        end
        function obj = GetEpochData(obj,epochData,fileno,test)
            %porovna poradi bloku v puvodnich souborej a novem souboru
            %pokud je jine, vytvori pole blokyprehazej, ktere se pak pouzije v PrehazejEpochy
            epochData =  epochData(test,:); %vyradim trening;
            if isempty(obj.epochData)                                                
                %bloky = obj.GetBlocks(epochData); 
                obj.epochData = epochData; %pouziju epochData prvniho souboru, ulozim bez treningu
                obj.blokyprehazej = []; %nebudu nic prehazovat
            else
                assert(size(obj.epochData,1)==size(epochData,1), 'pocet podminek v epochData musi byt stejny');                                                 
                bloky0 = obj.GetBlocks(obj.epochData); 
                bloky1 = obj.GetBlocks(epochData); 
                assert(size(bloky0,1)==size(bloky1,1),'pocet bloku musi byt stejny');
                prehazet = []; %v kterych blocich jsou jine podminky
                for b = 1:size(bloky0,1) %pres vsechny bloky
                   if bloky0(b,2) ~= bloky1(b,2)                               
                       assert(false,['ruzna velikost bloku ' b ]);
                   end
                   if bloky0(b,1) ~= bloky1(b,1), prehazet = [prehazet b]; end %#ok<AGROW> %budu muset poradi bloku prehazet, podminky v jinem poradi 
                end
                if numel(prehazet) > 0
                    obj.blokyprehazej = obj.ShuffleBlocks(bloky0, bloky1);                    
                else
                    obj.blokyprehazej = []; %nebudu nic prehazovat
                end
            end           
            obj.orig(fileno).epochData = epochData; % ulozim original taky                        
        end
        function obj = GetHeader(obj,H,fileno)
            %spoji header (kanaly) v puvodnich souborech a tom novem
            if isempty(obj.CH)
                %GetHHeader@CHilbert(obj,H);  
                els = obj.els;
                obj.GetHHeader(H);
                obj.els = els; %chci jina els,nez nacteno v GetHHeader
                obj.CH.chgroups = {1:numel(H.channels)};
                obj.CH.els = numel(H.channels);                
            else
                obj.CH.H.subjName = [obj.CH.H.subjName ',' H.subjName]; %spojim jmena subjektu                
                CHfields = fieldnames(obj.CH.H.channels);
                Hfields = fieldnames(H.channels); %nove pridavany header
                if numel(CHfields) == numel(Hfields)
                    obj.CH.H.channels = cat(2,obj.CH.H.channels,H.channels); %spojim udaje o kanalech
                elseif numel(CHfields) > numel(Hfields) %predchozi soubor mel vice poli v H.channels
                    rozdil = setdiff(CHfields,Hfields); %jaka pole chybi v Hfields
                    prazdnepole = cell(numel(H.channels),1); %cell array o tolika polich, ktere ma nove pridavany header
                    for r = 1:numel(rozdil)                       
                        [H.channels(:).(rozdil{r})] = prazdnepole{:}; %pridam kazde chybejici pole
                        %tomuhle kodu moc nerozumim, ale nasel jsem ho na webu
                        %(rozdil{r}) je dynamic field
                    end
                    obj.CH.H.channels = cat(2,obj.CH.H.channels,H.channels); %spojim udaje o kanalech
                else
                    error(['predchozi soubor mel mene poli v H.channels:' num2str(numel(CHfields))]);
                end
                
                obj.CH.H.selCh_H = [obj.CH.H.selCh_H  (1:numel(H.channels))+obj.CH.H.selCh_H(end) ];                
                obj.CH.chgroups{fileno} = (1:numel(H.channels)) + obj.CH.els(fileno-1);
                obj.CH.els(fileno) = numel(H.channels)+obj.CH.els(fileno-1);                
            end
            obj.subjNames{fileno,1} = H.subjName;
            obj.orig(fileno).H = H; %orignalni header ulozim, v kazdem pripade
        end
        function obj = GetEpochTime(obj,epochtime,baseline)
            %epochtime - cas epochy, napriklad -0.2 - 1.2, taky musi byt pro vsechny soubory stejne
            if ~isempty(obj.epochtime)
                assert(isequal(obj.epochtime,epochtime(1:2)),'epochtime musi byt stejne');
            else
                obj.epochtime = epochtime(1:2);
            end  
            %baseline - cas 
            if ~isempty(obj.baseline)
                assert(isequal(obj.baseline,baseline),'baseline musi byt stejne');
            else
                obj.baseline = baseline;
            end              
        end
        function obj = GetStat(obj) %zatim neimportuju
            obj.Wp = [];
        end
        function PAC = GetPAC(obj,filename)
            %vytvori strukturu, ktera se pak da nacist do CM.ExtractData pro vytvoreni stejnych extraktu z jineho souboru
            %filename je jmeno souboru, ze ktereho je tento CM soubor, napriklad 'Menrot CHilbert 50-150 -1.0-1.5 refBipo Ep2018-01_CHilb.mat'
            PAC = {}; iPAC = 1;
            previousNick = '';
            for ch = 1:numel(obj.CH.H.channels)                
                nick = obj.CH.H.channels(ch).name(1: find(obj.CH.H.channels(ch).name==' ')-1);
                if numel(nick)==3
                    nick = ['p0' nick(2:3)];
                end
                if ~strcmp(nick,previousNick) %musim nacist pacienta, E ji jine
                    E = pacient_load(nick,obj.PsyData.testname,filename);
                    previousNick = nick;
                end
                nick_ch = find([E.CH.H.channels.numberOnAmplifier]==obj.CH.H.channels(ch).numberOnAmplifier); % podle numberOnAmplifier vyhledam cislo kanalu
                PAC(iPAC).pacient = nick; %#ok<AGROW>
                PAC(iPAC).ch = nick_ch; %#ok<AGROW>
                PAC(iPAC).name = E.CH.H.channels(nick_ch).name; %#ok<AGROW>
                PAC(iPAC).neurologyLabel = E.CH.H.channels(nick_ch).neurologyLabel; %#ok<AGROW>
                PAC(iPAC).ass_brainAtlas = E.CH.H.channels(nick_ch).ass_brainAtlas;%#ok<AGROW>
                PAC(iPAC).ass_cytoarchMap = E.CH.H.channels(nick_ch).ass_cytoarchMap; %#ok<AGROW>
                iPAC = iPAC + 1;                
            end
        end
        %% SAVE AND LOAD FILE
        %dve funkce na ulozeni a nacteni dat tridy
        %uklada se vcetne dat parenta CHilbert a pres nej taky CiEEGData        
        function obj = Save(obj,filename)
            if ~exist('filename','var')
                filename = obj.mfilename;
                assert( ~isempty(filename), 'no filename given or saved before');
            else
                obj.mfilename = filename;
            end            
            Save@CHilbert(obj,CHilbert.filenameH(filename));  %ulozim do prvniho souboru data z nadrazene tridy          
            if ~isempty(obj.filenames)                
                filenames = obj.filenames;   %#ok<PROPLC,NASGU>
                orig = obj.orig;         %#ok<PROPLC,NASGU> 
                blokyprehazej = obj.blokyprehazej; %#ok<PROPLC,NASGU> 
                filesimported = obj.filesimported; %#ok<PROPLC,NASGU> %time x channel x frequency x epoch
                subjNames = obj.subjNames; %#ok<PROPLC,NASGU> 
                label = obj.label; %#ok<NASGU,PROPLC> 
                mfilename = obj.mfilename; %#ok<NASGU,PROPLC>              
                save(CHilbertMulti.filenameM(filename),'filenames','orig','blokyprehazej','filesimported','subjNames','label','mfilename','-v7.3'); %do druheho souboru data z teto tridy
            end
        end
         %pokud je druhy parametr 1, nenacitaji se data z parenta
        function obj = Load(obj,filename,onlyself)
            %parametr loadall se hodi pro FE data se vsemi ulozenymi epochami, ktere jsou giganticke
            
            if ~exist('onlyself','var') || onlyself == 0
                assert(exist(CHilbert.filenameH(filename),'file')==2, ['soubor s daty CHilbert neexistuje:' char(10) CHilbert.filenameE(filename) char(10) 'mozna se jedna o data tridy CiEEGData?']);    
                Load@CHilbert(obj,CHilbert.filenameH(filename));  
            end
            if exist(CHilbertMulti.filenameM(filename),'file')                
                load(CHilbertMulti.filenameM(filename),'orig','filenames','blokyprehazej','filesimported','subjNames','label','mfilename');
                obj.orig = orig;        %#ok<CPROPLC>
                obj.filenames = filenames;               %#ok<CPROPLC>                 
                obj.blokyprehazej = blokyprehazej; %#ok<CPROPLC> 
                obj.filesimported = filesimported; %#ok<CPROPLC> 
                obj.subjNames = subjNames; %#ok<CPROPLC> 
                obj.label = label; %#ok<CPROPLC> 
                obj.mfilename = mfilename;      %#ok<CPROPLC>            
                %vars = whos('-file',filename);               
                
                disp(['nacten soubor ' CHilbertMulti.filenameM(filename)]); 
            else
                warning(['soubor CHilbertMulti neexistuje ' CHilbert.filenameH(filename)]);
            end
            obj.hfilename = filename; 
        end 
    end
    methods (Static,Access = public)
        function filenames = ExtractData(PAC,testname,filename,label,overwrite)
            %filenames = ExtractData(PAC,testname,filename,label)
            %podle struct vystupu z funkce CBrainPlot.StructFind, jmena testu, a jmena souboru
            %u kazdeho pacienta vytvori extract pro danou strukturu se jmenemm label
            %vrati seznam filenames, ktery se pak da primo pouzit ve funkcich TestExtract a ImportExtract
            %overwrite urcuje, jestli se maji prepisovat existujici soubory
            if ~exist('overwrite','var'), overwrite = 0; end %defaultne se nemaji prepisovat existujici soubory
            
            pacienti = unique({PAC.pacient}); %trik po delsim usili vygooglovany, jak ziskat ze struct jedno pole
            filenames = cell(numel(pacienti),1);
            
            for p = 1:numel(pacienti)                                
                ipacienti = strcmp({PAC.pacient}, pacienti{p})==1; %indexy ve strukture PAC pro tohoto pacienta
                E = pacient_load(pacienti{p},testname,filename,[],[],[],0);  %pokud spatny testname, zde se vrati chyba
                if ~isempty(E) %soubor muze neexistovat, chci pokracovat dalsim souborem
                    [filename_extract,basefilename_extract] = E.ExtractData([PAC(ipacienti).ch],label,overwrite);
                    filenames{p,1} = filename_extract;
                    disp(['*** OK: ' pacienti{p} ': chns ' num2str([PAC(ipacienti).ch],'%i ') ', ' basefilename_extract]);
                    fprintf('\n'); 
                else
                    disp(['*** nenalezen: ' pacienti{p} ]);
                end                  
            end            
        end
        function filenames = FindExtract(testname,label,filename)
            %filenames = FindExtract(testname,label,filename) 
            %najde existujici extrakty dat podle label            
            if ~exist('filename','var') || isempty(filename) , filename = '*'; end %filename nemusim zadat, pak hledam cokoliv s timto label
            if ~exist('label','var') || isempty(label) , label = '*'; end 
            [pacienti,setup] = pacienti_setup_load( testname );
            filename = strrep(filename,'_CHilb.mat',''); %funkce ExtractData pro zmenu vyzaduje tuhle priponu
            filenames = cell(0,4); %budu vrace vsechny nalezene udaje, nejen filename, kvuli prehledu
            for p = 1:numel(pacienti)
               path = [setup.basedir pacienti(p).folder filesep setup.subfolder filesep];
               files = dir([path filename ' ' label '_Extract.mat']);
               for f = 1:numel(files)
                   files(f).name = [path files(f).name];                   
               end
               files_cell = struct2cell(files)';
               filenames = cat(1,filenames,files_cell(:,1:4)); %prevedu struct na cell array, jen prvni 4 sloupce
            end
            if numel(filenames) > 0
                disp(['OK: nalezeno ' num2str(size(filenames,1)) ' souboru']);
            else
                disp('zadne soubory nenalezeny');
            end
        end        
        function filename2 = filenameM(filename)
            %vraci jmeno souboru s daty tridy CHilbertMulti
           filename=strrep(filename,'_CHilb',''); %odstranim pripony vytvorene pri save
           filename=strrep(filename,'_CiEEG','');
           filename=strrep(filename,'_CHMult','');
           [pathstr,fname,ext] = CiEEGData.matextension(filename);         
           filename2 = fullfile(pathstr,[fname '_CHMult' ext]);
        end
    end
    methods (Static,Access = private)
        function bloky = GetBlocks(epochData)
            %ziska info o blocich 
            %data uz bez treningu
            typ = cell2mat(epochData(:,2)) ;
            b1 = [find(typ(2:end) ~= typ(1:end-1)); size(typ,1)]; %konce bloku
            b0 = [ 1; (b1(1:end-1)+1)]; %zacatky bloku
            bloky = [typ(b0),b1-b0+1,b0,b1];
        end
        function blokyprehazej = ShuffleBlocks(bloky0, bloky1)
            %zjisti indexy k prehazeni bloku, v bloky1=Cil podle bloky0=Zdroj
            blokyprehazej = [bloky0(:,1) zeros(size(bloky0,1),1) bloky0(:,3:4)]; %vzorove poradi bloku , k tomu budu pridavat cisla 
            lastkat = [unique(bloky1(:,1)) zeros(numel(unique(bloky1(:,1))),1)]; %tam si budu ukladat podledni nalezena cisla bloku
            for b = 1:size(blokyprehazej,1)
                
                blok = blokyprehazej(b,1); %cislo podminky z bloku 1 v Zdroj
                ilastkat = find(lastkat(:,1)==blok,1); %index podminky v kastkat
                temp = find(bloky1(lastkat(ilastkat,2)+1:end,1)==blok,1)+lastkat(ilastkat,2); %najde dalsi blok s touto podminkou v Cil
                %find hleda jen v te casti bloky1, takze musim pricist index, odkud jsem hledal
                lastkat(ilastkat,2) = temp;               
                blokyprehazej(b,2)=temp;
            end
            
        end       
        function [str]= basename(filename)
            % vraci filename s koncem path pro identifikaci pacienta
            %[path,basename,ext] = fileparts(filename); %takhle to nechci
            if isempty(filename)
                str = filename;
                return;
            end
            fslash = strfind(filename,'\');
            str = filename(fslash(end-2)+1:end); %dve casti path pred basename
        end
        
    end

    
end

