classdef CiEEGData < matlab.mixin.Copyable 
    
    properties (Access = public)
        d; %double matrix: time x channel, muze byt i time x channel x epoch
        tabs; 
        tabs_orig; %originalni tabs, ktere se zachovaji po epochaci. Downsamplovani se u nich dela; pouzivaji se jen pri hledani EpiEvents
        fs; %vzorkovaci frekvence
        mults; %nepovinne
        header; %nepovinne
        
        samples; %pocet vzorku v zaznamu = rozmer v case
        channels; %pocet kanalu v zaznamu
        epochs;   %pocet epoch
        epochData; %cell array informaci o epochach; epochy v radcich, sloupce: jmeno a cislo kategorie, tab(pondet/odpoved)
        PsyData; %objekt ve formatu CPsyData (PPA, AEDist aj) podle prezentace KISARG
            %pole PsyData.P.data, sloupce, strings, interval, eegfile, pacientid
        epochtime; %delka eventu pre a po event v sekundach , treti cislo je 1, pokud podle responses
        baseline; %delka baseline zacatek a konec v sekundach, vzhledem k podnetu/odpovedi podle epochtime(3)
        CH; %objekt formatu CHHeader s Hammer headerem 
        els; %cisla poslednich kanalu v kazde elektrode
        plotES; % current electrode, second of plot/epoch, range of y values, time range, allels, rangey all els
        plotH;  % handle to plot
        plotRCh = struct; %stavove udaje o grafu PlotResponseCh
        plotEp = struct; %stavove udaje o grafu PlotEpochs
        RjCh; %seznam cisel rejectovanych kanalu
        RjEpoch; %seznam vyrazenych epoch - pouzivam spolecne s RjEpochCh - pro vyrazeni celych epoch
        RjEpochCh; %seznam vyrazenych epoch x kanal - v kazdem kanalu zvlast
        epochTags; %seznam oznacenych epoch
        epochLast; %nejvyssi navstivena epocha
        filename;
        reference; %slovni popis reference original, avg perHeadbox, perElectrode, Bipolar
        yrange = [10 10 50 50]; %minimum y, krok y0, hranice y1, krok y1, viz funkce - a + v hybejPlot
        Wp = {}; %pole signifikanci pro jednotlive kanaly vuci baseline, vysledek  ResponseSearch     
        WpActive=1; %cislo aktivniho pole signifikanci - muze byt spocitano vic statistik
        DE = {}; %trida objektu CEpiEvents - epilepticke eventy ziskane pomoci skriptu spike_detector_hilbert_v16_byISARG
        DatumCas = {}; %ruzne casove udaje, kdy bylo co spocitano. Abych mel historii vypoctu pro zpetnou referenci
        PL = {}; %objekt CPlots
        CS = {}; %objekt CStat
    end
    
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS 
        function obj = CiEEGData(d,tabs,fs,mults,header)
            %konstruktor, parametry d,tabs,fs[,mults,header]
            if (nargin ~= 0 && ~isempty(d))  %konstruktor uplne bez parametru - kvuli CHilbertMulti                
            if ischar(d) && (~exist('fs','var') || isempty(fs)) %pokud je prvni parametr retezec, tak ho beru jako nazev souboru, ktery nactu                
                if ~exist('tabs','var') || isempty(tabs)
                    obj.Load(d);
                else
                    obj.Load(d,tabs); %pokud je druhy parametr, ktery je loadall
                end
            else
                assert(numel(fs)==1,'fs must be a single number');
                assert(size(d,1)== size(tabs,1),'d and tabs have to be the same length'); 
                assert(isfloat(d) || exist('mults','var'),'d neni float, musite zadat parametr mults');
                obj.tabs = tabs;
                obj.tabs_orig = tabs;
                obj.fs = fs;
                if exist('mults','var') && ~isempty(mults)
                    assert(size(mults,1)<= 1 || size(mults,1)==size(d,2),'d and mults have to have same number of channels');
                    assert(~isstruct(mults),'mults cant be structure');
                    obj.mults = mults;
                    obj.d = bsxfun(@times,double(d), mults); %rovnou to roznasobim mults, nechci to resit dodatecne - 24.6.2016
                else
                    obj.mults = ones(1,size(d,2)); %defaultove jednicky pro kazdy kanal
                    obj.d = d;
                end
                [obj.samples,obj.channels, obj.epochs] = obj.DSize();
                if exist('header','var')
                    obj.header = header;                    
                else
                    obj.header = [];
                end
                obj.plotES = [1 1 150 5 0]; %nastavim defaultni hodnoty grafy
                obj.epochLast = 1;
                obj.reference = 'original';
                obj.DatumCas.Created = datestr(now);
                obj.RjEpochCh = false(obj.channels,1); %zatim nejsou zadne epochy              
                disp('vytvoren objekt CiEEGData'); 
            end
            fprintf('epochs: %i, rejected %i (RjEpochCh %i), epochtime: [',obj.epochs,numel(obj.RjEpoch),sum(max(obj.RjEpochCh,[],1)));
            fprintf('%.1f ',obj.epochtime); 
            fprintf('], baseline: [');
            fprintf('%.1f ',obj.baseline);
            fprintf('] \n');
            fprintf('channels: %i, rejected: %i, fs: %i \n',obj.channels, numel(obj.RjCh),obj.fs);                        
            disp(['reference: ' obj.reference]);            
            if ~isempty(obj.DE) 
                    disp(['epievents: ' num2str(size(obj.DE.d,1))]);
            else
                    disp('no epievents');
            end
            if ~isempty(obj.CH)
                if ~isfield(obj.CH.H,'subjName')
                    obj.CH.H.subjName = [obj.CH.H.patientTag ' ' obj.CH.H.patientNick];
                end
                disp(['Hheader ' obj.CH.H.subjName ', triggerch: ' num2str(obj.CH.GetTriggerCh())]);
            else
                disp('no Hheader');
            end 
            if ~isempty(obj.Wp)
                for WpA = 1:numel(obj.Wp)
                    [katstr, opakstr] = obj.KatOpak2Str(WpA);                       
                    disp (['Wilcox stats ' num2str(WpA) ' done, kats: ' katstr ', opakovani: ' opakstr]);
                end
            else
                disp('no Wilcox stats');
            end
            obj.PL = CPlots(); %prazdny objekt na grafy
            obj.CS = CStat; %prazdy objekt na statistiku
            end %(nargin ~= 0) 
        end
        
        function [samples, channels, epochs] = DSize(obj)
            % vraci velikosti pole d - samples, channels, epochs
            samples = size(obj.d,1);
            channels = size(obj.d,2);
            epochs = size(obj.d,3);
        end    
              
        function obj = GetHHeader(obj,H)
            %nacte header z promenne H - 25.5.2016
            if isfield(H,'selCh_H'),  H_channels = size(H.selCh_H,2); else, H_channels = 0; end
            assert(H_channels == size(obj.d,2) || size(H.channels,2)==size(obj.d,2), ...
                 ['nesouhlasi pocet elektrod (data:' num2str(size(obj.d,2)) ',H_channels:' num2str(H_channels) ', header' num2str(size(H.channels,2)) ') - spatny header?']);
            obj.CH = CHHeader(H); %vypocita i selCh_H
            [~, ~, obj.els] = obj.CH.ChannelGroups();  
            assert(max(obj.els)<=size(obj.d,2),['nesouhlasi pocet elektrod (data:' num2str(size(obj.d,2)) ',header:' num2str(max(obj.els)) ') - spatny header?']);
            
            disp(['header nacten: ' obj.CH.PacientTag() ', triggerch: ' num2str(obj.CH.GetTriggerCh())]);
        end
         
        function obj = GetEpiEvents(obj,DE)
           if exist('DE','var') && isstruct(DE)
               obj.DE = CEpiEvents(DE, obj.tabs_orig, obj.fs);
           else
               assert(obj.epochs <= 1, 'nelze pouzit, data jiz jsou epochovana');
               obj.DE = CEpiEvents(obj.d, obj.tabs_orig, obj.fs);    %vytvorim instanci tridy      
           end
           disp(['nacteno ' num2str(size(obj.DE.d,1)) ' epileptickych udalosti']);
        end

        function obj = RejectChannels(obj,RjCh)
            %ulozi cisla vyrazenych kanalu - kvuli pocitani bipolarni reference 
            obj.RjCh = RjCh;
            obj.CH.RejectChannels(RjCh); %ulozim to i do headeru
            disp(['vyrazeno ' num2str(numel(RjCh)) ' kanalu']); 
        end
        
        function obj = RejectEpochs(obj,RjEpoch,RjEpochCh)
            %ulozi cisla vyrazenych epoch - kvuli prevodu mezi touto tridou a CHilbert
            if RjEpoch ~= 0  %muzu takhle vynechat vlozeni vyrazenych epoch              
                obj.RjEpoch = RjEpoch; 
                disp(['globalne vyrazeno ' num2str(numel(RjEpoch)) ' epoch']); 
            end
            if exist('RjEpochCh','var') 
                if ~isempty(RjEpochCh)                 
                    obj.RjEpochCh = RjEpochCh;                                                                         
                    if ~strcmp(obj.reference,'original') && ~isempty(obj.CH.filterMatrix)  %pokud to neni originalni reference                      
                        obj.ChangeReferenceRjEpochCh(obj.CH.filterMatrix); %prepocitam na jinou referenci i RjEpochCh
                    end
                    assert( size(obj.RjEpochCh,1)== size(obj.d,2), ['RjEpochCh ma jiny pocet kanalu (' num2str(size(obj.RjEpochCh,1)) ') nez data (' num2str(size(obj.d,2)) ')']);
                    assert( size(obj.RjEpochCh,2)== size(obj.d,3), ['RjEpochCh ma jiny pocet epoch (' num2str(size(obj.RjEpochCh,2)) ') nez data (' num2str(size(obj.d,3)) ')']);
                    
                    disp(['+ vyrazeno ' num2str(sum(max(RjEpochCh,[],1))) ' epoch s epi udalostmi podle jednotlivych kanalu']);   
                else %takhle muzu vyrazene epochy vymazat
                    obj.RjEpochCh = false(obj.channels,obj.epochs); %zadne vyrazene epochy
                end
            end
            
        end
        
        function [BadChannels,obj] = RjEpochsEpi(obj,NEpi,obrazek)
            %vyradi epochy podle poctu epileptickych udalosti - pokud >= NEpi
            %uz vyrazene epochy rucne nemeni - neoznaci jako spravne
            assert(obj.epochs > 1,'nejsou epochovana data');
            assert(~isempty(obj.DE),'nejsou zadna epi data');
            if ~exist('NEpi','var'), NEpi = []; end
            if ~exist('obrazek','var'), obrazek = 1; end %vykreslim obrazek o poctu vyrazenych epoch v kanalech
            [RjEpoch,RjEpochCh,vyrazeno] =  obj.DE.RejectEpochsEpi(NEpi,obj.CH,obj.epochs,obj.tabs,obj.tabs_orig); %#ok<PROP> 
            
            if isempty(NEpi) %jen pokud jsem RjEpochCh pocital
                obj.RjEpochCh = RjEpochCh; %#ok<PROP>  %prepisu puvodni epochy
                disp(['vyrazeno ' num2str(vyrazeno) ' epoch s epi udalostmi podle jednotlivych kanalu']);
            else
                obj.RjEpoch = unique( [obj.RjEpoch RjEpoch]); %#ok<PROP> %pridam k puvodnim epocham
                disp(['vyrazeno ' num2str(vyrazeno) ' epoch s vice epi udalostmi nez ' num2str(NEpi)]);
            end
            
            if obrazek &&  isempty(NEpi)
                BadChannels = obj.PL.EpochsEpi(obj.RjEpochCh,obj.els,obj.CH); %graf Rejected epochs in individual channels
            end            
            
        end
        
        function [RjEpoch,RjEpochCh]=GetRjEpoch(obj)
            %jen vraci vyrazene epochy pro jejich ulozeni
            RjEpoch = obj.RjEpoch;
            RjEpochCh = obj.RjEpochCh;                      
        end
        function [selCh] = GetSelCh(obj)
            %vraci cisla kanalu vybranych v grafu plotResponseCh, naprikla pro CBrainPLot
            if isprop(obj, 'plotRCh') && isfield(obj.plotRCh,'selCh')
                 selCh = obj.plotRCh.selCh;    %ukladam kvuli selected channels, bez file handelu            
            else
                 selCh = [];
            end
        end
        function obj = SetSelCh(obj,selCh,markno)
            %nastaveni vsechny vybrane kanaly najednou
            if ~exist('markno','var'), markno = 1; end
            if isempty(selCh)
                obj.plotRCh.selCh = zeros(obj.channels,6);                
            elseif selCh(end) >=999 %kdy zadam 999 nebo vyssi cislo, tak se vyberou vsechny kanaly
                obj.plotRCh.selCh = zeros(obj.channels,6);
                obj.plotRCh.selCh(:,1) = ones(obj.channels,1); %prvni mark nastavim vsude 1
            elseif selCh(1) <= -1
                obj.plotRCh.selCh = double(~ obj.plotRCh.selCh); %zeros dela taky double type
            elseif max(selCh) > 1
                obj.plotRch.selCh(selCh,1) = 1; %budu predpokladat ze to jsou cisla kanalu
            elseif size(selCh,2) == 1
                obj.plotRCh.selCh(:,markno) = selCh; %jen jedna znacka najednou
            else
                obj.plotRCh.selCh = selCh; %cele pole najednou
            end
        end
 
        function obj = ExtractEpochs(obj, psy,epochtime,baseline)
            % psy je struct dat z psychopy, 
            % epochtime je array urcujici delku epochy v sec pred a po podnetu/odpovedi: [pred pod podnet=0/odpoved=1]
            % epochuje data v poli d, pridava do objektu: 1. cell array epochData, 2. double(3) epochtime v sekundach, 3. struct psy na tvorbu PsyData
            % upravuje obj.mults, samples channels epochs
            if obj.epochs > 1
                disp('already epoched data');
                return;
            end
            assert(isa(psy,'struct'),'prvni parametry musi by struktura s daty z psychopy');
            if ~exist('baseline','var') || isempty(baseline), baseline = [epochtime(1) 0]; end %defaultni baseline je do 0 sec
            obj.PsyData = CPsyData(psy); %vytvorim objekt CPsyData
            if numel(epochtime)==2, epochtime(3) = 0; end %defaultne epochuji podle podnetu
            obj.epochtime = epochtime; %v sekundach cas pred a po udalosti, prvni cislo je zaporne druhe kladne
            obj.baseline = baseline; %v sekundach
            iepochtime = round(epochtime(1:2).*obj.fs); %v poctu vzorku cas pred a po udalosti
            ibaseline =  round(baseline.*obj.fs); %v poctu vzorku cas pred a po udalosti
            ts_events = obj.PsyData.TimeStimuli(epochtime(3)); %timestampy vsech podnetu/odpovedi
            de = zeros(iepochtime(2)-iepochtime(1), size(obj.d,2), size(ts_events,1)); %nova epochovana data time x channel x epoch  
            tabs = zeros(iepochtime(2)-iepochtime(1),size(ts_events,1)); %#ok<*PROPLC,PROP> %udelam epochovane tabs
            obj.epochData = cell(size(ts_events,1),3); % sloupce kategorie, cislo kategorie, timestamp
            for epoch = 1:size(ts_events,1) %pro vsechny eventy
                izacatek = find(obj.tabs<=ts_events(epoch), 1, 'last' ); %najdu index podnetu/odpovedi podle jeho timestampu
                    %kvuli downsamplovani Hilberta, kdy se mi muze ztratit presny cas zacatku
                    %epochy, beru posledni nizsi tabs nez je cas zacatku epochy
                assert(~isempty(izacatek), ['epocha podle psychopy' num2str(epoch) ' nenalezena v tabs']);
                [Kstring Knum] = obj.PsyData.Category(epoch);    %#ok<*NCOMMA> %jmeno a cislo kategorie
                obj.epochData(epoch,:)= {Kstring Knum obj.tabs(izacatek)}; %zacatek podnetu beru z tabs aby sedel na tabs pri downsamplovani
                for ch = 1:obj.channels %pro vsechny kanaly                    
                    if ibaseline(1)==ibaseline(2)
                        baseline_mean = 0; %baseline nebudu pouzivat, pokud jsem zadal stejny cas jejiho zacatku a konce
                    else
                        baseline_mean = mean(obj.d(izacatek+ibaseline(1) : izacatek+ibaseline(2)-1, ch)); %baseline toho jednoho kanalu, jedne epochy
                    end
                    de(:,ch,epoch) = obj.d( izacatek+iepochtime(1) : izacatek+iepochtime(2)-1,ch) - baseline_mean; 
                    tabs(:,epoch) = obj.tabs(izacatek+iepochtime(1) : izacatek+iepochtime(2)-1); %#ok<PROP>
                end
            end
            obj.d = de; %puvodni neepochovana budou epochovana            
            obj.tabs = tabs; %#ok<PROP>
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
            obj.DatumCas.Epoched = datestr(now);
            obj.RjEpochCh = false(obj.channels,obj.epochs); %zatim zadne vyrazene epochy
            disp(['rozdeleno na ' num2str(obj.epochs) ' epoch']); 
        end
        function obj = ResampleEpochs(obj,newepochtime)
            %resampluje epochy na -0.1 1, pricemz 0-1 je cas odpovedi
            %epochy s delsim casem odpovedi nez je puvodni delka epochy vyradi
            rt = obj.PsyData.ReactionTime(1);
            rtnew = 1; %v kolika sec bude reakcni cas            
            if ~exist('newepochtime','var'), newepochtime = [-0.1 1.1]; end %epocha muze koncit po odpovedi
            newepochlength = newepochtime(2)-newepochtime(1); %delka epochy v sec po resamplovani
            d2 = zeros(round(obj.fs*newepochlength),size(obj.d,2), size(obj.d,3)); %tam budu ukladat eeg data, plati fs=512
            tabs2 = zeros(round(obj.fs*newepochlength),size(obj.d,3)); %nove timestampy
            rjepoch = false(1,obj.epochs); %seznam epoch k vyrazeni, kdyz nesedi cas na konci nebo zacatku
            for ep = 1:obj.epochs
                dd = obj.d(:,:,ep); %jedna epocha, vsechny kanaly a cas
                dd = resample(dd,obj.fs*rtnew,floor(obj.fs*rt(ep))); %resampluju celou epochu, takze rt bude 1s;
                fsnew = size(dd,1)/diff(obj.epochtime(1:2)); %nove spocitana vzorkovaci frekvence dd podle vysledku resample
                    %dve ruzne fs: pro puvodni RT plati fsnew, pro normovane RT=1.0 plati fs=512
                sample0 = round(fsnew*(-obj.epochtime(1))); %cislo vzorku v dd kde je cas 0 (podle nove i stare fs = zarovnani)
                sample01 = round(sample0 - obj.fs*(-newepochtime(1))); %zacatek nove epochy v poli dd - takze normovane fs=512
                if sample01<1 %pokud je zacatek nove epochy newepochtime(1) pred zacatkem puvodni epochy
                    new01 = max(-sample01,1); %to je vlastne -sample1
                    sample01 = 1;
                    rjepoch(ep) = 1; %epochu musim vyradit, data na zacatku jsou jen 0                     
                else
                    new01 = 1; %zacatek nove epochy v poli d2, cas newepochtime(1)
                end                
                sample1 = round(sample0 + obj.fs*newepochtime(2)); %konec nove epochy v poli dd, cas newepochtime(2), normovane fs=512 pro RT=1
                if sample1>size(dd,1) %pokud byl reakcni cas*newepochtime(2) za koncem puvodni epochy                    
                    sample1 = size(dd,1); %zmensim, aby se veslo do epochy newepochtime(1 2)
                    %new1 = new01+size(dd,1)-1;
                    new1 = sample1-sample01+new01; %konec zmensene epochy v poli d2
                    rjepoch(ep) = 1; %epochu musim vyradit, data na konci jsou jen 0                    
                else
                    new1 = size(d2,1); %konec nezmensene epochy se shoduje s velikosti d2
                end
                if sample1-sample01 > new1-new01 %pokud chyba zaokrouhleni
                    sample1 = sample1-( (sample1-sample01)-(new1-new01)); %odectu rozdil delek
                end
                d2(new01:new1,:,ep) = dd(sample01:sample1,:);
                tablimits = [round((-newepochtime(1))*obj.fs) , round((-newepochtime(1))*obj.fs)+size(tabs2,1)-1  ]; %round((newepochlength-newepochtime(1))*obj.fs)                
                    %ten konec mi obcas nevychazel na velikost tabs2, tak to udelam takhle - 26.1.2018
                tabs2(:,ep) = obj.tabs(tablimits(1) : tablimits(2) );                
                %close all;
                %figure,plot(dd);
                %figure,plot(d2(:,:,ep))
            end
            obj.d = d2;
            obj.tabs = tabs2;
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
            obj.epochtime = newepochtime;
            obj.baseline = [newepochtime(1) 0];
            obj.RjEpoch = unique([obj.RjEpoch find(rjepoch)]); %pridam dalsi vyrazene epochy k dosud vyrazenym            
            disp(['resampled ' num2str(obj.epochs) ' epochs to ' num2str(newepochtime) ', rejected new epochs: ' num2str(numel(setdiff(find(rjepoch),obj.RjEpoch)))]);
        end
        function [d,psy_rt,RjEpCh,iEpochy]= CategoryData(obj, katnum,rt,opak,ch)
            %vraci eegdata epoch ve kterych podnet byl kategorie/podminky=katnum + reakcni casy - s uz globalne vyrazenymi epochami
            %Pokud rt>0, vraci epochy serazene podle reakcniho casu 
            %pokud opak>0, vraci jen jedno opakovani obrazku - hlavne kvuli PPA test 
            %vraci i epochy k vyrazeni pro kazdy kanal (uz s globalne vyrazenymi epochami)
            %  vyradit rovnou je nemuzu, protoze pocet epoch v d pro kazdy kanal musi by stejny
            %ch ovlivujen jen RjRpCh, d obsahuje vzdy vsechny kanaly, samply a prislusne epochy
            assert(obj.epochs > 1,'data not yet epoched'); %vyhodi chybu pokud data nejsou epochovana            
            assert(obj.channels == size(obj.RjEpochCh,1),'RjEpochCh: spatny pocet kanalu');
            if exist('opak','var') && ~isempty(opak)
                epochyopak = obj.PsyData.GetOpakovani(); %cislo opakovani pro kazdou epochu
                iOpak = ismember(epochyopak , opak); %epochy jen s timto opakovanim
            else
                iOpak = true(obj.epochs,1);  %vsechny epochy              
            end
            if ~exist('ch','var'), ch = 1:obj.channels; end 
            iEpCh = obj.GetEpochsExclude(ch); %seznam epoch k vyhodnoceni (bez chyb, treningu a rucniho vyrazeni=obj.RjEpoch),channels x epochs ; pro CM data to bude ruzne pro kazdy kanal, jinak stejne pro kazdy kanal
            iEpochy = [ ismember(cell2mat(obj.epochData(:,2)),katnum) , iOpak]; %seznam epoch pro tuto kategorii a toto opakovani - k vyhodnoceni
            d = obj.d(:,:,all(iEpochy,2)); %epochy z teto kategorie a tohoto opakovani = maji ve vsech sloupcich 1            
            RjEpCh = obj.RjEpochCh(ch,all(iEpochy,2)) | ~iEpCh(ch,all(iEpochy,2)); %epochy k vyrazeni u kazdeho kanalu - jen pro epochy teto kategorie katnum - odpovidaji poli d       
            
            if numel(ch)==1 %reakcni cas pocitam, jen kdyz vim pro jaky kanal - 8.6.2018 kvuli CPsyDataMulti
                obj.PsyData.SubjectChange(find(obj.els >= ch,1));
                [~,psy_rt,~,~] = obj.PsyData.GetResponses();
                iEpochyP = iEpochy(1:size(psy_rt,1),:); %psy_rt maji rozmer podle puvodniho poctu epoch pred sloucenim v CM. Uz jsou spravne prehazene. 
                %Kdezto iEpochy obsahuji i vyloucene epochy na konci (pokud u tohohle subjektu nebylo tolik epoch)
                psy_rt = psy_rt(all(iEpochyP,2)); %reakcni casy jen pro vybrane kategorie a opakovani a nevyrazene
            else
                psy_rt = zeros(size(d,3),1); %nulove reakcni casy
            end
            
            if exist('rt','var') && ~isempty(rt) %chci hodnoty serazene podle reakcniho casu               
                [psy_rt, isorted] = sort(psy_rt);
                d = d(:,:,isorted); 
            end   
            iEpochy = all(iEpochy,2) & ~obj.RjEpochCh(ch,:)' & iEpCh(ch,:)'; %jeste pripravim k vystup seznam validnich epoch pro tento kanal - bez vyrazenych 
            % epochy podle podminky &~ epochy s epiaktivitou & epochy bez treningu, chyb a rucniho vyrazeni (RjEpoch)
        end      
        
        function obj = ChangeReference(obj,ref)            
            assert(any(ref=='heb'),'neznama reference, mozne hodnoty: h e b');
            assert(isobject(obj.CH),'Hammer header not loaded');
            
            selCh_H = obj.CH.H.selCh_H; %kopie protoze se mi to zmeni v nasledujicim prikazu         
            obj.CH.ChangeReference(ref); %zmeni referenci u headeru - 18.1.2018            
            % zmena EEG dat v poli d
            if obj.epochs <= 1 %ne epochovana data
                filtData = obj.d(:,selCh_H) * obj.CH.filterMatrix;
                assert(size(filtData,1) == size(obj.d,1),'zmenila se delka zaznamu'); %musi zustat stejna delka zaznamu  
                obj.d=filtData;                
            else %epochovana data
                dd = zeros(obj.samples*obj.epochs,numel(selCh_H));
                for ch = 1:numel(selCh_H) %predelam matici 3D na 2D
                    dd(:,ch) = reshape(obj.d(:,ch,:),obj.samples*obj.epochs,1);
                end                
                filtData = dd(:,selCh_H) * obj.CH.filterMatrix;
                assert(size(filtData,1) == size(dd,1),'zmenila se delka zaznamu'); %musi zustat stejna delka zaznamu  
                obj.d = zeros(obj.samples,size(filtData,2),obj.epochs); %nove pole dat s re-referencovanymi daty
                for ch=1:size(filtData,2) %vratim puvodni 3D tvar matice
                    obj.d(:,ch,:) = reshape(filtData(:,ch),obj.samples,obj.epochs); % !! tohle strasne dlouho trva - ZRYCHLIT
                end
                obj.ChangeReferenceRjEpochCh(obj.CH.filterMatrix); %prepocitam na tuhle  referenci i RjEpochCh
                %pocet elektrod se meni jejen u bipolarni ref, kdyz jsou nektere kanaly na konci vyrazene
            end
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
            if ref=='b' 
                obj.DatumCas.ChangeReference = datestr(now);
                obj.header.RjChOriginal = obj.RjCh; %zazalohuju si puvodni rjch, kvuli referenci
                obj.RjCh = []; %rejectovane kanaly uz byly vyrazeny, ted nejsou zadne  
            end          
           
            obj.filename = []; %nechci si omylem prepsat puvodni data 
            switch ref
                case 'h', obj.reference = 'perHeadbox';
                case 'e', obj.reference = 'perElectrode'; 
                case 'b', obj.reference = 'Bipolar';                    
            end
           
            disp(['reference zmenena: ' obj.reference]); 
        end

        function [iEpCh]=GetEpochsExclude(obj,channels)
            %vraci iEpCh (ch x epoch) = index epoch k vyhodnoceni - bez chyb, treningu a rucniho vyrazeni - pro kazdy kanal zvlast            
            %muzu pouzit parametr channels, pokud chci jen nektere kanaly - ostatni jsou pak prazdne
            if ~exist('channels','var'), channels = 1:obj.channels; end
            iEpCh = zeros(obj.channels,obj.epochs);             
            PsyData = obj.PsyData.copy();  %nechci menit puvodni tridu
            for ch = 1:numel(channels) 
                if isa(obj.PsyData,'CPsyDataMulti') || ch==1 %pokud je to prvni kanal nebo se jedna o data CHilbertMulti s ruznymi subjektu a tedy ruznymi pocty chyb aj                                        
                    PsyData.SubjectChange(find(obj.els >= channels(ch),1)); 
                    chyby = PsyData.GetErrorTrials();                    
                    epochsEx = [chyby , zeros(size(chyby,1),1) ]; %pridam dalsi prazdny sloupec
                    epochsEx(obj.RjEpoch,5)=1; %rucne vyrazene epochy podle EEG  - pro tento kanal 
                    if size(epochsEx,1) < size(iEpCh,2) %pokud ruzny pocet epoch u ruznych subjektu
                        epochsEx = cat(1,epochsEx,ones(size(iEpCh,2) - size(epochsEx,1),5));
                    end
                    iEpCh(channels(ch),:) = all(epochsEx==0,2)'; %index epoch k pouziti                 
                    
                else
                    iEpCh(channels(ch),:) = iEpCh(channels(ch-1),:); %pokud se jedna o CPsyData s jednim subjektem, pro vsechny kanaly to bude stejne
                end
            end
            
        end
            
        function obj = ResponseSearch(obj,timewindow,kats,opakovani,method)
            %projede vsechny kanaly a hleda signif rozdil proti periode pred podnetem
            %timewindow - pokud dve hodnoty - porovnava prumernou hodnotu mezi nimi - sekundy relativne k podnetu/odpovedi
            % -- pokud jedna hodnota, je to sirka klouzaveho okna - maximalni p z teto delky
            %TODO - moznost spojit kategorie 
            assert(obj.epochs > 1,'only for epoched data');                       
            if ~exist('method','var') || isempty(method), method = {'wilcox'}; end  %defaultni metoda statistiky je wilcox test
            if ~iscell(method), method = {method,'chn1'}; end %predelam retezec na cell
            if numel(method) < 2, method{2} = 'chn1'; end %druha polozka bude urcovat, jestli se ma vyhodnocovat vsechny kanaly (chnall), nebo kazdy kanal zvlast (chn1)
            
            iEpCh = obj.GetEpochsExclude(); %ziska seznam Chs x Epochs k vyhodnoceni, neni v tom RjEpochCh
            iEp = true(obj.epochs,1); %musim predat nejaky parametr, ale uz ho ted nepotrebuju, kvuli iEpCh - 8.6.2018
            EEGStat = CEEGStat(obj.d,obj.fs);
            WpA = obj.WpActive; %jen zkratka
            %CELKOVA SIGNIFIKANCE VUCI BASELINE - bez ohledu na kategorie
            baseline = EEGStat.Baseline(obj.epochtime,obj.baseline);
            [Pbaseline,ibaseline,iepochtime,itimewindow] = EEGStat.WilcoxBaseline(obj.epochtime,baseline,timewindow,iEp,obj.RjEpochCh | ~iEpCh);   %puvodni baseline uz v epose nemam        
                %11.12.2017 - pocitam signifikanci hned po konci baseline
                %ibaseline je cast iepochtime pred koncem baseline nebo pred casem 0
            if numel(timewindow) <= 1 %chci maximalni hodnotu p z casoveho okna
                obj.Wp(WpA).D2 = Pbaseline; %pole 2D signifikanci si ulozim kvuli kresleni - cas x channels                
            else
                obj.Wp(WpA).D1 = Pbaseline; %pole 1D signifikanci - jedna hodnota pro kazdy kanal            
            end
            obj.Wp(WpA).Dparams = timewindow; %hodnoty pro zpetnou kontrolu
            obj.Wp(WpA).Dfdr = 1;
            obj.Wp(WpA).DiEp = iEp; %index zpracovanych epoch 
            obj.Wp(WpA).DiEpCh = iEpCh & ~obj.RjEpochCh; %index zpracovanych epochCh, pro ruzne kanaly budou ruzna data je u tridy CHilbertMulti 
            obj.Wp(WpA).epochtime = obj.epochtime;
            obj.Wp(WpA).baseline = obj.baseline; %pro zpetnou kontrolu, zaloha parametru
            
            if exist('opakovani','var') && ~isempty(opakovani) && iscell(opakovani)
                assert(numel(opakovani)<=4,'kategorie opakovani mohou byt maximalne ctyri');
                KATNUM = kats; %cisla kategorii, pro ktere se pocita efekt opakovani
                kats = opakovani;   %POZOR kats se meni na opakovani, abych mohl pouzit kod dole   
                disp('hodnotim opakovani');
            else
                KATNUM = []; %hodnotim rozdily mezi kategorimi, 
            end
            %STATISTIKA PRO JEDNOTLIVE KATEGORIE
            if exist('kats','var') && numel(kats)>1  && numel(timewindow)<= 1 
                %vygeneruju EEG data na statistiku pro kazdou kategorii zvlast
                responsekat = cell(numel(kats),1); %eeg response zvlast pro kazdou kategorii 
                baselinekat = cell(numel(kats),1); %baseline zvlast pro kazdou kategorii 
                rjepchkat = cell(numel(kats),1); %epochyxkanaly k vyrazeni u kazde kategori          
                for k = 1:numel(kats) %pro vsechny kategorie/opakovani
                    if ~isempty(KATNUM) 
                        [katdata,~,RjEpCh] = obj.CategoryData(KATNUM,[],kats{k}); %v kats jsou ted opakovani
                    else
                        [katdata,~,RjEpCh] = obj.CategoryData(cellval(kats,k)); %epochy time*channel*epochs jedne kategorie, uz jsou vyrazeny vyrazene epochy
                    end
                    responsekat{k,1} = katdata( ibaseline(2) - iepochtime(1)+1 :end,:,:); %jen cas po podnetu : cas x channel x epochs; 
                    baselinekat{k,1} = katdata( ibaseline(1) - iepochtime(1)+1 : ibaseline(2) - iepochtime(1),:,:); %jen cas po podnetu : cas x channel x epochs; 
                    rjepchkat{k,1} = RjEpCh;
                end
                %provedu statisticke testy
                if strcmp(method{2},'chnall'), Pbaseline = []; end %pokud chci delat statistiku pres vsechny kanaly, P vuci baseline smazu
                [obj.Wp(WpA).WpKat,obj.Wp(WpA).WpKatBaseline] = EEGStat.WilcoxCat(kats,responsekat,baselinekat,rjepchkat,itimewindow,method{1},Pbaseline);                
                %ulozim parametry
                if ~isempty(KATNUM) %pokud vyhodnocuju opakovani
                    obj.Wp(WpA).kats = KATNUM;    %puvodni kategorie
                    obj.Wp(WpA).opakovani = kats; %v kats jsou ted opakovani
                else
                    obj.Wp(WpA).kats = kats; %ulozim si cisla kategorii kvuli grafu PlotResponseCh
                    obj.Wp(WpA).opakovani = {}; %opakovani nedelam
                end  
            else
                obj.Wp(WpA).kats = kats;
                obj.Wp(WpA).WpKat = cell(0);
                obj.Wp(WpA).opakovani = {};
            end
            obj.DatumCas.ResponseSearch = datestr(now);
        end
        function obj = ResponseSearchMulti(obj,timewindow,stat_kats,opakovani,method)
            %vola ResponseSearch pro kazdy kontrast, nastavi vsechny statistiky
            if ~exist('opakovani','var'), opakovani = []; end
            if ~exist('method','var'), method = []; end
            if iscelldeep(stat_kats) %pokud mam nekolik ruznych statistik na spocitani
                for WpA = 1:numel(stat_kats)
                    obj.SetStatActive(WpA);
                    disp(['pocitam kontrast' cell2str(stat_kats{WpA}) ]);
                    obj.ResponseSearch(timewindow,stat_kats{WpA},opakovani,method);
                end
            else
                obj.E.ResponseSearch(timewindow,stat_kats, opakovani,method);
            end
        end
        function obj = SetStatActive(obj,WpActive)
            WpActive = max(1,min(size(obj.Wp,2)+1,WpActive)); %osetreni na prilis vysoke a nizke cislo            
            if WpActive > numel(obj.Wp)
                disp(['nastavena nova prazdna statistika ' num2str(WpActive)]);
            else
               [katstr, opakstr] = obj.KatOpak2Str(WpActive);                
                disp(['nastavena statistika ' num2str(WpActive) ' s kats: ' katstr ', opakovani: ' opakstr]);            
            end
            obj.WpActive = WpActive;
            
        end
        function katnum = Categories(obj)
            %funkce ktera jen vypise kategorie
            [katnum] = obj.PsyData.Categories(1);
        end
        function Fourier(obj,channels,freq,epochs,method)
           % perform FFT and plot, freq should be max and min to plot
            if ~exist('epochs','var') || isempty(epochs),  epochs = 1:obj.epochs;        end  
            if ~exist('method','var') 
                if numel(epochs) > 1 %'pwelch' nebo 'fft'
                    method = 'fft';
                else
                    method = 'pwelch'; %mensi rozliseni frekvenci, na epochovana data uz nema smysl
                end
            end
            figure('Name','Fourier'); 
            for ch = 1:numel(channels)
                for ep = 1:numel(epochs) %spocitam to pro kazdou epochu zvlast a pak z toho udelam prumer
                    [f,fft_d_abs] = CStat.Fourier(obj.d(:,channels(ch),ep),obj.fs,method);
                    if ep == 1 % u prvni epochy vytvorim pole frekvenci, kdyz uz vim jejich pocet
                      frexs = zeros(numel(epochs),length(f));  % epochy x frekvence - pro fft je pocet freq floor(size(obj.d,1)/2+1);, ale pro pwelch nejak jinak
                    end
                    frexs(ep,:) = f;
                    %frexs(ep,:) = movingmean(frexs(ep,:),1001,2); %nic nepomuze
                end
                frequencies = mean(frexs,1); %prumery pres epochy - kazdy sloupec zvlast                    
                plot(frequencies,fft_d_abs,'.-'); %neplotuju prvni frekvenci, cili DC
    %             loglog(frequencies,fft_d_abs,'.'); %log log plot, kde muze byt i prvni frekvence - viz zhang 2015
    %             semilogy(frequencies,fft_d_abs); % log linearni plot - y = log10 
                hold all;
            end
            xlim(freq);
            %set(gca,'xlim',[0 max(frex)*2])
            title(['Power spectral density by ' method ' in channels ' num2str(channels)]);
            ylabel('PSD [ 10*log10(uV ^{2} / Hz) ]');
            if numel(channels) > 1
                legendCell = cellstr(num2str(channels','Ch=%-d'));
                legend(legendCell);
            end
        end
        function Filter(obj,freq,channels,epoch,vykresli)            
            if ~exist('channels','var') || isempty(channels), channels  = 1:obj.channels;        end 
            if ~exist('epoch','var') || isempty(epoch)     ,  epoch = 1;        end 
            if ~exist('vykresli','var'),  vykresli = 1;        end  %defaultne se dela obrazek
            
            fprintf('channels to filter (z %i):',numel(channels));
            for ch = channels
                fprintf('%i, ',ch);
                dd = squeeze(obj.d(:,ch,epoch));
                dd2 = CStat.FIR(freq,dd,obj.fs); %vyfiltruju data
                obj.d(:,ch,epoch) = dd2;  %ulozim vysledek filtrovani do puvodnich dat              
                if vykresli %vykreslim jenom prvni kanal
                    %charakteristika kanalu pred filtrovanim
                    [frequencies,fft_d_abs] = CStat.Fourier(dd,obj.fs);  
                    figure('Name','Filter effects');   
                    plot(frequencies(2:end),fft_d_abs(2:end),'k.'); %neplotuju prvni frekvenci, cili DC
                    xlim([0 250]);            
                    hold on;                
                    %charakterisika kanalu po filtrovani
                    [frequencies,fft_d_abs] = CStat.Fourier(dd2,obj.fs);   %spocitam znova frekvencni charakteristiku
                    plot(frequencies(2:end),fft_d_abs(2:end),'b.'); %vykreslim modre
                    vykresli = 0;
                    title(['Channel ' num2str(ch)]);
                end
            end
            fprintf('... done\n');
        end
          
        function [obj]= Decimate(obj,podil,rtrim)
            %zmensi data na nizsi vzorkovaci frekvenci podle urceneho podilu, naprilkad 4x 512-128Hz, pokud podil=4
            dd = zeros(ceil(size(obj.d,1)/podil) , size(obj.d,2), size(obj.d,3)); % d uz obsahuje jen svoji prvni pulku a delka je delitelna podilem
            tabs =  zeros(ceil(size(obj.d,1)/podil) , size(obj.d,3));
            fprintf('channels to decimate (z %i):',numel(obj.channels));
            for ch = 1:obj.channels %musim decimovat kazdou elektrodu zvlast
                fprintf('%i, ',ch);
                if obj.epochs == 1
                    dd(:,ch) = decimate(obj.d(:,ch),podil); %na 500 Hz z 8000 Hz
                    if ch==1, tabs = downsample(obj.tabs,podil); end % to delam jen jednou, treba u prvniho kanalu                   
                else
                    for ep = 1:obj.epochs
                        dd(:,ch,ep) = decimate(obj.d(:,ch,ep),podil); %na 500 Hz z 8000 Hz
                        if ch==1, tabs(:,ep) = downsample(obj.tabs(:,ep),podil); end 
                    end
                end                
            end
            obj.tabs_orig = downsample(obj.tabs_orig,podil); % 
            obj.fs = obj.fs/podil;
            obj.d = dd;
            obj.tabs = tabs;
            fprintf('... done\n');
            disp(['decimated to ' num2str(obj.fs) ' Hz']);
            if exist('rtrim','var') && ~isempty(rtrim)
                obj.d = obj.d(1:rtrim,:,:);
                obj.tabs = obj.tabs(1:rtrim,:);
            end
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
        end
        function obj = Resample(obj,fsnew)
            assert(obj.epochs==1, 'muze resamplovat jen neepochovana data');
            obj.d = resample(obj.d,fsnew,obj.fs);
            dv = datevec([obj.tabs(1),obj.tabs(end)]); %rozlozi datetime na komponenty Y M D H M S
            secdif = etime(dv(2,:),dv(1,:)); %pocet vterin mezi zacatkem a koncem tabs
            dates = datenum(dv(1,1),dv(1,2),dv(1,3),dv(1,4),dv(1,5),dv(1,6):(1/fsnew):dv(1,6)+secdif+(1/fsnew)); %pridavam jednu hodnotu na konci, protoze to pocitam z rozdilu
            obj.tabs = dates'; %dates jsou v radku
            obj.tabs_orig = obj.tabs;
            obj.fs = fsnew;
            obj.DatumCas.Resampled = datestr(now);
            disp(['resampled to ' num2str(fsnew) 'Hz']);
        end
        function [katsnames,kombinace,kats] = GetKatsNames(obj)
            %vraci nazvy kategorii ve statistice v aktivnim kontrastu a jejich kombinaci, do intervalyResp aj
           if numel(obj.Wp) >= obj.WpActive
               kats = obj.Wp(obj.WpActive).kats; 
               [katnum, katstr] = obj.PsyData.Categories();
               kombinace = combinator(length(kats),2,'p'); %permutace bez opakovani z poctu kategorii
               kombinace = kombinace(kombinace(:,1)>kombinace(:,2),:); %vyberu jen perumtace, kde prvni cislo je vetsi nez druhe   
               katsnames =  cell(1,numel(kats)+ size(kombinace,1));
               for kat = 1: numel(kats)
                    if iscell(kats(kat)) %mame tu vic kategorii proti vice - na jedne strane kontrastu
                        kknames = cell(1,numel(kats{kat})); %jmena individualnich kategorii na jedne strane kontrastu
                        for kk = 1: numel(kats{kat})
                            kknames{kk}=katstr{kats{kat}(kk)+1}; %katnum jsou od 0, katstr indexovany od 1
                        end
                        katsnames{kat} = strjoin(kknames,'+'); %vice kategorii
                    else
                        katsnames{kat} = katstr{katnum==kats(kat)}; %jde to udelat najednou bez for cyklu?
                    end
               end
               for kat = 1:size(kombinace,1) %cyklus pres vsechny kombinace kategorii
                   katsnames{kat+numel(kats)} = [katsnames{kombinace(kat,1)} 'X' katsnames{kombinace(kat,2)} ];
               end
           else
               disp(['neni vypocitana statistika']);
           end
        end
        function [prumery, MNI,names,intervaly,katsnames,neurologyLabels] = IntervalyResp(obj, intervaly,channels,signum, dofig)
            %vypocita hodnoty v jednotlivych intervalech casu pro jednotlive kategorie i pro celkovy prumer       
            %vykresli graf pro kazdy interval do spolecneho plotu
            %vraci prumery [channels x intervaly x kategorie] a MNI(channels)    
            %signum  = jestli chci jen kat1>kat2 (1), nebo obracene (-1), nebo vsechny (0)
            assert(isfield(obj.Wp(obj.WpActive), 'kats'),'musi byt definovany kategorie podnetu');
            assert(isfield(obj.Wp(obj.WpActive), 'WpKatBaseline'),'musi byt spocitana statistika kategorii');
            if ~exist('intervaly','var') || isempty(intervaly), intervaly = [0.1 obj.epochtime(2)]; end %defaultni epocha je cely interval
            if ~exist('channels','var') || isempty(channels) , channels = 1:obj.channels; end
            if ~exist('signum','var') || isempty(signum) , signum = 0; end %defaultne vraci hodnoty vetsi i mensi v prvni kat
            if ~exist('dofig','var'), dofig = 1; end %defaultne delam obrazek
            [katsnames,kombinace,kats] = obj.GetKatsNames();                                

            %spocitam dynamicky permutace vsech kategorii, pro ktere mam spocitanou statistiku                       
            prumery = zeros(numel(channels),size(intervaly,1),numel(kats)+size(kombinace,1));   % channels x intervaly x kategorie - celkova data a jednotlive kategorie            
           
            if dofig, figure('Name','IntervalyResp'); end
            ploth = zeros(1,max(numel(kats),size(kombinace,1))); %handles na jednotlive ploty, kvuli legende
            for int = 1:size(intervaly,1) 
                legendstr = cell(1,max(numel(kats),size(kombinace,1)));
                if dofig, subplot(min(2,size(intervaly,1)),ceil(size(intervaly,1) /2),int);  end %pro kazdy interval jiny subplot
                %spocitam prumery celkove i za kazdou kategorii v kazdem casovem intervalu
                % dve cisla v kazdem sloupci - od do ve vterinach   
                iintervalyData = min(round((intervaly(int,:)-obj.epochtime(1)).*obj.fs),size(obj.d,1)); % pro katdata kde je na zacatku baseline             
                iintervalyStat = min(round(intervaly(int,:).*obj.fs) , size(obj.Wp(obj.WpActive).WpKat{1,2},1)); % pro statistiku obj.Wp.WpKat, kde na zacatku neni baseline
                iintervalyStat(1) = iintervalyStat(1) + (diff(iintervalyStat)-diff(iintervalyData)); %korekce zaokrouhlovani, posunu zacatek iintervalyStat aby stejne dlouhe jako iintervalyData
                if iintervalyStat(1)<=0, iintervalyStat = iintervalyStat + diff([iintervalyStat(1) 1]); end %aby od 1, kvuli vzorkovani d-size 51 vzorku je 0.7969s
                
                colorskat = {[0 0 0],[0 1 0],[1 0 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1]}; %barvy jako v PlotResponseCh, black, green, red, blue, yellow, aqua, fuchsia
                colorkombinace = {0,1,2,4;0 0 3 5;0 0 0 6};
                iChKats = false(2,numel(channels));  %dva radky pro rozdily vuci baselina a kategorii vuci sobe                                                                          
                
                %nejdriv samotne kategorie
                Pmax = zeros(numel(kats),1); %sbiram maxima kategorii kvuli tomu kde posadit kontrasty mezi kat
                for kat = 1: numel(kats) % cyklus pres kategorie - rozdil vuci baseline
                    [katdata,~,RjEpCh] = obj.CategoryData(cellval(kats,kat)); %time x channels x epochs
                    iCh = min(obj.Wp(obj.WpActive).WpKatBaseline{kat,1}(iintervalyStat(1):iintervalyStat(2),channels),[],1) < 0.05; %kanaly kde je signifikantni rozdil vuci baseline, alespon jednou
                    fiCh = find(iCh); %absolutni cisla kanalu
                    data = zeros(diff(iintervalyData)+1,sum(iCh)); % samples x vybrane kanaly 
                    sub = zeros(1,sum(iCh)); % indexy=cislo samplu maximalnich signif hodnot  
                    Wp = obj.Wp(obj.WpActive).WpKatBaseline{kat,1}(iintervalyStat(1):iintervalyStat(2),iCh); %statistika jen pro vyber kanalu, kde je neco signif
                    for ch = 1:sum(iCh) %musim jet po jednotlivych kanalech kvuli RjEpCh, ch je index v ramci je vybranych kanalu se signif rozdilem, takze fiCh
                        %ted vyberu data jen z nevyrazenych epoch:
                        data(:,ch) = mean(katdata(iintervalyData(1):iintervalyData(2),fiCh(ch),~RjEpCh(fiCh(ch),:)),3); %prumer pres epochy pro jeden kanal, pro nevyrazene epochy pro tento kanal
                        fitime = find(Wp(:,ch)<0.05); %indexy vzorku, kde je signif rozdil
                        %ted vyberu maximalni hodnotu jen ze signifikantnich vzorku - ziskam jeji index v data: 
                        [~,subitime] = max(abs(data(fitime,ch))); % index maximalni absolutni hodnoty se signif rozdilem - jen relativni indexy v ramci fitime
                        sub(ch) = fitime(subitime); %prevedu na absolutni indexy v ramci data(:,ch)
                    end                                      
                    %ted ziskam ty maximalni hodnoty pro vsechny kanaly:
                    ind = sub2ind(size(data),sub,1:size(data,2)); %predelam indexovani na absolutni = ne time x channels, ale 1-n
                    prumery(iCh,int,kat) = data(ind); %max nebo min hodnota z kazdeho kanalu                    
                    P = squeeze(prumery(:,int,kat));  %max/min z kazdeho kanalu                  
                    Pmax(kat) = max(P); %maximum pro kategorii pres vsechny kanaly
                    if dofig
                        ploth(kat) = plot(P','o-','Color',colorskat{kat}); %kreslim tuto kategorii                       
                        hold on;           
                        selCh = find(any(obj.plotRCh.selCh,2)); %indexy jakkoliv oznacenych kanalu
                        if ~isempty(selCh) %pokud existuji nejake vybrane kanaly, vykreslim je plnou barvou
                            plot(selCh,P(selCh)','o','Color',colorskat{kat},'MarkerFaceColor', colorskat{kat});
                        end
                    end
                    iChKats(1,:) = iChKats(1,:) | iCh; %pridam dalsi kanaly, kde je signif odpoved                                        
                    legendstr{kat}=katsnames{kat}; %pridam jmeno kategorie na zacatek [legendstr{k}]
                end                
                
                yKombinace = ceil(max(Pmax)+0.5);
                for kat = 1:size(kombinace,1) %cyklus pres vsechny kombinace kategorii
                    [katdata1, ~, RjEpCh1] = obj.CategoryData(cellval(kats,kombinace(kat,1))); %time x channels x epochs 
                    [katdata2, ~, RjEpCh2] = obj.CategoryData(cellval(kats,kombinace(kat,2))); %druha vyssi kategorie, ktera se bude odecitat od te prvni
                    WpK = obj.Wp(obj.WpActive).WpKat{kombinace(kat,2),kombinace(kat,1)}(iintervalyStat(1):iintervalyStat(2),:); %statistika vsech kanalu
                    WpB = obj.Wp(obj.WpActive).WpKatBaseline{kombinace(kat,1),1}(iintervalyStat(1):iintervalyStat(2),:); %rozdil prvni kategorie vuci baseline
                    dataK = zeros(diff(iintervalyData)+1,obj.channels); %tam budu davat rozdil mezi dvema kategoriemi
                    for ch = 1:obj.channels
                        dataK(:,ch) = mean(katdata1(iintervalyData(1):iintervalyData(2),ch,~RjEpCh1(ch,:)),3) - mean(katdata2(iintervalyData(1):iintervalyData(2),ch,~RjEpCh2(ch,:)),3);
                    end
                    idataK = iff(signum>0, dataK > 0, iff(signum < 0, dataK < 0, true(size(dataK)) ));  % jestli chci vetsi, mensi nebo jakekoliv
                    WpAll = cat(3, WpK<0.05 , WpB<0.05 , idataK); %time x channels x tyhle tri podminky, rozdil vuci baseline, rozdil kategorii a kat1 > kat2
                    iCh = any(all(WpAll,3),1); %vsechny tri podminky, alespon v jednom case                                        
                    fiCh = find(iCh); %absolutni cisla kanalu
                    data = zeros(diff(iintervalyData)+1,sum(iCh)); % samples x vybrane kanaly - tam budu ukladat rozdily mezi kategoriemi
                    sub = zeros(1,sum(iCh)); % indexy=cislo samplu maximalnich signif hodnot                      
                    for ch=1:sum(iCh)
                        %ted vyberu data z nevyrazenych epoch a vypocitam rozdil - time x channels
                        data(:,ch) = dataK(:,fiCh(ch));                        
                        idataK = iff(signum>0, dataK(:,fiCh(ch))>0, iff(signum < 0, dataK(:,fiCh(ch))<0, true(size(dataK,1))));  % jestli chci vetsi, mensi nebo jakekoliv
                        WpAll = [ WpK(:,fiCh(ch))<0.05 , WpB(:,fiCh(ch))<0.05 , idataK];
                        fitime = find(all(WpAll,2)); %indexy vzorku, kde je signif rozdil a prvni kat je vetsi
                        %ted vyberu maximalni hodnotu jen ze signifikantnich vzorku - ziskam jeji index v data: 
                        [~,subitime] = max(abs(data(fitime,ch))); % index maximalni absolutni hodnoty se signif rozdilem - jen relativni indexy v ramci fitime
                        sub(ch) = fitime(subitime); %prevedu na absolutni indexy v ramci data(:,ch)
                    end
                    %ted ziskam ty maximalni hodnoty pro vsechny kanaly:
                    ind = sub2ind(size(data),sub,1:size(data,2)); %predelam indexovani na absolutni = ne time x channels, ale 1-n
                    prumery(iCh,int,kat+numel(kats)) = data(ind); %max nebo min hodnota z kazdeho kanalu
                    P = squeeze(prumery(:,int,kat+numel(kats)));  %max/min z kazdeho kanalu    
                                
                    colorindex = colorkombinace{kombinace(kat,2),kombinace(kat,1)};
                    if dofig %kreslim rozdily mezi odpovedmi pro kategorie                        
                        ph = plot(P'+yKombinace,'o-','Color',colorskat{colorindex}); %kreslim tuto kombinaci kategorii nahoru                        
                        selCh = find(any(obj.plotRCh.selCh,2)); %indexy jakkoliv oznacenych kanalu
                        if ~isempty(selCh) %pokud existuji nejake vybrane kanaly, vykreslim je plnou barvou
                            plot(selCh,P(selCh)'+yKombinace,'o','Color',colorskat{colorindex},'MarkerFaceColor', colorskat{colorindex});
                        end
                        if kat>numel(kats), ploth(kat) = ph; end %pokud je kombinaci vic nez kategorii, ulozim si handle, budu ho potrebovat na legendu
                    end
                    iChKats(2,:) = iChKats(2,:) | iCh ;  %pridam dalsi kanaly, kde je signif odpoved                    
                    legendstr{colorindex}=[legendstr{colorindex} ';' katsnames{kombinace(kat,1)} ' X ' katsnames{kombinace(kat,2)} ];                    
                end  
                
                if dofig                                  
                    title(['interval: ' mat2str(intervaly(int,:))]);
                    xlim([-1 numel(channels)+1]);
                    set(gca, 'XTick', [-1 5:5:numel(channels)+1]); %kompatibilni s 2016a a nizsi %xticks([-1 5:5:numel(channels)+1]);
                    grid on;
                    %vykreslim jmena u signifikatnich kanalu
                    for ch = 1:numel(channels)                        
                        if((iChKats(1,ch) && (ch==1 || ~iChKats(1,ch-1))) || (iChKats(2,ch) && (ch==1 || ~iChKats(2,ch-1))) || ismember(ch-1,obj.els)) %pokud je kanal signif a predchozi neni nebo se jedna o zacatek elektrody
                            th = text(ch,yKombinace*0.6,[num2str(ch) ':' obj.CH.H.channels(ch).name]);
                            if ~verLessThan('matlab','9.0') 
                                th.Rotation = 90;
                            end
                        end
                        if ismember(ch,obj.els)
                            line([ch+.5 ch+.5],[0 yKombinace],'Color',[0 153 255]/255); %kreslim hranice elektrod
                        end
                    end
                    text(0,yKombinace*1.1,'kontrasty mezi kat');
                    text(0,0.1,'kat vuci baseline');
                    legend(ploth,legendstr,'Location','best'); %samo to nejak umisti legendu co nejlepe, temi handely dam legendu jen nekam
                    line([0 numel(channels)],[yKombinace yKombinace],'Color','yellow'); %cara rozsilu kategorii
                    line([0 numel(channels)],[0 0],'Color',[0.8 0.8 .8]); %cara rozsilu kategorii
                end                
               
            end 
            MNI = obj.CH.GetMNI(channels);
            [names,neurologyLabels] = obj.CH.GetChNames(channels);            
            assert(numel(MNI)==size(prumery,1),'MNI a prumery maji jiny pocet kanalu');
        end
        function SelChannelStat(obj,kategorie,marks,add,signum)
            %vybere kanaly podle statistiky, podle vysledku IntervalyResp
            %kategorie jsou cisla odpovidajici katsnames z IntervalyResp - kategorie podnetu a jejich kombinace - maximalne 6,muze byt cell array
            %marks jsou cisla znacek 1-6
            %add=1 znamena, ze se maji pridavat znacky k aktualnim znackam, jinak se puvodni smazou
            %signum se predava do IntervalyResp a znamena znamenko rozdilu - +1,0,-1
            if ~exist('signum','var') || isempty(signum) , signum = 0; end %defaultne vraci hodnoty vetsi i mensi v prvni kat
            if ~exist('kategorie','var') || isempty(kategorie) , kategorie = 1:6; end %ktere kategorie chci oznacit
            if ~exist('marks','var') || isempty(marks) , marks = 1:numel(kategorie); end %kterym maji kategorie odpovidat znackam
            if ~exist('add','var') || isempty(add) , add = 0; end %defaultne prepise stare znaceni 
            assert(numel(kategorie)==numel(marks), 'CiEEGData.SelChannelStat: pocet kategorie a marks musi byt stejny');
            [prumery, ~,~,~,katsnames,~] = obj.IntervalyResp([],[],signum, 0);
            selCh = iff(add,obj.GetSelCh(),zeros(size(prumery,1),6));
            pocty = zeros(1,numel(kategorie));
            katname = cell(1,numel(kategorie));
            
            for kat = 1:numel(kategorie)
                if iscell(kategorie)
                    K = kategorie{kat};    %cisla kategorii
                    KN = cell(1,numel(K)); %KategoryNames - bude tu vic kategorii
                else
                    K = kategorie(kat);   
                    KN = katsnames{kategorie(kat)}; %KategoryNames - jmeno tehle kategorie
                end
                for iK = 1:numel(K) %pro vsechny prvky tehle kategorie - muze jich byt vic pokud kategorie je cellarray
                    %K(iK) je ted cislo kategorie
                    if(K(iK)<=size(prumery,3)) && marks(kat) <= 6
                        iP = prumery(:,1,K(iK))~=0;                    
                        selCh(iP,marks(kat)) = 1;
                        pocty(kat)=pocty(kat) + sum(iP); %kolik vybrano v teto kategorii kanalu
                        if iscell(KN) 
                            KN(iK) = katsnames(K(iK));
                        end
                    end
                end  
                katname(kat) = iff(numel(KN) > 1,{KN},KN);                
            end
            obj.SetSelCh(selCh);
            marks_str = 'fghjkl';    
            [marks,im] = sort(marks); %seradim znadky
            katname = katname(im); %seracim kategorie podle znacek
            %jeste popis marks ulozim, abych mohl pozdeji pouzit do grafu
            if ~add
                obj.plotRCh.selChNames = cell(1,6); %prazdna jmena
                obj.plotRCh.selChNames(marks) = katname;
            else
                for m = 1:numel(marks)  
                    if isempty(obj.plotRCh.selChNames{marks(m)})
                        obj.plotRCh.selChNames{marks(m)} = katname{m};
                    else
                        obj.plotRCh.selChNames{marks(m)} = union(obj.plotRCh.selChNames(marks(m)), katname{m});                    
                    end
                end
            end
            disp([ marks_str(marks) ' = ' cell2str(katname) ', (' num2str(pocty) ')']);            
        end
        %% PLOT FUNCTIONS
        function PlotChannels(obj)  
            %vykresli korelace kazdeho kanalu s kazdym
            if obj.epochs <= 1
                CC = corrcoef(obj.d); %vypocitam a zobrazim korelacni matici kanalu
            else
                dd = zeros(obj.samples*obj.epochs,obj.channels);
                for ch = 1:obj.channels %predelam matici 3D na 2D
                    dd(:,ch) = reshape(obj.d(:,ch,:),obj.samples*obj.epochs,1);
                end                
                CC = corrcoef(dd); 
            end
            figure('Name','Channel Correlations');
            imagesc(CC); 
            
            for j = 1:numel(obj.els)
                line([obj.els(j)+0.5 obj.els(j)+0.5],[1 size(CC,1)],'Color','black');
                line([1 size(CC,1)],[obj.els(j)+0.5 obj.els(j)+0.5],'Color','black');
            end  
            for j = 1:numel(obj.RjCh) %oznacim rejektovane kanaly
                line([obj.RjCh(j) obj.RjCh(j)],[1 size(CC,1)],'Color','white', 'LineWidth',3);
                line([1 size(CC,1)],[obj.RjCh(j) obj.RjCh(j)],'Color','white', 'LineWidth',3);
            end    
            colorbar;
        end
        
        function obj = PlotEpochs(obj,ch,kategories,sortrt)
            %uchovani stavu grafu, abych ho mohl obnovit a ne kreslit novy
            assert(obj.epochs > 1,'only for epoched data');
            if ~exist('ch','var')
                if isfield(obj.plotEp,'ch'), ch =  obj.CH.sortorder(obj.plotEp.ch); %vytahnu cislo kanalu podle ulozeneho indexu
                else,  obj.plotEp.ch = 1; ch =  obj.CH.sortorder(1); end
            else
                obj.plotEp.ch = ch; %tady bude ulozeny index sortorder, parametr ch urcuje index v sortorder
                ch =  obj.CH.sortorder(ch); %promenna ch uz urcuje skutecne cislo kanalu
            end
            if ~exist('kategories','var')
                if isfield(obj.plotEp,'kategories'), kategories = obj.plotEp.kategories;
                else kategories = obj.PsyData.Categories(); obj.plotEp.kategories = kategories; end
            else
                obj.plotEp.kategories = kategories;
            end
            if ~exist('sortrt','var')
                if isfield(obj.plotEp,'sortrt'), sortrt = obj.plotEp.sortrt;
                else, sortrt = 1; obj.plotEp.sortrt = sortrt; end %defaultne radim podle reakcniho casu
            else
                obj.plotEp.sortrt = sortrt;
            end
            if isfield(obj.plotEp,'imgsc') %jestli se ma kreslit obrazek pomoci imagesc, bude se menit mezernikem
                if isempty(obj.plotEp.imgsc) || obj.plotEp.imgsc == 0
                    imgsc = 0; 
                    obj.plotEp.imgsc = 0;
                else
                    imgsc = 1; 
                end
            else
                imgsc = 1; %jeste nic nenastaveno, default je imagesc
                obj.plotEp.imgsc = 1;
            end
            if isfield(obj.plotEp,'fh') && ishandle(obj.plotEp.fh)
                figure(obj.plotEp.fh); %pouziju uz vytvoreny graf
                %clf(obj.plotEp.fh); %graf vycistim
            else
                obj.plotEp.fh = figure('Name','All Epochs','Position', [20, 100, 1200, 300]);
                colormap jet; %aby to bylo jasne u vsech verzi matlabu - i 2016
            end
            clf; 
            T = linspace(obj.epochtime(1),obj.epochtime(2),size(obj.d,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            
            maxy = 0; %budu pocitat za vsechny kategorie
            miny = 0;
            for k=1:numel(kategories)
                katnum = cellval(kategories,k);
                subplot(1,numel(kategories),k);
                [katdata,psy_rt] = obj.CategoryData(katnum,sortrt,[],ch);
                E = 1:size(katdata,3); %cisla epoch - kazdou kategorii muze byt jine                
                D = squeeze(katdata(:,ch,:)); %cas x epochs
                if imgsc
                    imagesc(T,E,D'); %barevny colormap epoch
                else
                    plot(T,D); %normalni plot s epochami pres sebe
                end
                maxy = max([maxy max(max( D ))]);
                miny = min([miny min(min( D ))]);                
                xlabel('Time [s]');                 
                title(obj.PsyData.CategoryName(katnum));
                hold on; 
                if(max(psy_rt)>0) %pokud jsou nejake reakcni casy, u PPA testu nejsou
                    if numel(obj.epochtime)<3 || obj.epochtime(3)==0
                        plot(psy_rt,E,'-k','LineWidth',1); %cara reakcnich casu, nebo podnetu, pokud zarovnano podle reakce      
                    else
                        plot(-psy_rt,E,'-k','LineWidth',1); %cara reakcnich casu, nebo podnetu, pokud zarovnano podle reakce      
                    end
                end                
                if imgsc
                    plot(zeros(size(E,2),1),E,'-k','LineWidth',1); %cara podnetu
                end
            end    
            if isfield(obj.plotEp,'ylim') && numel(obj.plotEp.ylim)>=2 %nactu nebo ulozim hodnoty y
                miny = obj.plotEp.ylim(1); maxy = obj.plotEp.ylim(2);
            else
                obj.plotEp.ylim = [miny maxy];
            end
            for k=1:numel(kategories)
                subplot(1,numel(kategories),k);
                if imgsc
                    caxis([miny,maxy]);
                else
                    ylim([miny,maxy]);
                    line([0 0],[miny maxy ],'Color','black','LineWidth',1);
                end
                if k == 1
                    chstr = iff(isempty(obj.CH.sortedby),num2str(ch), [ num2str(ch) '(' obj.CH.sortedby  num2str(obj.plotRCh.ch) ')' ]);
                    ylabel([ 'Epochs - channel ' chstr]); 
                end %ylabel jen u prniho obrazku
                if k == numel(kategories), colorbar('Position',[0.92 0.1 0.02 0.82]); end
            end
            methodhandle = @obj.hybejPlotEpochs;
            set(obj.plotEp.fh,'KeyPressFcn',methodhandle); 
        end
        
        function PlotCategory(obj,katnum,channel)
            %vykresli vsechny a prumernou odpoved na kategorii podnetu
            %nahrazeno funkcemi PlotResponseCh a PlotEpochs
            d1=obj.CategoryData(katnum,[],[],channel); %epochy jedne kategorie
            d1m = mean(d1,3); %prumerne EEG z jedne kategorie
            T = (0 : 1/obj.fs : (size(obj.d,1)-1)/obj.fs) + obj.epochtime(1); %cas zacatku a konce epochy
            E = 1:size(d1,3); %cisla epoch
            h1 = figure('Name','Mean Epoch'); %#ok<NASGU> %prumerna odpoved na kategorii
            plot(T,d1m(:,channel));
            xlabel('Time [s]'); 
            title([ 'channel ' num2str(channel) ' - ' obj.PsyData.CategoryName(katnum)]);
            h2 = figure('Name','All Epochs');  %#ok<NASGU> % vsechny epochy v barevnem obrazku
            imagesc(T,E,squeeze(d1(:,channel,:))');
            colorbar;
            xlabel('Time [s]');
            ylabel('Epochs');
            title([ 'channel ' num2str(channel) ' - ' obj.PsyData.CategoryName(katnum)]);
            
        end
        
        function obj = PlotElectrode(obj,e,s,range,time)
            %vykresli data (2 sekundy ) z jedne elektrody e od vteriny zaznamu s
            %osa y je v rozmezi [-r +r]
            %zatim jen neepochovana data
            assert(~isempty(obj.els) || ~isempty(obj.CH.els),'je nutne nacist header pomoci GetHHeader');
            if ~exist('e','var') || isempty(e), e= obj.plotES(1); end %cislo elektrody
            if ~exist('s','var') || isempty(s), s= obj.plotES(2); end %cislo epochy nebo vteriny zaznamu
            if ~exist('range','var') || isempty(range)
                range = obj.plotES(3); %150 defaultni rozsah osy y
            end
            if ~exist('time','var') || isempty(time)
                time = obj.plotES(4); %5 sekund defaultni casovy rozsah
            end            
            allels = obj.plotES(5); %jestli se maji zobrazovat vsechny kanaly
            
            if isempty(obj.plotH) || ~ishandle(obj.plotH)
                obj.plotH = figure('Name','Electrode Plot'); %zatim zadny neni, novy obrazek                 
            else
                figure(obj.plotH);  %kreslim do existujiciho plotu
                clf; %smazu graf - jinak mi to ted blbne pri posunu
            end
            
            % -------- nastavim rozsah elektrod k zobrazeni -----------------
            [~,els2plot] = obj.CH.ElsForPlot();
            if  allels==1  %chci zobrazit vsechny elektrody
                elektrodvsade = iff(obj.channels/numel(els2plot) > 6, 5, 8);  %31.8.2016 - chci zobrazovat vzdy pet elektrod, indexy v els jsou tedy 1 6 11
                elsmax = 0; %kolik zobrazim kontaktu - rozliseni osy y v poctu kontaktu
                elsdelsi = [0,els2plot]; %pridam jen nulu na zacatek, kvuli pocitani rozdilu 
                for n = 1 : elektrodvsade : numel(elsdelsi)-elektrodvsade
                    elsmax = max (elsmax , elsdelsi(n+elektrodvsade) - elsdelsi(n)); % pocitam jako maximum z petic elektrod
                end
                pocetsad = ceil(numel(els2plot)/elektrodvsade); %kolik ruznych sad petic elektrod budu zobrazovat, 2 pokud <= 10 els, jinak 3 pokud <=15 els%                 
                emod = mod(e-1,pocetsad);
                if emod==0, elmin=1; else elmin=els2plot(emod*elektrodvsade)+1; end                
                elmaxmax = elmin + elsmax -1 ; % horni cislo el v sade, i kdyz bude pripadne prazdne
                ielmax = find(els2plot <= min(elmaxmax,els2plot(end)) , 1, 'last') ; %horni cislo skutecne elektrody v sade
                elmax = els2plot(ielmax);
                els = els2plot( find(els2plot > elmin, 1,'first' )  : ielmax ); %#ok<PROP> %vyber z els2plot takze horni hranice cisel kontaktu
                els(2,1) = elmin;%#ok<PROP>
                els(2,2:end) = els(1,1:end-1)+1; %#ok<PROP> %doplnim dolni radku - zacatky kazde elektrody
            else
                if e==1, elmin = 1; else elmin = els2plot(e-1)+1; end %index prvni elektrody kterou vykreslit
                elmax = els2plot(e);            % index posledni elektrody kterou vykreslit
                els = [elmax; elmin]; %#ok<PROP>
                elmaxmax = elmax;
            end
            
            % -------- ziskam data k vykresleni do promenne dd -----------------
            time_n = time*obj.fs; %kolik vzorku v case chci zobrazit
            if obj.epochs <= 1 %pokud data jeste nejsou epochovana
                iD = [ (s-1)*obj.fs + 1, min(size(obj.d,1), (s-1)*obj.fs + time_n) ]; %indexy eeg, od kdy do kdy vykreslit
                dd = zeros(elmaxmax-elmin+1,time_n);%data k plotovani - prehodim poradi, prvni jsou kanaly
                dd(1:elmax-elmin+1,1 : iD(2)-iD(1)+1 ) = obj.d(iD(1) : iD(2), elmin:elmax)';
                %dd = obj.d( iD(1) : iD(2),elmin: elmaxmax)' ;   
                t = linspace(iD(1)/obj.fs, iD(2)/obj.fs, iD(2)-iD(1)+1); %casova osa  v sekundach          
            else %pokud data uz jsou epochovana 
                iD = []; %potrebuju to predat jako parametr
                assert(s>=1,'cislo epochy musi byt alespon 1');                
                dd = zeros(elmaxmax-elmin+1,time_n);
                time_nsum = 0; %kolik vzorku v case uz mam v poli dd
                ss = s-1; %cislo kreslene epochy
                while time_nsum < time_n %pridavam hodnoty do pole dd z nekolika epoch, dokud nenaplnim pozadovanou delku
                    ss = ss+1;
                    if time_n - time_nsum >= obj.samples 
                        if ss <= obj.epochs
                            dd(1:elmax-elmin+1,1+time_nsum : time_nsum+obj.samples) = squeeze(obj.d(:, elmin: elmax,ss))';  %data k plotovani - prehodim poradi, prvni jsou kanaly
                        end
                        time_nsum = time_nsum + obj.samples;
                    else %pokud mi nezbyva cela epocha do konce pozadovaneho casoveho rozsahu
                        if ss <= obj.epochs
                            dd(1:elmax-elmin+1,1+time_nsum : time_n ) = squeeze(obj.d(1:time_n - time_nsum, elmin: elmax,ss))';
                        end
                        time_nsum = time_n ;
                    end
                    
                end
                t = linspace(obj.epochtime(1), time+obj.epochtime(1), time_n); %casova osa pres nekolik epoch
            end
            % -------- KRESLIM -----------------
            %kod viz navod zde https://uk.mathworks.com/matlabcentral/newsreader/view_thread/294163            
            mi = repmat(-range,[size(dd,1) 1]); % rozsahu osy y mi:ma - repmat jen zopakuje hodnotu -range podle velikosti dd
            ma = repmat(+range,[size(dd,1) 1]);                        
            shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]); %soucet maxim s nasledujicimi minimy. kumulativne 0-nn
            shift = repmat(shift,1,size(dd,2));
            colors = [ 'b' 'k'];
            c = 0;
            h_els = cell(size(els,2)); %#ok<PROP> %budu si uklada handle plotu, abych je pak dal nahoru
            iel = 1;
            for el = els              %#ok<PROP>
                rozsahel = (el(2):el(1))-els(2,1)+1;  %#ok<PROP>
                rozsahel1 = setdiff(rozsahel, obj.RjCh-els(2,1)+1); %#ok<PROP> %nerejectovane kanaly  - zde se pocitaji od 1 proto odecitam els               
                h_els{iel} = plot(t, bsxfun(@minus,shift(end,:),shift( rozsahel1,:)) + dd( rozsahel1,:) ,colors(c+1) );  
                %lepsi bude je nezobrazovat
                %rozsahel0 = intersect(rozsahel,obj.RjCh); %rejectovane kanaly                
                %if numel(rozsahel0)>0 
                    %hrj = plot(t, bsxfun(@minus,shift(end,:),shift( rozsahel0,:)) + dd( rozsahel0,:) ,'color',[.8 .8 .8] );  
                    %uistack(hrj,'bottom');
                %end
                hold on;
                c = 1-c;
                iel = iel + 1 ;
            end
            %hold off;
            set(gca,'ytick',shift(elmaxmax-elmax+1:elmaxmax-elmin+1,1),'yticklabel',elmax:-1:elmin); %znacky a popisky osy y
            grid on;
            
            ylim([min(min(shift))-range max(max(shift))+range]); %rozsah osy y
            ylabel(['Electrode ' num2str(e) '/' num2str(numel(els2plot)) ]);
            xlabel(['Seconds of ' num2str( round(obj.samples*obj.epochs/obj.fs)) ]);
            if allels==1, ty = -shift(4,1); else ty = -shift(2,1); end %jak muze byt size(shift)=[1,2560] - 119Bucko
            text(t(1),ty,[ 'resolution +/-' num2str(range) 'uV']);         
            xlim([t(1) t(end)]);
            
            % -------- ulozim  handle na obrazek a nastaveni grafu -----------------
            methodhandle = @obj.hybejPlot;
            set(obj.plotH,'KeyPressFcn',methodhandle); 
            
            obj.plotES = [e s range time allels ]; %ulozim hodnoty pro pohyb klavesami
            obj.epochLast = max([s obj.epochLast]); %oznaceni nejvyssi navstivene epochy
            
            for j = 1:elmax-elmin+1
                yshift = shift(end,1)-shift(j,1);
                text(t(end),yshift,[ ' ' obj.CH.H.channels(1,elmin+j-1).neurologyLabel ',' obj.CH.H.channels(1,elmin+j-1).ass_brainAtlas]);
                text(t(1)-size(dd,2)/obj.fs/10,yshift,[ ' ' obj.CH.H.channels(1,elmin+j-1).name]);
                if find(obj.RjCh==elmin-1+j) %oznacim vyrazene kanaly
                    text(t(1),yshift+20,' REJECTED');
                end
            end  
            
            % -------- popisy epoch -----------------
            if obj.epochs > 1
                line([0 0],[shift(1,1)-range  shift(end,1)+range],'Color','m'); %cas 0 - stimulus
                for timex = obj.epochtime(2) : obj.epochtime(2)-obj.epochtime(1)  :time+obj.epochtime(1)
                    line([timex timex],[shift(1,1) shift(end,1)],'Color',[.5 .5 .5],'LineWidth',2); %oddelovac epoch
                    line([timex timex]-obj.epochtime(1),[shift(1,1)-range  shift(end,1)+range],'Color','m'); %cas 0 - stimulus
                end
                titul = ['Epoch ' num2str(s) '/' num2str(obj.epochs)];
                for sj = s:ss %pro vsechny zobrazene epochy
                    if find(obj.RjEpoch==sj) 
                        if sj == s, titul = [titul ' - EXCLUDED'];  end %#ok<AGROW>
                        line([obj.epochtime(1) obj.epochtime(2)]+(sj-s)*(obj.epochtime(2)-obj.epochtime(1)),[shift(1,1) shift(end,1)],'Color','r','LineWidth',2);
                    end
                    if find(obj.epochTags==sj)
                        if sj == s, titul = [titul ' - TAGGED']; end   %#ok<AGROW>
                        line([0 0]+(sj-s)*(obj.epochtime(2)-obj.epochtime(1)),[shift(1,1) shift(end,1)],'Color','g','LineWidth',4);
                    end
                    for el = 1:elmax-elmin+1
                        if ss <= obj.epochs && obj.RjEpochCh(el+elmin-1, sj) %pokud je u tohoto kanalu epocha vyrazena
                            yshift = shift(end,1)-shift(el,1);
                            line([obj.epochtime(1) obj.epochtime(2)]+(sj-s)*(obj.epochtime(2)-obj.epochtime(1)),[yshift yshift], ...
                                'Color',[255 91 71]./255,'LineWidth',1,'LineStyle','-')
                        end
                    end
                end
                title(titul);   
                text(t(end)-((t(end)-t(1))/10),ty,[ 'excluded ' num2str(numel(obj.RjEpoch))]); 
            end
            
            %vykresleni epileptickych eventu
            if ~isempty(obj.DE)
                hold on;
                obj.PL.PlotElectrodeEpiEvents(elmin:elmax,obj.RjCh,obj.DE,obj.tabs,obj.tabs_orig,obj.epochs,obj.samples,obj.epochtime,t,ty,s,time_n,elmaxmax,shift,iD);                
                hold off;
            end
            
%             tohle uplne nejvic zdrzuje z cele funkce            
%             for k= 1 : size(els,2) %#ok<PROP> %
%                     uistack(h_els{k}, 'top'); %dam krivky eeg ulne dopredu
%             end
            
        end
        
        function [responses] = PlotResponses(obj)
            %vykresli uspesnost odpovedi spolu s chybami a vyrazenymi epochami
            obj.PsyData.PlotResponses();
            figure(obj.PsyData.fhR); %kreslim dal do stejneho obrazku
            [resp,rt,kategorie,test] = obj.PsyData.GetResponses();            
            plot(obj.RjEpoch,rt(obj.RjEpoch),'*','MarkerSize',10,'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]); %vykreslim vyrazene epochy
            plot(obj.RjEpoch,kategorie(obj.RjEpoch),'*r','MarkerSize',5); %vykreslim vyrazene epochy
            responses = [resp rt kategorie test];
        end
        
        function obj = PlotResponseCh(obj,ch,kategories,pvalue,opakovani)
            %vykresli odpovedi pro jednotlivy kanal
            %opakovani je cell - maximalne tri hodnoty nebo arrays 
            %kategories 
            assert(obj.epochs > 1,'only for epoched data');
            if ~exist('pvalue','var') || isempty(pvalue) || numel(pvalue)>1 %0 neni isempty
                if isfield(obj.plotRCh,'pvalue'), pvalue = obj.plotRCh.pvalue;
                else pvalue = 0; obj.plotRCh.pvalue = pvalue; end %defaulne se NEzobrazuje krivka p value, ale je mozne ji zobrazit
            else
                obj.plotRCh.pvalue = pvalue;
            end
            if ~exist('ch','var')
                if isfield(obj.plotRCh,'ch')
                    ch = obj.CH.sortorder(obj.plotRCh.ch); %vytahnu cislo kanalu podle ulozeneho indexu
                else
                    ch = obj.CH.sortorder(1); %prvni kanal podle sortorder
                    obj.plotRCh.ch = 1; 
                end
            else
                obj.plotRCh.ch = ch; %tady bude ulozeny index sortorder, parametr ch urcuje index v sortorder
                ch = obj.CH.sortorder(ch); %promenna ch uz urcuje skutecne cislo kanalu
                
            end
            WpA = obj.WpActive; %jen zkratka
            if ~exist('kategories','var') || isempty(kategories) 
                if ~isempty(obj.Wp) && isfield(obj.Wp(WpA), 'kats')
                    kategories = obj.Wp(WpA).kats; %pokud jsou kategorie v parametru, prvni volba je pouzit je ze statistiky
                elseif isfield(obj.plotRCh,'kategories')
                    kategories = obj.plotRCh.kategories; %hodnoty drive pouzite v grafu, ty maji prednost pred statistikou
                elseif ~isempty(obj.Wp) && isfield(obj.Wp(WpA),'kats')
                    kategories = obj.Wp(WpA).kats; %hodnoty pouzite ve statistice
                else
                   if numel(obj.PsyData.Categories())<=4 %uz muzu pouzivat 4 kategorie, kvuli Menrot
                     kategories = obj.PsyData.Categories(); %pokud neni vic nez 3 kategorie, vezmu vsechny
                     obj.plotRCh.kategories = kategories;
                   end %pokud je kategorii vic nez tri, neberu je v uvahu a zobrazim pouze prumer
                end                
            else                
                assert(numel(kategories)<=3,'kategorie mohou byt maximalne tri');
                if ~isempty(obj.Wp) && ~isempty(obj.Wp(WpA).WpKat) && (isempty(obj.Wp(WpA).kats) || ~isequal(obj.Wp(WpA).kats, kategories))
                    disp('Statistika spocitana bez kategorii nebo pro jine kategorie')
                end
                obj.plotRCh.kategories = kategories;    %hodnoty zadane parametrem, ty maji absolutni prednost
            end
            %opakovani obrazku kvuli PPA - 28.9.2016
            if ~exist('opakovani','var') || isempty(opakovani)     
                if ~isempty(obj.Wp) && isfield(obj.Wp(WpA),'opakovani')
                    opakovani = obj.Wp(WpA).opakovani;
                elseif isfield(obj.plotRCh,'opakovani')  %neni zadne drive ulozene
                    opakovani = obj.plotRCh.opakovani; %hodnoty drive pouzite v grafu, ty maji prednost pred statistikou    
                else 
                    opakovani = {};                    
                end
            elseif ~iscell(opakovani) && opakovani == 0 %nulou vyresetuju opakovani, ze se nebude pouzivat
                opakovani = {};
                obj.plotRCh.opakovani = opakovani;  
            else
                assert(numel(opakovani)<=3,'kategorie opakovani mohou byt maximalne tri');
                if ~isempty(obj.Wp) && ~isempty(obj.Wp(WpA).WpKat) && (isempty(obj.Wp(WpA).opakovani) || ~isequal(obj.Wp(WpA).opakovani, opakovani))
                    disp('Statistika spocitana bez opakovani nebo pro jina opakovani')
                end
                obj.plotRCh.opakovani = opakovani;    %hodnoty zadane parametrem, ty maji absolutni prednost
            end
            KATNUM = kategories; % kategorie, ktere chci vykreslovat - vsechny dohromady     
            if ~isempty(opakovani)                
                kategories = opakovani; %POZOR - misto kategorii jsou nyni opakovane - cell array
            end 
            T = linspace(obj.epochtime(1),obj.epochtime(2),size(obj.d,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            if isfield(obj.plotRCh,'fh') && (verLessThan('matlab','9.0') || isvalid(obj.plotRCh.fh)) %isvalid je od verze 2016
                figure(obj.plotRCh.fh); %pouziju uz vytvoreny graf
                clf(obj.plotRCh.fh); %graf vycistim
            else
                if isprop(obj,'label') && ~isempty(obj.label)
                    figurename = ['PlotResponseCh - ' obj.label];
                else
                    figurename = 'PlotResponseCh';
                end
                obj.plotRCh.fh = figure('Name',figurename);
            end
            [ymin ymax] = obj.responseChYLim(iff(~isempty(opakovani),KATNUM,kategories));
            
            %TODO - popisky vic vlevo u zarovnani podle odpovedi
            %TODO vypsat i '( - )' jako neurology label
            %TODO trosku vetsi fonty - i do naseho corelu se bude hodit
            obj.PsyData.SubjectChange(find(obj.els >= ch,1)); %to je tu jen kvuli CHilbertMulti a tedy CPsyDataMulti
            rt = obj.PsyData.ReactionTime(); %reakcni casy podle kategorii, ve sloupcich
            
            %ZACINAM VYKRESLOVAT - NEJDRIV MEAN VSECH KATEGORII
            if ~exist('kategories','var') && ~exist('opakovani','var') %26.5.2017 - jen kdyz neexistuji kategorie
                katdata =  obj.CategoryData(KATNUM,[],[],ch);
                M = mean(katdata(:,ch,:),3);             
                E = std(katdata(:,ch,:),[],3)/sqrt(size(katdata,3)); %std err of mean          
                h_errbar = errorbar(T,M,E,'.','Color',[.6 .6 1]); %nejdriv vykreslim errorbars aby byly vzadu [.8 .8 .8]
                hold on;
                h_mean = plot(T,M,'LineWidth',2,'Color',[0 0 1]);  %prumerna odpoved, ulozim si handle na krivku          
                xlim(obj.epochtime(1:2));
               
                obj.plotRCh.range = [min(M)-max(E) max(M)+max(E)]; %zjistim a ulozim rozsah hodnot pro moznost nastaveni osy y
                if ~isempty(obj.Wp) && isfield(obj.Wp(WpA),'D2') %krivka p hodnot z W testu
                    Tr = linspace(0,obj.epochtime(2),size(obj.Wp(WpA).D2,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
                    if pvalue %pokud chci zobrazovat hodnotu p value jako krivku
                        plot(Tr,obj.Wp(WpA).D2(:,ch),'b:');  %carkovana modra cara oznacuje signifikanci prumeru
                    end
                    y = ymin + (ymax-ymin)*0.2;
                    iWp = obj.Wp(WpA).D2(:,ch) <= 0.05;
                    plot(Tr(iWp),ones(1,sum(iWp))*y,'b.'); %tecky jsou p < 0.05                
                    iWpfirst = find(iWp,1,'first');                 
                    if(numel(iWpfirst)>0) 
                        text(-0.01,y,[ num2str( round(Tr(iWpfirst)*1000)) 'ms']); %cas zacatku signifikance
                        text(-0.18,y,[ 'p=' num2str(CStat.round(min(obj.Wp(WpA).D2(:,ch)),3))]);  %cas zacatku signifikance 
                        line([Tr(iWpfirst) Tr(iWpfirst)],obj.plotRCh.ylim,'Color','blue'); %modra svisla cara u zacatku signifikance
                        text(0.05,y, 'Mean');
                    end
                    iWp = obj.Wp(WpA).D2(:,ch) <= 0.01;
                    plot(Tr(iWp),ones(1,sum(iWp))*y,'b*'); %hvezdicky jsou p < 0.01
                end               
                ylim( [ymin ymax].*1.1);
                katlinewidth = 1;
            else
                obj.plotRCh.range = [0 0]; %zjistim a ulozim rozsah hodnot pro moznost nastaveni osy y
                h_mean = [];
                h_errbar = []; %prazdne handle na obrazky
                katlinewidth = 2;
            end
            
            % POTOM JEDNOTLIVE KATEGORIE
            if exist('kategories','var') || exist('opakovani','var') %kategorie vykresluju jen pokud mam definovane karegorie                   
                hue = 0.8;
                colorskat = {[0 0 0],[0 1 0],[1 0 0],[0 0 1]; [hue hue hue],[hue 1 hue],[1 hue hue],[hue hue 1]}; % prvni radka - prumery, druha radka errorbars = svetlejsi
                h_kat = zeros(numel(kategories),2); 
               
                for k= 1 : numel(kategories) %index 1-3 (nebo 4)
                    if ~isempty(obj.Wp)
                        if iscell(obj.Wp(WpA).kats), kk = obj.Wp(WpA).kats{k}(1);  else,   kk = obj.Wp(WpA).kats(k);    end % aby barvy odpovidaly kategoriim podnetum spis nez kategoriim podle statistiky
                    else                        
                        kk =k-1;
                    end
                    %to se hodi zvlast, kdyz se delaji jen dve kategorie vuci sobe ane vsechny, nebo dve vuci jedne nebo dve vuci dvema
                    colorkatk = [colorskat{1,kk+1} ; colorskat{2,kk+1}]; %dve barvy, na caru a stderr plochu kolem
                    if exist('opakovani','var') && ~isempty(opakovani)
                        opaknum = kategories{k}; %v kategories jsou opakovani k vykresleni, a je to cell array
                        [katdata,~,RjEpCh] = obj.CategoryData(KATNUM,[],opaknum,ch); %eegdata - epochy pro tato opakovani
                    elseif iscell(kategories) %iff tady nefunguje, to by bylo samozrejme lepsi85858
                        katnum = kategories{k}; %cislo kategorie, muze byt cell, pokud vice kategorii proti jedne
                        [katdata,~,RjEpCh] = obj.CategoryData(katnum,[],[],ch); %eegdata - epochy jedne kategorie
                    else
                        katnum = kategories(k); %cislo kategorie, muze byt cell, pokud vice kategorii proti jedne
                        [katdata,~,RjEpCh] = obj.CategoryData(katnum,[],[],ch); %eegdata - epochy jedne kategorie                       
                    end   %7.8.2018 - RjEpCh obsahuje jen aktualni kanal, takze rozmer 1x samples 
                    M = mean(katdata(:,ch,~RjEpCh(1,:)),3);
                    E = std(katdata(:,ch,~RjEpCh(1,:)),[],3)/sqrt(size(katdata,3)); %std err of mean
                    %h_kat(k,2) = errorbar(T,M,E,'.','color',colorskat{2,k}); %nejdriv vykreslim errorbars aby byly vzadu[.8 .8 .8]
                    %h_kat(k,2) = plotband(T, M, E, colorskat{2,k}); %nejlepsi, je pruhledny, ale nejde kopirovat do corelu
                    h_kat(k,2) = ciplot(M+E, M-E, T, colorkatk(2,:)); %funguje dobre pri kopii do corelu, ulozim handle na barevny pas
                    xlim(obj.epochtime(1:2)); 
                    hold on;
                    h_kat(k,1) = plot(T,M,'LineWidth',katlinewidth,'Color',colorkatk(1,:));  %prumerna odpoved,  ulozim si handle na krivku  
                    obj.plotRCh.range = [ min(obj.plotRCh.range(1),min(M)-max(E)) max(obj.plotRCh.range(2),max(M)+max(E))]; %pouziju to pak pri stlaceni / z obrazku                    
                    
                    if ~isempty(obj.Wp) && isfield(obj.Wp(WpA),'WpKat') %signifikance mezi kategoriemi
                        Tr = linspace(obj.Wp(WpA).baseline(2),obj.Wp(WpA).epochtime(2),size(obj.Wp(WpA).D2,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
                        for l = k+1:numel(kategories) %katnum jde od nuly 
                            if iscell(obj.Wp(WpA).kats), colorkatl = obj.Wp(WpA).kats{l}(1)+1; else, colorkatl = obj.Wp(WpA).kats(l)+1; end
                            y = ymin + (ymax-ymin)*(0.3 - (k+l)*0.05)  ; %pozice na ose y
                            if k==1, color=colorskat{1,colorkatl}; else color = colorskat{1,1}; end %green a red jsou proti kategorii 0, cerna je kat 1 vs kat 2
                            if pvalue %pokud chci zobrazovat hodnotu p value jako krivku                                
                                plot(Tr,obj.Wp(WpA).WpKat{k,l}(:,ch), ':','Color',color); %carkovana cara oznacuje signifikanci kategorie vuci jine kategorii
                            end
                            %nejdriv p < 0.05                            
                            iWp = obj.Wp(WpA).WpKat{k,l}(:,ch)  <= 0.05; 
                            plot(Tr(iWp),ones(1,sum(iWp))*y, '*','Color',color); %                        
                            iWpfirst = find(iWp,1,'first');                        
                            if(numel(iWpfirst)>0)                                
                                text(-0.025+Tr(1),y,[ num2str(round(Tr(iWpfirst)*1000)) 'ms']);  %cas zacatku signifikance 
                                text(-0.06+Tr(1),y,[ 'p=' num2str(CStat.round(min(obj.Wp(WpA).WpKat{k,l}(:,ch)),3))]);  %cas zacatku signifikance 
                                line([Tr(iWpfirst) Tr(iWpfirst)],obj.plotRCh.ylim,'Color',color); %modra svisla cara u zacatku signifikance                                
                            end                            
                            %potom jeste p < 0.01
                            iWp = obj.Wp(WpA).WpKat{k,l}(:,ch)  <= 0.01;                          
                            plot(Tr(iWp),ones(1,sum(iWp))*y,  '*','Color',color); %
                            % jmena kategorii vypisuju vzdy
                            if exist('opakovani','var') && ~isempty(opakovani)   %pokud vyhodnocuju opakovani
                                kat1name =  obj.PsyData.OpakovaniName(kategories{l});
                                kat2name =  obj.PsyData.OpakovaniName(kategories{k});
                                kat3name =  [ ' (' obj.PsyData.CategoryName(obj.Wp(WpA).kats) ')' ]; %jmeno kategorie obrazku, ze ktere se opakovani pocitalo
                            elseif iscell(kategories)
                                kat1name =  obj.PsyData.CategoryName(kategories{l});
                                kat2name =  obj.PsyData.CategoryName(kategories{k});
                                kat3name = '';
                            else
                                kat1name =  obj.PsyData.CategoryName(kategories(l));
                                kat2name =  obj.PsyData.CategoryName(kategories(k));
                                kat3name = '';
                            end
                            text(0.04+obj.Wp(WpA).epochtime(1),y, ['\color[rgb]{' num2str(colorskat{1,colorkatl}) '}' kat1name ...
                                    '\color[rgb]{' num2str(color) '} *X* '  ...
                                    '\color[rgb]{' num2str(colorkatk(1,:)) '}' kat2name kat3name]);                            
                        end                                              
                    end
                   
                    
                    if ~isempty(obj.Wp) && isfield(obj.Wp(WpA),'WpKatBaseline') %signifikance vuci baseline
                            iWpB = obj.Wp(WpA).WpKatBaseline{k,1}(:,ch)  <= 0.05; %nizsi signifikance
                            y = ymin + (ymax-ymin)*(0.28 - (k+2)*0.05)  ;
                            plot(Tr(iWpB),ones(1,sum(iWpB))*y, '.','Color',colorkatk(1,:),'MarkerSize',5); % 
                            iWpB = obj.Wp(WpA).WpKatBaseline{k,1}(:,ch)  <= 0.01; % vyssi signifikance
                            %y = ymin + (ymax-ymin)*(0.28 - (k+2)*0.05)  ;
                            plot(Tr(iWpB),ones(1,sum(iWpB))*y, 'p','Color',colorkatk(1,:),'MarkerSize',5); % 
                            if exist('opakovani','var') && ~isempty(opakovani)
                                kat2name =  obj.PsyData.OpakovaniName(kategories{k}); %pokud vyhodnocuju opakovani
                            elseif iscell(kategories)
                                kat2name =  obj.PsyData.CategoryName(kategories{k});
                            else
                                kat2name =  obj.PsyData.CategoryName(kategories(k));
                            end
                            text(0.04+obj.Wp(WpA).epochtime(1), y, ['\color[rgb]{' num2str(colorkatk(1,:)) '}' kat2name ' vs.baseline'] );
                            line([Tr(1) Tr(end)],[y y]+(ymax-ymin)*0.03 ,'Color',[0.5 0.5 0.5]);
                                %kazde jmeno kategorie jinou barvou
                            if pvalue %pokud chci zobrazovat hodnotu p value jako krivku
                               plot(Tr,obj.Wp(WpA).WpKatBaseline{k,1}(:,ch), '-.','Color',colorskat{1,k}); %teckovana cara oznacuje signifikanci kategorie vuci baseline
                            end
                    end
                    %cara reakcnich casu pro tuhle kategorii
                    y=ymax-(ymax-ymin)*0.07*k;
                    line([quantile(rt(:,k),0.25) quantile(rt(:,k),0.75)],[y y],'Color',colorkatk(1,:)); %cara kvantilu 
                    plot(median(rt(:,k)),y,'o','Color',colorkatk(1,:)); %median
                end
                y = (ymax-ymin)*0.1  ; %pozice na ose y
                text(0.04+obj.Wp(WpA).epochtime(1),y,['stat ' num2str(obj.WpActive) '/' num2str(numel(obj.Wp)) '-'  cell2str(obj.PsyData.CategoryName(kategories,[])) ]); %vypisu cislo aktivni statistiky a jmena kategorii
                for k= 1 : numel(kategories) %index 1-3
                    uistack(h_kat(k,1), 'top'); %dam krivky prumeru kategorii uplne dopredu
                    uistack(h_kat(k,2), 'bottom'); %dam krivky errorbars uplne dozadu
                end
                ylim( [ymin ymax].*1.1);
            end
            if ~isempty(h_mean)
                uistack(h_errbar, 'top');
                uistack(h_mean, 'top'); %uplne nahoru dam prumer vsech kategorii
            end 
            
            chstr = iff(isempty(obj.CH.sortedby),num2str(ch), [ num2str(ch) '(' obj.CH.sortedby  num2str(obj.plotRCh.ch) ')' ]);
            title(['channel ' chstr '/' num2str(obj.channels) ' - ' obj.PacientID()], 'Interpreter', 'none'); % v titulu obrazku bude i pacientID napriklad p132-VT18
            text(-0.1,ymax*.95,[ obj.CH.H.channels(1,ch).name ' : ' obj.CH.H.channels(1,ch).neurologyLabel ',' obj.CH.H.channels(1,ch).ass_brainAtlas]);
            if  isfield(obj.CH.H.channels,'MNI_x') %vypisu MNI souradnice
                text(-0.1,ymax*.90,[ 'MNI: ' num2str(obj.CH.H.channels(1,ch).MNI_x) ', ' num2str(obj.CH.H.channels(1,ch).MNI_y ) ', ' num2str(obj.CH.H.channels(1,ch).MNI_z)]);
            else
                text(-0.1,ymax*.90,'no MNI');
            end
            if isfield(obj.CH.H.channels,'seizureOnset') %vypisu epilepticke info
                seizureOnset    = iff(obj.CH.H.channels(1,ch).seizureOnset==1,'seizureOnset','-');
                interictalOften = iff(obj.CH.H.channels(1,ch).interictalOften==1,'interictalOften','-');
                if isfield(obj.CH.H.channels,'rejected')
                    rejected = iff( ~isempty(obj.CH.H.channels(1,ch).rejected==1),'rejected','-');
                else
                    rejected = '';
                end
                text(-0.1,ymax*.85,['epiinfo: ' seizureOnset ',' interictalOften ',' rejected ]);
            else
                text(-0.1,ymax*.85,['no epiinfo']);
            end
            if isprop(obj,'plotRCh') && isfield(obj.plotRCh,'selCh') && any(obj.plotRCh.selCh(ch,:),2)==1                
                klavesy = 'fghjkl'; %abych mohl vypsat primo nazvy klaves vedle hvezdicky podle selCh
                text(-0.18,ymax*.95,['*' klavesy(logical(obj.plotRCh.selCh(ch,:)))], 'FontSize', 12,'Color','red');
            end
            if isprop(obj,'label') && ~isempty(obj.label)
                text(-0.1,ymax*.78,strrep(obj.label,'_','\_'), 'FontSize', 10,'Color','blue'); 
            end            
            methodhandle = @obj.hybejPlotCh;
            set(obj.plotRCh.fh,'KeyPressFcn',methodhandle);          
        end        
            
        function obj = PlotResponseP(obj)
            %vykresli signifikanci odpovedi u vsech kanalu EEG vypocitanou pomoci ResponseSearch                                
            assert(obj.epochs > 1,'only for epoched data');
            T = 0:0.1:obj.epochtime(2); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            WpA = obj.WpActive; %jen zkratka
            if isfield(obj.Wp(WpA),'D1') %první 2D plot
                figure('Name','W plot 1D');
                isignif = obj.Wp(WpA).D1<0.05;
                plot(find(~isignif),obj.Wp(WpA).D1(~isignif),'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','b');
                hold on;
                plot(find(isignif),obj.Wp(WpA).D1(isignif),'o','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r');
                ylim([0 0.1]);
                view(-90, 90); %# Swap the axes
                set(gca, 'ydir', 'reverse'); %# Reverse the y-axis 
                set(gca, 'xdir', 'reverse'); %# Reverse the x-axis 
                for e = 1:numel(obj.els) %hranice elektrod a jmeno posledniho kontaktu
                    line([obj.els(e)+0.5 obj.els(e)+0.5],[0 0.1],'color',[.5 0.5 0.5]);
                    text(obj.els(e)-1,-0.01,obj.CH.H.channels(1,obj.els(e)).name);
                end
                for ch=1:obj.channels
                    if obj.Wp(WpA).D1(ch)<0.1 %anatomicka jmena u signif kontaktu
                        text(ch,0.102,obj.CH.H.channels(1,ch).neurologyLabel);
                    end
                end
            end
            if isfield(obj.Wp(WpA),'D2') %isprop(obj,'Wp') && isfield(obj.Wp,'D2')
                figure('Name','W map 2D');
                imagesc(T,1:obj.channels,1 - obj.Wp(WpA).D2', [0.95 1]); %mapa, od p>0.05 bude modra barva 
                axis ij;
                ylabel('channels');
                xlabel('time [s]');
                colorbar;
                for e = 1:numel(obj.els) %hranice elektrod a jmeno posledniho kontaktu
                    line([T(1) T(end)],[obj.els(e)+0.5 obj.els(e)+0.5],'color','w');
                    text(-T(end)/10,obj.els(e)-1,obj.CH.H.channels(1,obj.els(e)).name);
                end
                for ch=1:obj.channels %anatomicka jmena u signif kontaktu
                    if any(obj.Wp(WpA).D2(:,ch)<0.05)
                        text(T(end)*1.07,ch,[ ' ' obj.CH.H.channels(1,ch).neurologyLabel ',' obj.CH.H.channels(1,ch).ass_brainAtlas],'FontSize',8);
                        text(T(end)/10,ch,num2str(ch),'color','w');
                    end
                end
                text(0,-3,['rejected ' num2str(numel(obj.RjEpoch)) ' epochs']);
            end
        end 
        function PlotEpiEvents(obj)
            %vykresli pocty epileptickych events u jednotlivych kanalu. U epochovanych dat i pocty epoch s epi udalostmi
            %since 27.4.2017
            %TODO husteji popisovat elektrody - u tech co nad 10% vypisovat cislo a jmeno u vrcholu
            assert(isobject(obj.CH),'Hammer header not loaded');
            assert(isobject(obj.DE),'Epievents not loaded');
            [evts,names,epochs,evts_nonseeg] = obj.DE.CountEpiEvents(obj.CH,obj.epochs,obj.tabs,obj.tabs_orig); %#ok<PROP>            
            [seizureOnset,interIctal]=obj.CH.GetSeizures(); %indexy interictalOften a seizureOnset kanalu 
            figure('Name','Epievents in individual channels');
            if obj.epochs > 1
                subplot(2,1,1);                
            end
            
            plot(evts,'.-');
            hold on;
            plot(seizureOnset,repmat(40,1,numel(seizureOnset)),'o','Color','red','MarkerFaceColor', 'red'); %cervene o jsou seizure onset
            plot(interIctal,repmat(40,1,numel(interIctal)),'o','Color','magenta','MarkerSize', 10);
            set(gca,'xtick',obj.CH.H.selCh_H ,'xticklabel',names); %znacky a popisky osy y
            for el = 1:numel(obj.els)-1
                line([obj.els(el) obj.els(el)]+1,[0 max(evts)],'Color',[0.5 0.5 0.5]);
            end            
            xlabel('channels');
            title('pocet epi eventu celkove');
            disp(['celkem vykresleno epieventu: ' num2str(sum(evts)) ', + nevykresleno ' num2str(evts_nonseeg) ' v non eeg kanalech']);
            
            if obj.epochs > 1 %epochovana data - druhy graf
                subplot(2,1,2);
               
                plot(epochs./obj.epochs,'.-'); %#ok<PROP>
                set(gca,'xtick',obj.CH.H.selCh_H ,'xticklabel',names); %znacky a popisky osy y
                for el = 1:numel(obj.els)-1
                    line([obj.els(el) obj.els(el)]+1,[0 1],'Color',[0.5 0.5 0.5]);
                end
                ylim([0 1]);
                line([1 obj.CH.H.selCh_H(end)],[0.30 0.30],'Color','red');
                xlabel('channels');
                spatne = find(epochs./obj.epochs >= 0.30); %#ok<PROP>
                disp('kanaly s pocet epi epoch >= 0.3');
                disp(spatne');
                title('podil epoch s epi eventy');
            end
            
        end            
                
        %% SAVE AND LOAD FILE
        function obj = Save(obj,filename)   
            %ulozi veskere promenne tridy do souboru
            if ~exist('filename','var')
                filename = obj.filename;
                assert( ~isempty(filename), 'no filename given or saved before');
            else
                obj.filename = filename;
            end
            d = obj.d;                      %#ok<PROP,NASGU>            
            tabs = obj.tabs;                %#ok<PROP,NASGU>
            tabs_orig = obj.tabs_orig;      %#ok<PROP,NASGU>
            fs = obj.fs;                    %#ok<PROP,NASGU>            
            header = obj.header;            %#ok<PROP,NASGU>
            sce = [obj.samples obj.channels obj.epochs]; %#ok<NASGU>
            if isobject(obj.PsyData)
                if isa(obj.PsyData,'CPsyDataMulti')
                    PsyData = obj.PsyData; %#ok<NASGU> %v tomhle pripade budu ukladat cely objekt
                    PsyDataP = []; testname = ''; %#ok<NASGU>
                else
                    PsyDataP = obj.PsyData.P;       %#ok<NASGU>         %ulozim pouze strukturu P
                    testname = obj.PsyData.testname; %#ok<NASGU>
                    PsyData = [];%#ok<NASGU>
                end
            else
                PsyDataP = []; %#ok<NASGU>
                testname = ''; %#ok<NASGU>
            end
            epochtime = obj.epochtime;      %#ok<PROP,NASGU>
            baseline = obj.baseline;        %#ok<PROP,NASGU>
            CH_H=obj.CH.H;                  %#ok<NASGU>            
            CH_filterMatrix = obj.CH.filterMatrix; %#ok<NASGU>  
            els = obj.els;                  %#ok<PROP,NASGU>
            plotES = obj.plotES;            %#ok<PROP,NASGU>
            selCh = obj.GetSelCh();      %#ok<NASGU>
            %plotH = obj.plotH;             %#ok<PROP,NASGU> %plotH je blbost ukladat, vytvori se novy, jen to brani vice grafum - 14.6.2016
            RjCh = obj.RjCh;                %#ok<PROP,NASGU>
            RjEpoch = obj.RjEpoch;          %#ok<PROP,NASGU>
            RjEpochCh = obj.RjEpochCh;      %#ok<PROP,NASGU>
            epochTags = obj.epochTags;      %#ok<PROP,NASGU>
            epochLast = obj.epochLast;      %#ok<PROP,NASGU>
            reference = obj.reference;      %#ok<PROP,NASGU>
            epochData = obj.epochData;      %#ok<PROP,NASGU>
            Wp = obj.Wp;                    %#ok<PROP,NASGU>
            DE = obj.DE;                    %#ok<PROP,NASGU>
            DatumCas = obj.DatumCas;        %#ok<PROP,NASGU>
            if isa(obj,'CHilbertMulti'), label = obj.label; else label = []; end %#ok<NASGU>
            [pathstr,fname,ext] = CiEEGData.matextension(filename);        
            filename2 = fullfile(pathstr,[fname ext]);
            save(filename2,'d','tabs','tabs_orig','fs','header','sce','PsyDataP','PsyData','testname','epochtime','baseline','CH_H','els',...
                    'plotES','selCh','RjCh','RjEpoch','RjEpochCh','epochTags','epochLast','reference','epochData','Wp','DE','DatumCas', 'label', ...
                    'CH_filterMatrix','-v7.3');  
            disp(['ulozeno do ' filename2]); 
        end
        function obj = Load(obj,filename,~,~)
            % nacte veskere promenne tridy ze souboru
            assert(exist(filename,'file')==2, 'soubor s daty neexistuje, nejde o data tridy CHilbert?');
            vars = whos('-file',filename) ;
            assert(ismember('d', {vars.name}), 'soubor neobsahuje promennou d, nejde o data tridy CHilbert?'); 
            load(filename,'d','tabs','tabs_orig','fs','header','sce','epochtime','els','plotES','RjCh','RjEpoch','epochTags','epochLast','reference');            
            obj.d = d;                      %#ok<CPROPLC,CPROP,PROP> 
            obj.tabs = tabs;                %#ok<CPROPLC,CPROP,PROP> 
            obj.tabs_orig = tabs_orig;      %#ok<CPROPLC,CPROP,PROP> 
            obj.fs = fs;                    %#ok<CPROPLC,CPROP,PROP>          
            obj.mults = ones(1,size(d,2));  %#ok<CPROPLC,CPROP,PROP> 
            obj.header = header;            %#ok<CPROPLC,CPROP,PROP> 
            obj.samples = sce(1); obj.channels=sce(2); obj.epochs = sce(3); %sumarni promenna sce
            vars = whos('-file',filename);
            if ismember('PsyDataP', {vars.name}) %ulozena pouze struktura P z PsyData
                load(filename,'PsyDataP'); 
                if ~isempty(PsyDataP)
                    obj.PsyData = CPsyData(PsyDataP); %vytvorim objekt psydata ze struktury                    
                end
            end
            if (~isprop(obj,'PsyData') || isempty(obj.PsyData)) && ismember('PsyData', {vars.name}) %pokud jsem nevytvoril objekt v predchozim if
                load(filename,'PsyData');                 
                obj.PsyData = PsyData ; %#ok<CPROPLC>  %  %drive ulozeny objekt, nez jsem zavedl ukladani struct nebo CPsyDataMulti                
            end
            if ismember('testname', {vars.name})
                load(filename,'testname');
                obj.PsyData.GetTestName(testname); %#ok<CPROPLC> %  %zjisti jmeno testu
            else
                obj.PsyData.GetTestName(''); %#ok<CPROPLC> %  %zjisti jmeno testu
            end
            if obj.epochs > 1
                if ismember('epochData', {vars.name}), load(filename,'epochData');  obj.epochData = epochData;   end  %#ok<CPROPLC,CPROP,PROP> 
                if ismember('baseline',  {vars.name}), load(filename,'baseline');   obj.baseline = baseline;   end    %#ok<CPROPLC,CPROP,PROP>      
                load(filename,'epochtime');                
                obj.epochtime = epochtime;      %#ok<CPROPLC,CPROP,PROP>               
            else
                obj.epochtime = [];
                obj.baseline = [];
                obj.epochData = [];
            end
            if ismember('CH_H', {vars.name})
                load(filename,'CH_H');      obj.CH = CHHeader(CH_H);
                [~, ~, obj.els] = obj.CH.ChannelGroups();  
            else
                load(filename,'CH');
                obj.CH = CH; %#ok<CPROPLC,CPROP,PROP> %  %drive ulozeny objekt, nez jsem zavedl ukladani struct
            end 
            if ismember('CH_filterMatrix', {vars.name})
                load(filename,'CH_filterMatrix');      obj.CH.filterMatrix = CH_filterMatrix;                              
            end 
            
            if ismember('Wp', {vars.name})
                load(filename,'Wp');      obj.Wp = Wp; %#ok<CPROPLC,CPROP,PROP>
            else
                obj.Wp = struct;
            end
            if ismember('DE', {vars.name}) %1.9.2016
                load(filename,'DE');      obj.DE = DE; %#ok<CPROPLC,CPROP,PROP>
            else
                obj.Wp = struct;
            end
            if ismember('DatumCas', {vars.name}) %7.4.2017
                load(filename,'DatumCas');      obj.DatumCas = DatumCas; %#ok<CPROPLC,CPROP,PROP>
            else
                obj.DatumCas = {};
            end
            if ismember('RjEpochCh', {vars.name}) %17.7.2017
                load(filename,'RjEpochCh');      obj.RjEpochCh = RjEpochCh; %#ok<CPROPLC,CPROP,PROP> 
            else
                obj.RjEpochCh = false(obj.channels,obj.epochs); %zatim zadne vyrazene epochy
            end
            obj.els = els;                  %#ok<CPROPLC,CPROP,PROP> 
            obj.plotES = plotES;            %#ok<CPROPLC,CPROP,PROP> 
            if ismember('selCh', {vars.name}) %nastaveni grafu PlotResponseCh
                load(filename,'selCh'); obj.plotRCh.selCh = selCh;          %#ok<CPROPLC,CPROP,PROP> 
            end
            if isempty(obj.plotRCh.selCh) %kdyz to je prazdne, tak to pak zlobi, musi byt zeros
                obj.SetSelCh([]); %nastavim prazdne - zadne vybrane kanaly
            end
            %obj.plotH = plotH;             %#ok<CPROPLC,CPROP,PROP> 
            obj.RjCh = RjCh;                %#ok<CPROPLC,CPROP,PROP>     
            obj.RjEpoch = RjEpoch;          %#ok<CPROPLC,CPROP,PROP> 
            if exist('epochTags','var'),  obj.epochTags = epochTags;   else obj.epochTags = []; end         %#ok<CPROPLC,CPROP,PROP>     
            if exist('epochLast','var'),  obj.epochLast = epochLast;   else obj.epochLast = []; end         %#ok<CPROPLC,CPROP,PROP> 
            if exist('reference','var'),  obj.reference = reference;   else obj.reference = 'original'; end  %#ok<CPROPLC,CPROP,PROP>  %14.6.2016            
            obj.filename = filename;
            if isa(obj,'CHilbertMulti') && ismember('label', {vars.name}), load(filename,'label'); obj.label = label; end %#ok<NASGU>
            disp(['nacten soubor ' filename]); 
        end
    end
    %% staticke metody
    methods (Static,Access = public)
        function [pathstr,fname,ext] = matextension(filename)
            [pathstr,fname,ext] = fileparts(filename);
            if strcmp(ext,'.mat')==false || numel(ext)<1
               fname = [fname ext]; %pokud pripona neni mat, pridam ji na konec jmena a vytvorim priponu mat
               ext = '.mat';
            end 
        end
    end
    %% privatni metody
    methods  (Access = private)
        function obj = hybejPlot(obj,~,eventDat)    %pohybuje grafem PlotElectrode       
           switch eventDat.Key
               case 'rightarrow' 
                   if obj.epochs == 1
                       rightval = obj.plotES(2)+obj.plotES(4);
                       maxval = round(obj.samples/obj.fs);                       
                   else
                       rightval = obj.plotES(2); 
                       maxval = obj.epochs; %maximalni cislo epochy k vykresleni
                   end
                   if( rightval < maxval)   %pokud je cislo vteriny vpravo mensi nez celkova delka                        
                        obj.PlotElectrode(obj.plotES(1),obj.plotES(2)+1,obj.plotES(3),obj.plotES(4));
                   end
               case 'leftarrow'
                   if(obj.plotES(2))>1 %pokud je cislo vteriny vetsi nez 1
                        obj.PlotElectrode(obj.plotES(1),obj.plotES(2)-1,obj.plotES(3),obj.plotES(4));
                   end
               case 'pagedown' %posunuti o vetsi usek doprava
                   if obj.epochs == 1
                       step = 5; %5 sekund, o kolik se chci posunout doprava
                       rightval = obj.plotES(2)+obj.plotES(4);
                       maxval = round(obj.samples/obj.fs)-step;
                   else
                       step = obj.PlottedEpochs(); %pocet zobrazenych celych epoch 
                       rightval = obj.plotES(2)+step-1; %cislo zobrazene epochy vpravo
                       maxval = obj.epochs-step; %pocet epoch
                        
                   end
                   if( rightval < maxval)   %pokud je cislo vteriny/epochy vpravo mensi nez celkova delka
                        obj.PlotElectrode(obj.plotES(1),obj.plotES(2)+step,obj.plotES(3),obj.plotES(4));
                   end
               case 'pageup' %posunuti o vetsi usek doleva
                   if obj.epochs == 1
                       step = 5; %5 sekund, o kolik se chci posunout doleva
                   else
                       step = obj.PlottedEpochs();
                   end
                   if(obj.plotES(2))>step %pokud je cislo vteriny vetsi nez 1
                        obj.PlotElectrode(obj.plotES(1),obj.plotES(2)-step,obj.plotES(3),obj.plotES(4));
                   end
               case 'home'     % na zacatek zaznamu              
                        obj.PlotElectrode(obj.plotES(1),1,obj.plotES(3),obj.plotES(4)); 
               case 'end'     % na konec zaznamu        
                   if obj.epochs == 1
                        obj.PlotElectrode(obj.plotES(1),obj.samples/obj.fs - obj.plotES(4),obj.plotES(3),obj.plotES(4));  
                   else
                        obj.PlotElectrode(obj.plotES(1),obj.epochs - obj.PlottedEpochs()+1,obj.plotES(3),obj.plotES(4));     
                   end                       
               case 'downarrow'
                   [~,els2plot] = obj.CH.ElsForPlot();
                   if(obj.plotES(1))<numel(els2plot) %pokud je cislo elektrody ne maximalni
                        obj.PlotElectrode(obj.plotES(1)+1,obj.plotES(2),obj.plotES(3),obj.plotES(4));
                   end                   
               case 'uparrow'
                   if(obj.plotES(1))>1 %pokud je cislo elektrody vetsi nez 1
                        obj.PlotElectrode(obj.plotES(1)-1,obj.plotES(2),obj.plotES(3),obj.plotES(4));
                   end               
               case {'add' ,  'equal'}     %signal mensi - vetsi rozliseni %u terezy na notebooku 
                   if obj.plotES(3)>=obj.yrange(3), pricist = obj.yrange(4);
                   else pricist = obj.yrange(2);                   
                   end
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3)+pricist,obj.plotES(4));                
               case {'subtract' , 'hyphen'} %signal vetsi - mensi rozliseni   %u terezy na notebooku  
                   if obj.plotES(3)>obj.yrange(3), odecist = obj.yrange(4);
                   elseif obj.plotES(3)>obj.yrange(1), odecist = obj.yrange(2);
                   else odecist = 0;
                   end
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3)-odecist,obj.plotES(4));
               
               case 'delete' %epoch exclusion
                   s = obj.plotES(2);
                   if find(obj.RjEpoch== s)                      
                        obj.RjEpoch = obj.RjEpoch(obj.RjEpoch~=s); %vymazu hodnoty s                        
                   else
                        obj.RjEpoch = [obj.RjEpoch  obj.plotES(2)]; %pridam hodnotu s
                   end   
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3),obj.plotES(4));
               case 'space' %epoch tag  - oznaceni jednolivych epoch 
                   obj.AddTag();
                   obj.PlotElectrode();
               case {'numpad4','a'}  %predchozi oznacena epocha
                   s = obj.plotES(2);
                   prevTag = obj.epochTags(obj.epochTags < s);
                   prevDel = obj.RjEpoch(obj.RjEpoch < s);
                   if numel(prevTag) > 0 || numel(prevDel)>0
                     obj.PlotElectrode(obj.plotES(1),max([prevTag prevDel]),obj.plotES(3),obj.plotES(4));
                   end                   
               case {'numpad6','d'} %dalsi oznacena epocha
                   s = obj.plotES(2);
                   nextTag = obj.epochTags(obj.epochTags > s);
                   nextDel = obj.RjEpoch(obj.RjEpoch > s);
                   nextLast = obj.epochLast(obj.epochLast > s);
                   if numel(nextTag) > 0 || numel(nextDel)>0 || numel(nextLast)>0
                    obj.PlotElectrode(obj.plotES(1),min([nextTag nextDel nextLast]),obj.plotES(3),obj.plotES(4));                    
                   end
                   
               case 'return'  %prehazuje mezi zobrazeni jednotlivych elektrod a cele poloviny elektrod                   
                   obj.plotES(5) = 1-obj.plotES(5);
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3),obj.plotES(4));
               otherwise
                   disp(['You just pressed: ' eventDat.Key]);                      
           end
        end
        function obj = hybejPlotCh(obj,~,eventDat)  
           %reaguje na udalosti v grafu PlotResponseCh
           switch eventDat.Key
               case {'rightarrow','c'} %dalsi kanal
                   obj.PlotResponseCh( min( [obj.plotRCh.ch + 1 , obj.channels]));
               case 'pagedown' %skok o 10 kanalu dopred
                   obj.PlotResponseCh( min( [obj.plotRCh.ch + 10 , obj.channels]));
               case {'leftarrow','z'} %predchozi kanal
                   obj.PlotResponseCh( max( [obj.plotRCh.ch - 1 , 1]));
               case 'pageup' %skok 10 kanalu dozadu
                   obj.PlotResponseCh( max( [obj.plotRCh.ch - 10 , 1]));
               case 'home' %skok na prvni kanal
                   obj.PlotResponseCh( 1);
               case 'end' %skok na posledni kanal
                   obj.PlotResponseCh( obj.channels);
               case {'multiply','8'} %hvezdicka na numericke klavesnici, nebo hvezdicka nad osmickou
                   %dialog na vlozeni minima a maxima osy y
                   answ = inputdlg('Enter ymax and min:','Yaxis limits', [1 50],{num2str(obj.plotRCh.ylim)});
                   if numel(answ)>0  %odpoved je vzdy cell 1x1 - pri cancel je to cell 0x0
                       if isempty(answ{1}) || any(answ{1}=='*') %pokud vlozim hvezdicku nebo nic, chci znovy spocitat max a min
                           obj.plotRCh.ylim = [];
                       else %jinak predpokladam dve hodnoty
                           data = str2num(answ{:});  %#ok<ST2NM>
                           if numel(data)>= 2 && data(1)< data(2) %pokud nejsou dve hodnoty, nedelam nic
                             obj.plotRCh.ylim = [data(1) data(2)];
                           end
                       end
                   end
                   obj.PlotResponseCh( obj.plotRCh.ch); %prekreslim grafy
               case {'divide','slash'} %lomeno na numericke klavesnici - automaticke meritko na ose y
                   obj.plotRCh.ylim = obj.plotRCh.range; %spocitalo se pri volani PlotResponseCh
                   obj.PlotResponseCh( obj.plotRCh.ch); %prekreslim grafy
               case 'space' %zobrazi i prumerne krivky - vsechny epochy a vsechny frekvence
                   if isa(obj,'CHilbert'), obj.PlotResponseFreq(obj.plotRCh.ch,obj.Wp(obj.WpActive).kats); end %vykreslim vsechna frekvencni pasma
                   obj.PlotEpochs(obj.plotRCh.ch,obj.Wp(obj.WpActive).kats); %vykreslim prumery freq u vsech epoch
                   figure(obj.plotRCh.fh); %dam puvodni obrazek dopredu
               case 'return' %zobrazi obrazek mozku s vybranych kanalem                      
                   if isprop(obj,'label') && ~isempty(obj.label), label = obj.label; else, label = ''; end
                   obj.CH.ChannelPlot2D(obj.plotRCh.ch,obj.plotRCh,@obj.PlotResponseCh,label);  %vykreslim obrazek mozku s vybranym kanalem
                   figure(obj.plotRCh.fh); %dam puvodni obrazek dopredu
               case {'add' ,  'equal','f'}     % + oznaceni kanalu
                   obj.SelChannel(obj.CH.sortorder(obj.plotRCh.ch));
                   obj.PlotResponseCh( obj.plotRCh.ch); %prekreslim grafy
               case {'g','h'}     % + oznaceni kanalu Mark 2-6
                   obj.SelChannel(obj.CH.sortorder(obj.plotRCh.ch),eventDat.Key - 'f' +1 ); %g je 2, f je 1
                   obj.PlotResponseCh( obj.plotRCh.ch); %prekreslim grafy
               case {'j','k','l'}     % + oznaceni kanalu Mark 2-6
                   obj.SelChannel(obj.CH.sortorder(obj.plotRCh.ch),eventDat.Key - 'f' ); %g je 2, f je 1
                   obj.PlotResponseCh( obj.plotRCh.ch); %prekreslim grafy
               case {'numpad6','d'}     % skok na dalsi oznaceny kanal   
                   if isfield(obj.plotRCh,'selCh')
                       selCh = find(any(obj.plotRCh.selCh,2)); %seznam cisel vybranych kanalu
                       iselCh = find(ismember(obj.CH.sortorder,selCh)); %indexy vybranych kanalu v sortorder
                       chn2 = iselCh(find(iselCh>obj.plotRCh.ch,1)); %dalsi vyznaceny kanal
                       obj.PlotResponseCh( iff(isempty(chn2),obj.plotRCh.ch,chn2) ); %prekreslim grafy                        
                   end                   
               case {'numpad4','a'}     % skok na predchozi oznaceny kanal
                   if isfield(obj.plotRCh,'selCh')
                       selCh = find(any(obj.plotRCh.selCh,2)); %seznam cisel vybranych kanalu
                       iselCh = find(ismember(obj.CH.sortorder,selCh)); %indexy vybranych kanalu v sortorder
                       chn2 =  iselCh(find(iselCh < obj.plotRCh.ch,1,'last')) ;
                       obj.PlotResponseCh( iff(isempty(chn2),obj.plotRCh.ch,chn2) ); %prekreslim grafy
                   end
               case {'numpad9','e'}     % skok na dalsi kanal s nejakou signifikanci
                   chsignif = obj.ChannelsSignif();
                   chn2 = chsignif(find(chsignif>obj.plotRCh.ch,1)); %nasledujici signif kanaly
                   obj.PlotResponseCh( iff(isempty(chn2),obj.plotRCh.ch,chn2) ); %prekreslim grafy  
               case {'numpad7','q'}     % skok na predchozi kanal s nejakou signifikanci
                   chsignif = obj.ChannelsSignif(); %seznam  kanalu s nejakou signifikanci
                   chn2 = chsignif(find(chsignif<obj.plotRCh.ch,1,'last')); %nasledujici signif kanaly
                   obj.PlotResponseCh( iff(isempty(chn2),obj.plotRCh.ch,chn2) ); %prekreslim grafy 
               case {'numpad8','w'}     % zvyseni cisla aktivni statistiky  
                   if numel(obj.Wp)> obj.WpActive
                       obj.WpActive = obj.WpActive + 1;
                       obj.PlotResponseCh();                       
                   end
               case {'numpad5','s'}     % snizeni cisla aktivni statistiky 
                   if obj.WpActive > 1
                       obj.WpActive = obj.WpActive - 1;
                       obj.PlotResponseCh();
                   end
               case 'period'     % prepinani razeni kanalu
                   sortorder0 = obj.CH.sortorder; %musi si ulozit stare razeni, abych potom nasel ten spravny kanal
                   obj.CH.NextSortChOrder();                   
                   obj.PlotResponseCh(find(obj.CH.sortorder==sortorder0(obj.plotRCh.ch))); %#ok<FNDSB> %takhle zustanu na tom stejnem kanale 
               case 'r' %roc krivka
                   obj.CS.AUCPlot(obj.plotRCh.ch,obj);
                   figure(obj.plotRCh.fh); %dam puvodni obrazek dopredu
               otherwise
                   disp(['You just pressed: ' eventDat.Key]);
           end
        end
        
        function obj = hybejPlotEpochs(obj,~,eventDat)
            %reaguje na klavesy v PlotEpochs
            switch eventDat.Key 
                case 'multiply' %hvezdicka na numericke klavesnici
                   %dialog na vlozeni minima a maxima osy y
                   answ = inputdlg('Enter ymax and min:','Yaxis limits', [1 50],{num2str(obj.plotEp.ylim)});
                   if numel(answ)>0  %odpoved je vzdy cell 1x1 - pri cancel je to cell 0x0
                       if isempty(answ{1}) || any(answ{1}=='*') %pokud vlozim hvezdicku nebo nic, chci znovy spocitat max a min
                           obj.plotEp.ylim = [];
                       else %jinak predpokladam dve hodnoty
                           data = str2num(answ{:});  %#ok<ST2NM>
                           if numel(data)>= 2 %pokud nejsou dve hodnoty, nedelam nic
                             obj.plotEp.ylim = [data(1) data(2)];
                           end
                       end
                   end
                   obj.PlotEpochs( obj.plotEp.ch); %prekreslim grafy
               case {'divide','slash'} %lomeno na numericke klavesnici - automaticke meritko na ose y 
                   obj.plotEp.ylim = [];
                   obj.PlotEpochs( obj.plotEp.ch); %prekreslim grafy
               case {'subtract' , 'hyphen'} %signal minus - razeni epoch podle rt on/off   %u terezy na notebooku  
                   obj.plotEp.sortrt = 1 - obj.plotEp.sortrt;  %zmenim sortrt
                   obj.PlotEpochs( obj.plotEp.ch); %prekreslim grafy  
               case 'space' %prepinani mezi color plot a plot pres sebe                   
                   if isfield(obj.plotEp,'imgsc')
                       obj.plotEp.imgsc = 1 - obj.plotEp.imgsc; %prepinam druh grafu
                   else
                       obj.plotEp.imgsc = 0;
                   end
                   obj.PlotEpochs( obj.plotEp.ch); %prekreslim grafy                     
            end
                   
        end
        function obj = AddTag(obj,s)
           if ~exist('s','var')
                s = obj.plotES(2);
           end
           if find(obj.epochTags== s)                      
                obj.epochTags = obj.epochTags(obj.epochTags~=s); %vymazu hodnoty s
           else
                obj.epochTags = [obj.epochTags  obj.plotES(2)]; %pridam hodnotu s                
           end  
        end
        
        function epochs = PlottedEpochs(obj)
            %vraci pocet zobrazenych celych epoch 
            epochs = floor(obj.plotES(4) / (obj.epochtime(2) - obj.epochtime(1))); 
        end
        function [ymin, ymax, obj] =  responseChYLim(obj,kategories)
            %nastavi max a min grafu PlotResponseCh();
            if isfield(obj.plotRCh,'ylim') && numel(obj.plotRCh.ylim)>=2 && obj.plotRCh.ylim(1)~=obj.plotRCh.ylim(2) %pokud mam drive ulozene ylim
                    ylim( obj.plotRCh.ylim .*1.1); %udelam rozsah y o 10% vetsi
                    ymax = obj.plotRCh.ylim(2);
                    ymin = obj.plotRCh.ylim(1);
                else
                    ymax = 0; ymin = 0;
                    for k=1:numel(kategories)
                        if iscell(kategories) %tady iff nefunguje, vraci mi to chybu
                            katdata = obj.CategoryData(kategories{k}); %pokud je kategorii spolecne vic
                        else
                            katdata = obj.CategoryData(kategories(k)); %epochy jedne kategorie
                        end
                        channels = 1:obj.channels; %#ok<PROP>
                        channels(ismember(channels, [obj.RjCh obj.CH.GetTriggerCh()]))=[]; %#ok<PROP> %vymazu rejectovana a triggerovane channels 
                        ymax = max([ ymax max(mean(katdata(:,channels,:),3))]); %#ok<PROP>
                        ymin = min([ ymin min(mean(katdata(:,channels,:),3))]); %#ok<PROP>
                    end  
                    ymin = ymin - 0.15*(ymax-ymin); %pridam patnact procent na napisy dole
                    %ylim( [ymin ymax].*1.1); %udelam rozsah y o 10% vetsi
                    assert(ymin~=ymax,'nemuzu urcit rozsah osy y - prazdna data?');
                    obj.plotRCh.ylim = [ymin ymax];
            end
        end
        function id = PacientID(obj)
            %vraci oznaceni pacienta, bud z CPsyData nebo z CHHeader
            id= obj.PsyData.PacientID();
            if isempty(id) || numel(id)<=1
                id = obj.CH.PacientTag();
            end
        end
        function [obj] = ChangeReferenceRjEpochCh(obj,filterMatrix)
            %kod Nada 2017-12-07 - prepocitani RjEpochCh na bipolarni referenci            
            RjEpochCh = obj.RjEpochCh(1:size(filterMatrix,1),:)';  %u zadneho z pacientu jsem nenasel trigger channel uprostred kanalu, vzdy je na konci. To by jinak byl problem            
            filterMatrix(filterMatrix<0) = 0; %oprava pro bipolarni referenci - chci mit v kazdem slouci je jednu 1
            filterMatrix(filterMatrix>0) = 1; %pridano kvuli jine = ele a head referenci
            RjEpochCh = RjEpochCh * filterMatrix; 
            RjEpochCh(RjEpochCh >= 2) = 1;
            obj.RjEpochCh = RjEpochCh'; %vyrazeni kazdeho kanalu puvodni reference znamena vyrazeni dvou kanalu bipolarni reference 
        end
        function [katstr, opakstr] = KatOpak2Str(obj,WpA)
            if ~exist('WpA','var'), WpA = 1; end
            if isfield(obj.Wp(WpA),'opakovani')
                opakstr = cell2str(obj.Wp(WpA).opakovani);
%                 if iscell(obj.Wp(WpA).opakovani)
%                     strjoin(cellfun(@num2str,obj.Wp(WpA).opakovani,'un',0));  %tohle jsem nejak vygooglil, ale nevypisuje to slozene zavorky
%                 else
%                     opakstr = num2str(obj.Wp(WpA).opakovani); 
%                 end
            else 
                opakstr = 'no'; 
            end
            if isfield(obj.Wp(WpA),'kats')
               katstr = cell2str(obj.Wp(WpA).kats); %moje nova funkce 6.3.2018
%                if iscell(obj.Wp(WpA).kats)
%                    katstr = strjoin(cellfun(@num2str,obj.Wp(WpA).kats,'un',0)); %chci nejak vypsat kategorie v cell array
%                else
%                    katstr = num2str(obj.Wp(WpA).kats);  
%                end
            else
                katstr = 'no';
            end
        end        
        function chsignif = ChannelsSignif(obj)
           %vrati seznam kanalu s nejakou signifikanci z WpKatBaseline nebo WpKat
           chsignif = []; %seznam kanalu se signifikanci           
           for kat = 1:numel(obj.Wp(obj.WpActive).kats)
               chsignif = sort(unique([chsignif find(min(obj.Wp(obj.WpActive).WpKatBaseline{kat,1},[],1)<0.05) ]));                       
           end
           for kat1 = 1:numel(obj.Wp(obj.WpActive).kats)
               for kat2 = kat:numel(obj.Wp(obj.WpActive).kats)
                   chsignif = sort(unique([chsignif find(min(obj.Wp(obj.WpActive).WpKat{kat1,kat2},[],1)<0.05) ]));                                 
               end               
           end
        end
    end
    methods  (Access = protected)
        function obj = SelChannel(obj,ch,markno)
            %vybere nebo odebere jeden kanal
            if ~exist('markno','var'), markno = 1; end
            assert(markno <= 6,'moc vysoke cislo znacky');
            if ~isfield(obj.plotRCh,'selCh')
                obj.plotRCh.selCh = zeros(obj.channels,6); %6 ruznych znacek moznych
                obj.plotRCh.selCh(ch,markno) = 1; %prvni vybrany kanal
            else
                obj.plotRCh.selCh(ch,markno) = 1 - obj.plotRCh.selCh(ch,markno); %pridam kanal k vyberu , nebo odeberu             
            end
        end
        
        function cpObj = copyElement(obj)
            % Override copyElement method: to copy also property objects
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Make a deep copy of the copyable classes
            if isobject(obj.PsyData), cpObj.PsyData = copy(obj.PsyData); end
            if isobject(obj.CH), cpObj.CH = copy(obj.CH); end
            if isobject(obj.DE), cpObj.DE = copy(obj.DE); end
            if isobject(obj.PL), cpObj.PL = copy(obj.PL); end
      end
    end
    
end

