classdef CiEEGData < handle
    %CEEGDATA Trida na praci s datama ve formatu ISARG od Petra Jezdika
    %   Kamil Vlcek, FGU AVCR, since 2016 04
    
    properties (Access = public)
        d; %double matrix: time x channel, muze byt i time x channel x epoch
        tabs; 
        tabs_orig; %originalni tabs, ktere se zachovaji po epochaci. Downsamplovani se u nich dela
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
        RjEpoch; %seznam vyrazenych epoch
        epochTags; %seznam oznacenych epoch
        epochLast; %nejvyssi navstivena epocha
        filename;
        reference; %slovni popis reference original, avg perHeadbox, perElectrode, Bipolar
        yrange = [10 10 50 50]; %minimum y, krok y0, hranice y1, krok y1, viz funkce - a + v hybejPlot
        Wp = {}; %pole signifikanci pro jednotlive kanaly vuci baseline, vysledek  ResponseSearch     
        DE = {}; %trida objektu CEpiEvents - epilepticke eventy ziskane pomoci skriptu spike_detector_hilbert_v16_byISARG
        DatumCas = {}; %ruzne casove udaje, kdy bylo co spocitano. Abych mel historii vypoctu pro zpetnou referenci
        PL = {}; %objekt CPlots
    end
    
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS 
        function obj = CiEEGData(d,tabs,fs,mults,header)
            %konstruktor, parametry d,tabs,fs[,mults,header]
            if ischar(d) && (~exist('tabs','var') || isempty(tabs)) %pokud je prvni parametr retezec, tak ho beru jako nazev souboru, ktery nactu                
                obj.Load(d);
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
                disp('vytvoren objekt CiEEGData'); 
            end
            disp(['epochs: ' num2str(obj.epochs) ', rejected: ' num2str(numel(obj.RjEpoch)), '; channels: ' num2str(obj.channels) ', rejected: ' num2str(numel(obj.RjCh)) ...
                ', fs: ' num2str(obj.fs)]); 
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
                if isfield(obj.Wp,'opakovani'), opakstat = num2str(obj.Wp.opakovani); else opakstat = 'no'; end
                if isfield(obj.Wp,'kats'), kats = num2str(obj.Wp.kats); else kats = 'no'; end
                disp (['Wilcox stats done, kats: ' kats ', opakovani: ' opakstat]);
            else
                disp('no Wilcox stats');
            end
            obj.PL = CPlots();
        end
        
        function [samples, channels, epochs] = DSize(obj)
            % vraci velikosti pole d - samples, channels, epochs
            samples = size(obj.d,1);
            channels = size(obj.d,2);
            epochs = size(obj.d,3);
        end    
              
        function obj = GetHHeader(obj,H)
            %nacte header z promenne H - 25.5.2016
            obj.CH = CHHeader(H);
            [~, ~, obj.els] = obj.CH.ChannelGroups();  
            assert(max(obj.els)<=size(obj.d,2),'nesouhlasi pocet elektrod - spatny header?');
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
            disp(['vyrazeno ' num2str(numel(RjCh)) ' kanalu']); 
        end
        
        function obj = RejectEpochs(obj,RjEpoch)
            %ulozi cisla vyrazenych epoch - kvuli prevodu mezi touto tridou a CHilbert
            obj.RjEpoch = RjEpoch;
            disp(['vyrazeno ' num2str(numel(RjEpoch)) ' epoch']); 
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
                [Kstring Knum] = obj.PsyData.Category(epoch);    %jmeno a cislo kategorie
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
            disp(['rozdeleno na ' num2str(obj.epochs) ' epoch']); 
        end
        
        function [d,psy_rt]= CategoryData(obj, katnum,rt,opak)
            %vraci eegdata epoch ve kterych podnet byl kategorie/podminky=katnum + reakcni casy 
            %Pokud rt>0, vraci epochy serazene podle reakcniho casu 
            %pokud opak>0, vraci jen jedno opakovani obrazku - hlavne kvuli PPA test 
            assert(obj.epochs > 1,'data not yet epoched'); %vyhodi chybu pokud data nejsou epochovana
            if exist('opak','var') && ~isempty(opak)
                epochyopak = obj.PsyData.GetOpakovani();
                iOpak = ismember(epochyopak , opak);
            else
                iOpak = true(obj.epochs,1);                
            end
            iEpochy = [ ismember(cell2mat(obj.epochData(:,2)),katnum) , obj.GetEpochsExclude() , iOpak]; %seznam epoch v ramci kategorie ve sloupci + epochy, ktere nejsou excludovane
            d = obj.d(:,:,all(iEpochy,2)); %epochy z kategorie, ktere nejsou excludovane
            [~,psy_rt,psy_katnum,~] = obj.PsyData.GetResponses();           
            psy_rt = psy_rt(ismember(psy_katnum,katnum),:);
            iEp = iEpochy( iEpochy(:,1)==1,2); %neexcludovane epochy v ramci kategorie
            psy_rt = psy_rt(iEp); %doufam, ze to jsou stejne pocty
            
            if exist('rt','var') && ~isempty(rt) %chci hodnoty serazene podle reakcniho casu               
                [psy_rt, isorted] = sort(psy_rt);
                d = d(:,:,isorted); 
            end   
        end      
        
        function obj = ChangeReference(obj,ref)            
            assert(any(ref=='heb'),'neznama reference, mozne hodnoty: h e b');
            assert(isobject(obj.CH),'Hammer header not loaded');
            H = obj.CH.H; %kopie headeru
            switch ref %jaky typ reference chci
                case 'h'  %headbox          
                    filterSettings.name = 'car'; % options: 'car','bip','nan'
                    filterSettings.chGroups = 'perHeadbox';        % options: 'perHeadbox' (=global CAR), OR 'perElectrode' (=local CAR per el. shank)
                case 'e'  %elektroda
                    filterSettings.name = 'car'; % options: 'car','bip','nan'
                    filterSettings.chGroups = 'perElectrode';        % options: 'perHeadbox' (=global CAR), OR 'perElectrode' (=local CAR per el. shank)
                case 'b'  %bipolarni
                    filterSettings.name = 'bip'; % options: 'car','bip','nan'           
            end
            filterMatrix = createSpatialFilter_kisarg(H, numel(H.selCh_H), filterSettings,obj.RjCh); %ve filterMatrix uz nejsou rejectovane kanaly
                %assert(size(rawData,2) == size(filterMatrix,1));
            % apply spatial filter
            if ref=='b' %u bipolarni reference se mi meni pocet kanalu
                H.channels = struct;
                for ch = 1:size(filterMatrix,2)
                    oldch = find(filterMatrix(:,ch)==1);
                    fnames = fieldnames(obj.CH.H.channels(oldch)); %jmena poli struktury channels
                    for f = 1:numel(fnames); %postupne zkopiruju vsechny pole struktury, najednou nevim jak to udelat
                        fn = fnames{f};
                        H.channels(ch).(fn) = obj.CH.H.channels(oldch).(fn); 
                    end
                    H.channels(ch).name = [H.channels(ch).name '-' obj.CH.H.channels(filterMatrix(:,ch)==-1).name]; %pojmenuju kanal jako rozdil
                end                
            end
            
            if obj.epochs <= 1 %ne epochovana data
                filtData = obj.d(:,H.selCh_H) * filterMatrix;
                assert(size(filtData,1) == size(obj.d,1),'zmenila se delka zaznamu'); %musi zustat stejna delka zaznamu  
                obj.d=filtData;                
            else %epochovana data
                dd = zeros(obj.samples*obj.epochs,numel(H.selCh_H));
                for ch = 1:numel(H.selCh_H) %predelam matici 3D na 2D
                    dd(:,ch) = reshape(obj.d(:,ch,:),obj.samples*obj.epochs,1);
                end                
                filtData = dd(:,H.selCh_H) * filterMatrix;
                assert(size(filtData,1) == size(dd,1),'zmenila se delka zaznamu'); %musi zustat stejna delka zaznamu  
                obj.d = zeros(obj.samples,size(filtData,2),obj.epochs); %nove pole dat s re-referencovanymi daty
                for ch=1:size(filtData,2) %vratim puvodni 3D tvar matice
                    obj.d(:,ch,:) = reshape(filtData(:,ch),obj.samples,obj.epochs);
                end
            end
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
            obj.RjCh = []; %rejectovane kanaly uz byly vyrazeny, ted nejsou zadne
            obj.GetHHeader(H); %novy header s vyrazenymi kanaly
            obj.filename = []; %nechci si omylem prepsat puvodni data 
            switch ref
                case 'h', obj.reference = 'perHeadbox';
                case 'e', obj.reference = 'perElectrode'; 
                case 'b', obj.reference = 'Bipolar';                    
            end
            disp(['reference zmenena: ' obj.reference]); 
        end
        
        function [iEp,epochsEx]=GetEpochsExclude(obj)
            %vraci iEp=index epoch k vyhodnoceni - bez chyb, treningu a rucniho vyrazeni
            %epochsEx=seznam vsech epoch s 1 u tech k vyrazeni
            chyby = obj.PsyData.GetErrorTrials();
            epochsEx = [chyby , zeros(size(chyby,1),1) ]; %pridam dalsi prazdny sloupec
            epochsEx(obj.RjEpoch,5)=1; %rucne vyrazene epochy podle EEG
            iEp = all(epochsEx==0,2); %index epoch k pouziti
        end
            
        function obj = ResponseSearch(obj,timewindow,kats,opakovani)
            %projede vsechny kanaly a hleda signif rozdil proti periode pred podnetem
            %timewindow - pokud dve hodnoty - porovnava prumernou hodnotu mezi nimi - sekundy relativne k podnetu/odpovedi
            % -- pokud jedna hodnota, je to sirka klouzaveho okna - maximalni p z teto delky
            assert(obj.epochs > 1,'only for epoched data');           
            iepochtime = round(obj.epochtime(1:2).*obj.fs); %v poctu vzorku cas pred a po udalosti, pred je zaporne cislo
            if ~ismember(obj,'baseline')
                obj.baseline = [obj.epochtime(1) 0];
            end
            ibaseline = round(obj.baseline.*obj.fs); %zaporna cisla pokud pred synchro eventem
            ibaselineA =  [ibaseline(1) ,    floor((ibaseline(2)-ibaseline(1))/2)+ibaseline(1)]; %prvni pulka baseline - 28.3.2017
            ibaselineB =  [ibaselineA(2)+1 , ibaseline(2) ]; %druha pulka baseline
            itimewindow = round(timewindow.*obj.fs); %
            iEp = obj.GetEpochsExclude(); %ziska seznam epoch k vyhodnoceni
            %proc driv?  mean(obj.d( abs(iepochtime(1)-ibaselineA(1))+1 : abs(iepochtime(1)-ibaselineA(2))  - 28.3.2017            
            baselineA = mean(obj.d( ibaselineA(1)-iepochtime(1)+1 : ibaselineA(2)-iepochtime(1) , : , iEp),1); %#ok<PROP>
            baselineB = mean(obj.d( ibaselineB(1)-iepochtime(1) : ibaselineB(2)-iepochtime(1) , : , iEp),1); %#ok<PROP>
                % cas x kanaly x epochy - prumer za cas pred podnetem, pro vsechny kanaly a nevyrazene epochyt
            if numel(itimewindow) == 2  %chci prumernou hodnotu d od do
                response = mean(obj.d( (itimewindow(1) : itimewindow(2)) - iepochtime(1) , : , iEp ),1); %prumer v case                
            else
                response = obj.d(abs(iepochtime(1)-ibaseline(2))+1 : end , : , iEp ); %hodnoty po konci baseline                 
            end
            WpA = CStat.Wilcox2D(response,baselineA,1,[],'mean vs baseline A'); %#ok<PROP> %1=mene striktni pdep, 2=striktnejsi dep;
            WpB = CStat.Wilcox2D(response,baselineB,1,[],'mean vs baseline B'); %#ok<PROP> %1=mene striktni pdep, 2=striktnejsi dep;
            Wp = max(WpA,WpB); %#ok<PROP> %vyssi hodnota z kazde poloviny baseline
            if numel(itimewindow) == 1 %chci maximalni hodnotu p z casoveho okna
                Wp = CStat.Klouzaveokno(Wp,itimewindow(1),'max',1); %#ok<PROP>
                obj.Wp.D2 = Wp; %#ok<PROP> %pole 2D signifikanci si ulozim kvuli kresleni - cas x channels
                obj.Wp.D2params = timewindow;
                obj.Wp.D2fdr = 1;
                obj.Wp.D2iEp = iEp; %index zpracovanych epoch pro zpetnou kontrolu
            else
                obj.Wp.D1 = Wp; %#ok<PROP>%pole 1D signifikanci - jedna hodnota pro kazdy kanal
                obj.Wp.D1params = timewindow;
                obj.Wp.D1iEp = iEp; %index zpracovanych epoch pro zpetnou kontrolu
            end
            if exist('opakovani','var') && ~isempty(opakovani)
                if iscell(opakovani) 
                    assert(numel(opakovani)<=3,'kategorie opakovani mohou byt maximalne tri');
                    KATNUM = kats;
                    kats = opakovani;   %POZOR kats se meni na opakovani, abych mohl pouzit kod dole             
                end
            end
            
            if exist('kats','var') && numel(kats)>1  && numel(timewindow)==1                                      
                %ziskam eeg data od jednotlivych kategorii
                responsekat = cell(numel(kats),1); %eeg response zvlast pro kazdou kategorii 
                baselinekat = cell(numel(kats),1); %baseline zvlast pro kazdou kategorii 
                for k = 1:numel(kats)
                    if exist('KATNUM','var') 
                        katdata = obj.CategoryData(KATNUM,[],kats{k}); %v kats jsou ted opakovani
                    else
                        katdata = obj.CategoryData(kats(k)); %epochy time*channel*epochs jedne kategorie, uz jsou vyrazeny vyrazene epochy
                    end
                    responsekat{k,1} = katdata( ibaseline(2) - iepochtime(1)+1 :end,:,:); %jen cas po podnetu : cas x channel x epochs; 
                    baselinekat{k,1} = katdata( ibaseline(1) - iepochtime(1)+1 : ibaseline(2) - iepochtime(1),:,:); %jen cas po podnetu : cas x channel x epochs; 
                end
                %rozdily kategorii vuci sobe
                WpKat = cell(numel(kats));
                for k = 1:numel(kats) %budu statisticky porovnavat kazdou kat s kazdou, bez ohledu na poradi
                    for j = k+1:numel(kats)
                        Wp = CStat.Wilcox2D(responsekat{k}, responsekat{j},1,[],['kat ' num2str(k) ' vs ' num2str(j)]); %#ok<PROP> % -------- WILCOX kazda kat s kazdou 
                        WpKat{k,j} = Wp; %#ok<PROP> %CStat.Klouzaveokno(Wp,itimewindow(1),'max',1); %#ok<PROP>
                    end
                end
                %rozdily kategorii vuci baseline - 28.3.2017
                WpKatBaseline = cell(numel(kats),1);
                for k =  1: numel(kats)
                        baselineall = baselinekat{k};
                        baselineA = mean(baselineall(1:floor(size(baselineall,1)/2)      ,:,:));
                        baselineB = mean(baselineall(  floor(size(baselineall,1)/2)+1:end,:,:));
                        WpBA = CStat.Wilcox2D(responsekat{k},baselineA,1,[],['kat ' num2str(k) ' vs baseline A']);
                        WpBB = CStat.Wilcox2D(responsekat{k},baselineB,1,[],['kat ' num2str(k) ' vs baseline B']);
                        WpKatBaseline{k,1} = max (WpBA,WpBB);
                end
                %ukladam vysledky
                obj.Wp.WpKat = WpKat; %rozdily mezi kategorieme
                obj.Wp.WpKatBaseline = WpKatBaseline; %rozdily kategorii vuci baseline
                if exist('KATNUM','var') %pokud vyhodnocuju opakovani
                    obj.Wp.kats = KATNUM;    %puvodni kategorie
                    obj.Wp.opakovani = kats; %v kats jsou ted opakovani
                else
                    obj.Wp.kats = kats; %ulozim si cisla kategorii kvuli grafu PlotResponseCh
                    obj.Wp.opakovani = {}; %opakovani nedelam
                end                    
            else
                obj.Wp.kats = kats;
                obj.Wp.WpKat = cell(0);
                obj.Wp.opakovani = {};
            end
            obj.DatumCas.ResponseSearch = datestr(now);
        end
        
        function Categories(obj)
            %funkce ktera jen vypise kategorie
            obj.PsyData.Categories(1);
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
          
        function [obj]= Decimate(obj,podil)
            %zmensi data na nizsi vzorkovaci frekvenci podle urceneho podilu, naprilkad 4x 512-128Hz, pokud podil=4
            dd = zeros(ceil(size(obj.d,1)/podil) , size(obj.d,2), size(obj.d,3)); % d uz obsahuje jen svoji prvni pulku a delka je delitelna podilem
            fprintf('channels to decimate (z %i):',numel(obj.channels));
            for ch = 1:obj.channels %musim decimovat kazdou elektrodu zvlast
                fprintf('%i, ',ch);
                if obj.epochs == 1
                    dd(:,ch) = decimate(obj.d(:,ch),podil); %na 500 Hz z 8000 Hz
                    if ch==1, obj.tabs = downsample(obj.tabs,podil); end % to delam jen jednou, treba u prvniho kanalu                   
                else
                    for ep = 1:obj.epochs
                        dd(:,ch,ep) = decimate(obj.d(:,ch,ep),podil); %na 500 Hz z 8000 Hz
                        if ch==1, obj.tabs(:,ep) = downsample(obj.tabs(:,ep),podil); end 
                    end
                end                
            end
            obj.tabs_orig = downsample(obj.tabs_orig,podil); % 
            obj.fs = obj.fs/podil;
            obj.d = dd;
            fprintf('... done\n');
            disp(['decimated to ' num2str(obj.fs) ' Hz']);
        end
        
        function [prumery, MNI] = IntervalyResp(obj, intervaly,channels)
            %vypocita hodnoty v jednotlivych intervalech casu pro jednotlive kategorie i pro celkovy prumer       
            %vykresli graf pro kazdy interval do spolecneho plotu
            %vraci prumery[channels x intervaly x kategorie] a MNI(channels)
            assert(isfield(obj.Wp, 'kats'),'musi byt definovany kategorie podnetu');
            assert(isfield(obj.Wp, 'WpKatBaseline'),'musi byt spocitana statisika kategorii');
            if ~exist('channels','var') , channels = 1:obj.channels; end
            kats = obj.Wp.kats; 
            prumery = zeros(numel(channels),size(intervaly,1),1+numel(kats));   % casy x kategorie x channels - celkova data a jednotlive kategorie
            figure('Name','IntervalyResp');
            for j = 1:size(intervaly,1) 
                subplot(2,ceil(size(intervaly,1) /2),j); %pro kazdy interval jiny subplot
                %spocitam prumery celkove i za kazdou kategorii v kazdem casovem intervalu
                % dve cisla v kazdem sloupci - od do ve vterinach   
                iintervalyData = round((intervaly(j,:)-obj.epochtime(1)).*obj.fs); % pro data kde je na zacatku baseline             
                iintervalyStat = round(intervaly(j,:).*obj.fs); % pro statistiku, kde na zacatku neni baseline              
                katdata = obj.CategoryData(kats); 
                iCh = min(obj.Wp.D2(iintervalyStat(1):iintervalyStat(2),channels),[],1) < 0.05; %kanaly kde je signifikantni rozdil vuci baseline, alesponjednou
                prumery(iCh,j,1) = mean(mean(katdata(iintervalyData(1):iintervalyData(2),iCh,:),3),1); %prumer za vsechy epochy a cely casovy interval
                for k = 1: numel(kats)
                    katdata = obj.CategoryData(k); 
                    iCh = min(obj.Wp.WpKatBaseline{k,1}(iintervalyStat(1):iintervalyStat(2),channels),[],1) < 0.05; %kanaly kde je signifikantni rozdil vuci baseline, alespon jednou
                    prumery(iCh,j,1+k) = mean(mean(katdata(iintervalyData(1):iintervalyData(2),iCh,:),3),1);
                    P = squeeze(prumery(:,j,1+k));
                    plot(P','.-'); %kreslim tuto kategorii
                    hold all;
                end
                P = squeeze(prumery(:,j,1)); %nakonec vykreslim prumer vsech katevorii, aby byl nejvic videt
                plot(P','.-');
                legend('kat1','kat2','kat3','mean','Location','NorthWest');
            end 
            MNI = obj.CH.GetMNI(channels);
            assert(numel(MNI)==size(prumery,1),'MNI a prumery maji jiny pocet kanalu');
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
        
        function obj = PlotEpochs(obj,ch,kategories)
            %uchovani stavu grafu, abych ho mohl obnovit a ne kreslit novy
            assert(obj.epochs > 1,'only for epoched data');
            if ~exist('ch','var')
                if isfield(obj.plotEp,'ch'), ch = obj.plotEp.ch;
                else ch = 1; obj.plotEp.ch = ch; end
            else
                obj.plotEp.ch = ch;
            end
            if ~exist('kategories','var')
                if isfield(obj.plotEp,'kategories'), kategories = obj.plotEp.kategories;
                else kategories = obj.PsyData.Categories(); obj.plotEp.kategories = kategories; end
            else
                obj.plotEp.kategories = kategories;
            end
            if isfield(obj.plotEp,'fh') && ishandle(obj.plotEp.fh)
                figure(obj.plotEp.fh); %pouziju uz vytvoreny graf
                %clf(obj.plotEp.fh); %graf vycistim
            else
                obj.plotEp.fh = figure('Name','All Epochs','Position', [20, 100, 1200, 300]);
            end
            T = linspace(obj.epochtime(1),obj.epochtime(2),size(obj.d,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            
            maxy = 0; %budu pocitat za vsechny kategorie
            miny = 0;
            for k=1:numel(kategories)
                katnum = kategories(k);
                subplot(1,numel(kategories),k);
                [dkat,rt] = obj.CategoryData(katnum,1);
                E = 1:size(dkat,3); %cisla epoch - kazdou kategorii muze byt jine                
                D = squeeze(dkat(:,ch,:));
                imagesc(T,E,D');
                maxy = max([maxy max(max( D ))]);
                miny = min([miny min(min( D ))]);                
                xlabel('Time [s]');                
                title(obj.PsyData.CategoryName(katnum));
                hold on; 
                if numel(obj.epochtime)<3 || obj.epochtime(3)==0
                    plot(rt,E,'-k','LineWidth',1); %cara reakcnich casu, nebo podnetu, pokud zarovnano podle reakce      
                else
                    plot(-rt,E,'-k','LineWidth',1); %cara reakcnich casu, nebo podnetu, pokud zarovnano podle reakce      
                end
                plot(zeros(size(E,2),1),E,'-k','LineWidth',1); %cara podnetu
            end    
            if isfield(obj.plotEp,'ylim') && numel(obj.plotEp.ylim)>=2 %nactu nebo ulozim hodnoty y
                miny = obj.plotEp.ylim(1); maxy = obj.plotEp.ylim(2);
            else
                obj.plotEp.ylim = [miny maxy];
            end
            for k=1:numel(kategories)
                subplot(1,numel(kategories),k);
                caxis([miny,maxy]);
                if k == 1, ylabel([ 'Epochs - channel ' num2str(ch)]); end %ylabel jen u prniho obrazku
                if k == numel(kategories), colorbar('Position',[0.92 0.1 0.02 0.82]); end
            end
            methodhandle = @obj.hybejPlotEpochs;
            set(obj.plotEp.fh,'KeyPressFcn',methodhandle); 
        end
        
        function PlotCategory(obj,katnum,channel)
            %vykresli vsechny a prumernou odpoved na kategorii podnetu
            %nahrazeno funkcemi PlotResponseCh a PlotEpochs
            d1=obj.CategoryData(katnum); %epochy jedne kategorie
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
            
            if isempty(obj.plotH)
                obj.plotH = figure('Name','Electrode Plot'); %zatim zadny neni, novy obrazek                 
            else
                figure(obj.plotH);  %kreslim do existujiciho plotu
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
            hold off;
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
                for sj = s:ss %po vsechny zobrazene epochy
                    if find(obj.RjEpoch==sj) 
                        if sj == s, titul = [titul ' - EXCLUDED'];  end %#ok<AGROW>
                        line([obj.epochtime(1) obj.epochtime(2)]+(sj-s)*(obj.epochtime(2)-obj.epochtime(1)),[shift(1,1) shift(end,1)],'Color','r','LineWidth',2);
                    end
                    if find(obj.epochTags==sj)
                        if sj == s, titul = [titul ' - TAGGED']; end   %#ok<AGROW>
                        line([0 0]+(sj-s)*(obj.epochtime(2)-obj.epochtime(1)),[shift(1,1) shift(end,1)],'Color','g','LineWidth',4);
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
            
            for k= 1 : size(els,2) %#ok<PROP> %
                    uistack(h_els{k}, 'top'); %dam krivky eeg ulne dopredu
            end
            
        end
        
        function obj = PlotResponses(obj)
            %vykresli uspesnost odpovedi spolu s chybami a vyrazenymi epochami
            obj.PsyData.PlotResponses();
            figure(obj.PsyData.fhR); %kreslim dal do stejneho obrazku
            [~,rt,kategorie] = obj.PsyData.GetResponses();            
            plot(obj.RjEpoch,rt(obj.RjEpoch),'*','MarkerSize',10,'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]); %vykreslim vyrazene epochy
            plot(obj.RjEpoch,kategorie(obj.RjEpoch),'*r','MarkerSize',5); %vykreslim vyrazene epochy
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
                if isfield(obj.plotRCh,'ch'), ch = obj.plotRCh.ch;
                else ch = 1; obj.plotRCh.ch = ch; end
            else
                obj.plotRCh.ch = ch;
            end
            if ~exist('kategories','var') || isempty(kategories)
                if isfield(obj.Wp, 'kats')
                    kategories = obj.Wp.kats; %pokud jsou kategorie v parametru, prvni volba je pouzit je ze statistiky
                elseif isfield(obj.plotRCh,'kategories')
                    kategories = obj.plotRCh.kategories; %hodnoty drive pouzite v grafu, ty maji prednost pred statistikou
                elseif isfield(obj.Wp,'kats')
                    kategories = obj.Wp.kats; %hodnoty pouzite ve statistice
                else
                   if numel(obj.PsyData.Categories())<=3
                     kategories = obj.PsyData.Categories(); %pokud neni vic nez 3 kategorie, vezmu vsechny
                     obj.plotRCh.kategories = kategories;
                   end %pokud je kategorii vic nez tri, neberu je v uvahu a zobrazim pouze prumer
                end                
            else                
                assert(numel(kategories)<=3,'kategorie mohou byt maximalne tri');
                if ~isempty(obj.Wp.WpKat) && (isempty(obj.Wp.kats) || ~isequal(obj.Wp.kats, kategories))
                    disp('Statistika spocitana bez kategorii nebo pro jine kategorie')
                end
                obj.plotRCh.kategories = kategories;    %hodnoty zadane parametrem, ty maji absolutni prednost
            end
            %opakovani obrazku kvuli PPA - 28.9.2016
            if ~exist('opakovani','var') || isempty(opakovani)     
                if isfield(obj.Wp,'opakovani')
                    opakovani = obj.Wp.opakovani;
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
                if ~isempty(obj.Wp.WpKat) && (isempty(obj.Wp.opakovani) || ~isequal(obj.Wp.opakovani, opakovani))
                    disp('Statistika spocitana bez opakovani nebo pro jina opakovani')
                end
                obj.plotRCh.opakovani = opakovani;    %hodnoty zadane parametrem, ty maji absolutni prednost
            end
            KATNUM = kategories; % kategorie, ktere chci vykreslovat - vsechny dohromady     
            if ~isempty(opakovani)                
                kategories = opakovani; %POZOR - misto kategorii jsou nyni opakovane - cell array
            end 
            T = linspace(obj.epochtime(1),obj.epochtime(2),size(obj.d,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            if isfield(obj.plotRCh,'fh')
                figure(obj.plotRCh.fh); %pouziju uz vytvoreny graf - tohle nefunguje v matlab 2016b
                clf(obj.plotRCh.fh); %graf vycistim
            else
                obj.plotRCh.fh = figure('Name','W plot channel');
            end
            [ymin ymax] = obj.responseChYLim(iff(~isempty(opakovani),KATNUM,kategories));
            
            %ZACINAM VYKRESLOVAT - NEJDRIV MEAN VSECH KATEGORII
            if ~exist('kategories','var') && ~exist('opakovani','var') %26.5.2017 - jen kdyz neexistuji kategorie
                katdata =  obj.CategoryData(KATNUM);
                M = mean(katdata(:,ch,:),3);             
                E = std(katdata(:,ch,:),[],3)/sqrt(size(katdata,3)); %std err of mean          
                h_errbar = errorbar(T,M,E,'.','Color',[.6 .6 1]); %nejdriv vykreslim errorbars aby byly vzadu [.8 .8 .8]
                hold on;
                h_mean = plot(T,M,'LineWidth',2,'Color',[0 0 1]);  %prumerna odpoved, ulozim si handle na krivku          
                xlim(obj.epochtime(1:2));
               
                obj.plotRCh.range = [min(M)-max(E) max(M)+max(E)]; %zjistim a ulozim rozsah hodnot pro moznost nastaveni osy y
                if isfield(obj.Wp,'D2') %krivka p hodnot z W testu
                    Tr = linspace(0,obj.epochtime(2),size(obj.Wp.D2,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
                    if pvalue %pokud chci zobrazovat hodnotu p value jako krivku
                        plot(Tr,obj.Wp.D2(:,ch),'b:');  %carkovana modra cara oznacuje signifikanci prumeru
                    end
                    y = ymin + (ymax-ymin)*0.2;
                    iWp = obj.Wp.D2(:,ch) <= 0.05;
                    plot(Tr(iWp),ones(1,sum(iWp))*y,'b.'); %tecky jsou p < 0.05                
                    iWpfirst = find(iWp,1,'first');                 
                    if(numel(iWpfirst)>0) 
                        text(-0.01,y,[ num2str( round(Tr(iWpfirst)*1000)) 'ms']); %cas zacatku signifikance
                        text(-0.18,y,[ 'p=' num2str(CStat.round(min(obj.Wp.D2(:,ch)),3))]);  %cas zacatku signifikance 
                        line([Tr(iWpfirst) Tr(iWpfirst)],obj.plotRCh.ylim,'Color','blue'); %modra svisla cara u zacatku signifikance
                        text(0.05,y, 'Mean');
                    end
                    iWp = obj.Wp.D2(:,ch) <= 0.01;
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
                colorskat = {[0 0 0],[0 1 0],[1 0 0]; [hue hue hue],[hue 1 hue],[1 hue hue]}; % prvni radka - prumery, druha radka errorbars = svetlejsi
                h_kat = zeros(numel(kategories),1); 
                for k= 1 : numel(kategories) %index 1-3
                    if exist('opakovani','var') && ~isempty(opakovani)
                        opaknum = kategories{k}; %v kategories jsou opakovani k vykresleni, a je to cell array
                        katdata = obj.CategoryData(KATNUM,[],opaknum); %eegdata - epochy pro tato opakovani
                    else
                        katnum = kategories(k); %cislo kategorie
                        katdata = obj.CategoryData(katnum); %eegdata - epochy jedne kategorie
                    end    
                    M = mean(katdata(:,ch,:),3);
                    E = std(katdata(:,ch,:),[],3)/sqrt(size(katdata,3)); %std err of mean
                    errorbar(T,M,E,'.','color',colorskat{2,k}); %nejdriv vykreslim errorbars aby byly vzadu[.8 .8 .8]
                    xlim(obj.epochtime(1:2)); 
                    hold on;
                    h_kat(k) = plot(T,M,'LineWidth',katlinewidth,'Color',colorskat{1,k});  %prumerna odpoved,  ulozim si handle na krivku  
                    obj.plotRCh.range = [ min(obj.plotRCh.range(1),min(M)-max(E)) max(obj.plotRCh.range(2),max(M)+max(E))]; %pouziju to pak pri stlaceni / z obrazku
                    Tr = linspace(0,obj.epochtime(2),size(obj.Wp.D2,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
                    if isfield(obj.Wp,'WpKat') 
                        for l = k+1:numel(kategories) %katnum jde od nuly
                            y = ymin + (ymax-ymin)*(0.3 - (k+l)*0.05)  ; %pozice na ose y
                            if k==1, color=colorskat{1,l}; else color = colorskat{1,1}; end %green a red jsou proti kategorii 0, cerna je kat 1 vs kat 2
                            if pvalue %pokud chci zobrazovat hodnotu p value jako krivku                                
                                plot(Tr,obj.Wp.WpKat{k,l}(:,ch), ':','Color',color); %carkovana cara oznacuje signifikanci kategorie vuci jine kategorii
                            end
                            %nejdriv p < 0.05                            
                            iWp = obj.Wp.WpKat{k,l}(:,ch)  <= 0.05; 
                            plot(Tr(iWp),ones(1,sum(iWp))*y, '.','Color',color); %                        
                            iWpfirst = find(iWp,1,'first');                        
                            if(numel(iWpfirst)>0)                                
                                text(-0.01,y,[ num2str(round(Tr(iWpfirst)*1000)) 'ms']);  %cas zacatku signifikance 
                                text(-0.18,y,[ 'p=' num2str(CStat.round(min(obj.Wp.WpKat{k,l}(:,ch)),3))]);  %cas zacatku signifikance 
                                line([Tr(iWpfirst) Tr(iWpfirst)],obj.plotRCh.ylim,'Color',color); %modra svisla cara u zacatku signifikance                                
                            end                            
                            %potom jeste p < 0.01
                            iWp = obj.Wp.WpKat{k,l}(:,ch)  <= 0.01;                          
                            plot(Tr(iWp),ones(1,sum(iWp))*y,  '*','Color',color); %
                            % jmena kategorii vypisuju vzdy
                            if exist('opakovani','var') && ~isempty(opakovani)   %pokud vyhodnocuju opakovani
                                kat1name =  obj.PsyData.OpakovaniName(kategories{l});
                                kat2name =  obj.PsyData.OpakovaniName(kategories{k});
                                kat3name =  [ ' (' obj.PsyData.CategoryName(obj.Wp.kats) ')' ]; %jmeno kategorie obrazku, ze ktere se opakovani pocitalo
                            else
                                kat1name =  obj.PsyData.CategoryName(kategories(l));
                                kat2name =  obj.PsyData.CategoryName(kategories(k));
                                kat3name = '';
                            end
                            text(0.05,y, ['\color[rgb]{' num2str(colorskat{1,l}) '}' kat1name ...
                                    '\color[rgb]{' num2str(color) '} *X* '  ...
                                    '\color[rgb]{' num2str(colorskat{1,k}) '}' kat2name kat3name]);                            
                        end                      
                        
                    end
                    if isfield(obj.Wp,'WpKatBaseline') %signifikance vuci baseline
                            iWpB = obj.Wp.WpKatBaseline{k,1}(:,ch)  <= 0.05; %nizsi signifikance
                            y = ymin + (ymax-ymin)*(0.28 - (k+2)*0.05)  ;
                            plot(Tr(iWpB),ones(1,sum(iWpB))*y, '.','Color',colorskat{1,k},'MarkerSize',5); % 
                            iWpB = obj.Wp.WpKatBaseline{k,1}(:,ch)  <= 0.01; % vyssi signifikance
                            %y = ymin + (ymax-ymin)*(0.28 - (k+2)*0.05)  ;
                            plot(Tr(iWpB),ones(1,sum(iWpB))*y, 'p','Color',colorskat{1,k},'MarkerSize',5); % 
                            if exist('opakovani','var') && ~isempty(opakovani)
                                kat2name =  obj.PsyData.OpakovaniName(kategories{k}); %pokud vyhodnocuju opakovani
                            else
                                kat2name =  obj.PsyData.CategoryName(kategories(k));
                            end
                            text(0.05, y, ['\color[rgb]{' num2str(colorskat{1,k}) '}' kat2name ' vs.baseline'] );
                            line([0 1],[y y]+(ymax-ymin)*0.03 ,'Color',[0.5 0.5 0.5]);
                                %kazde jmeno kategorie jinou barvou
                            if pvalue %pokud chci zobrazovat hodnotu p value jako krivku
                               plot(Tr,obj.Wp.WpKatBaseline{k,1}(:,ch), '.','Color',colorskat{1,k}); %tesckovana cara oznacuje signifikanci kategorie vuci baseline
                            end
                     end
                end
                for k= 1 : numel(kategories) %index 1-3
                    uistack(h_kat(k), 'top'); %dam krivky prumeru kategorii pred jejich errorbars
                end
                
                ylim( [ymin ymax].*1.1);
            end
            if ~isempty(h_mean)
                uistack(h_errbar, 'top');
                uistack(h_mean, 'top'); %uplne nahoru dam prumer vsech kategorii
            end            
            title(['channel ' num2str(ch) ' - ' obj.PsyData.PacientID()]); % v titulu obrazku bude i pacientID napriklad p132-VT18
            text(-0.1,ymax*.95,[ obj.CH.H.channels(1,ch).name ' : ' obj.CH.H.channels(1,ch).neurologyLabel ',' obj.CH.H.channels(1,ch).ass_brainAtlas]);
            if  isfield(obj.CH.H.channels,'MNI_x') %vypisu MNI souradnice
                text(-0.1,ymax*.90,[ 'MNI:' num2str(obj.CH.H.channels(1,ch).MNI_x) ',' num2str(obj.CH.H.channels(1,ch).MNI_y ) ',' num2str(obj.CH.H.channels(1,ch).MNI_z)]);
            else
                text(-0.1,ymax*.90,'no MNI');
            end
            methodhandle = @obj.hybejPlotCh;
            set(obj.plotRCh.fh,'KeyPressFcn',methodhandle);          
        end        
            
        function obj = PlotResponseP(obj)
            %vykresli signifikanci odpovedi u vsech kanalu EEG vypocitanou pomoci ResponseSearch                                
            assert(obj.epochs > 1,'only for epoched data');
            T = 0:0.1:obj.epochtime(2); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            if isfield(obj.Wp,'D1') %první 2D plot
                figure('Name','W plot 1D');
                isignif = obj.Wp.D1<0.05;
                plot(find(~isignif),obj.Wp.D1(~isignif),'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','b');
                hold on;
                plot(find(isignif),obj.Wp.D1(isignif),'o','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r');
                ylim([0 0.1]);
                view(-90, 90); %# Swap the axes
                set(gca, 'ydir', 'reverse'); %# Reverse the y-axis 
                set(gca, 'xdir', 'reverse'); %# Reverse the x-axis 
                for e = 1:numel(obj.els) %hranice elektrod a jmeno posledniho kontaktu
                    line([obj.els(e)+0.5 obj.els(e)+0.5],[0 0.1],'color',[.5 0.5 0.5]);
                    text(obj.els(e)-1,-0.01,obj.CH.H.channels(1,obj.els(e)).name);
                end
                for ch=1:obj.channels
                    if obj.Wp.D1(ch)<0.1 %anatomicka jmena u signif kontaktu
                        text(ch,0.102,obj.CH.H.channels(1,ch).neurologyLabel);
                    end
                end
            end
            if isfield(obj.Wp,'D2') %isprop(obj,'Wp') && isfield(obj.Wp,'D2')
                figure('Name','W map 2D');
                imagesc(T,1:obj.channels,1 - obj.Wp.D2', [0.95 1]); %mapa, od p>0.05 bude modra barva 
                axis ij;
                ylabel('channels');
                xlabel('time [s]');
                colorbar;
                for e = 1:numel(obj.els) %hranice elektrod a jmeno posledniho kontaktu
                    line([T(1) T(end)],[obj.els(e)+0.5 obj.els(e)+0.5],'color','w');
                    text(-T(end)/10,obj.els(e)-1,obj.CH.H.channels(1,obj.els(e)).name);
                end
                for ch=1:obj.channels %anatomicka jmena u signif kontaktu
                    if any(obj.Wp.D2(:,ch)<0.05)
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
            assert(isobject(obj.CH),'Hammer header not loaded');
            evts = zeros(numel(obj.CH.H.selCh_H),1); %pocet epievents pro kazdy kanal
            names = cell(numel(obj.CH.H.selCh_H));  %jmena kanalu
            epochs = evts; %#ok<PROP>
            namelast = ''; %tam budu ukladat pismeno jmena elektrody
            evts_nonseeg = 0; %pocet epi udalost v nonseeg kanalech
            obj.DE.Clear_iDEtabs();
            if obj.epochs > 1 %epochovana data
                tabs = [ obj.tabs(1,:)' obj.tabs(end,:)' ]; %#ok<PROP> %zacatky a konce vsech epoch
            end 
            for ch = 1:obj.channels; %obj.CH.H.selCh_H                
                if obj.epochs==1
                    [epitime ~] = obj.DE.GetEvents([obj.tabs(1) obj.tabs(end)],ch);
                else %epochovana data                   
                    [epitime ~] = obj.DE.GetEvents(tabs,ch,obj.tabs_orig(1)); %#ok<PROP>                     
                    if numel(epitime) > 0 && ~isempty(find(obj.CH.H.selCh_H==ch, 1))
                        epochs(ch) = numel( unique(epitime(:,2))); %#ok<PROP>                        
                    end
                end
                if ~isempty(find(obj.CH.H.selCh_H==ch, 1)) %pokud je kanal ve vyjmenovanych SEEG kanalech podle headeru
                    evts(ch) = size(epitime,1); 
                    if strcmp(namelast,obj.CH.H.channels(ch).name(1))
                       names{ch} = obj.CH.H.channels(ch).name(end); %cislo elektrody bez pismene, u druhe elektrody stejneho jmena abych usetril misto
                    else
                       names{ch} = obj.CH.H.channels(ch).name(1);
                       namelast = obj.CH.H.channels(ch).name(1);                     
                    end
                else
                   evts_nonseeg = evts_nonseeg + size(epitime,1);
                end
                
            end
            figure('Name','Epievents in individual channels');
            if obj.epochs > 1
                subplot(2,1,1);                
            end
            
            plot(evts,'.-');
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
                    line([obj.els(el) obj.els(el)]+1,[0 max(evts)],'Color',[0.5 0.5 0.5]);
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
                PsyDataP = obj.PsyData.P;       %#ok<NASGU>         %ulozim pouze strukturu 
            else
                PsyDataP = []; %#ok<NASGU>
            end
            epochtime = obj.epochtime;      %#ok<PROP,NASGU>
            baseline = obj.baseline;        %#ok<PROP,NASGU>
            CH_H=obj.CH.H;                  %#ok<NASGU>            
            els = obj.els;                  %#ok<PROP,NASGU>
            plotES = obj.plotES;            %#ok<PROP,NASGU>
            %plotH = obj.plotH;              %#ok<PROP,NASGU> %plotH je blbost ukladat, vytvori se novy, jen to brani vice grafum - 14.6.2016
            RjCh = obj.RjCh;                %#ok<PROP,NASGU>
            RjEpoch = obj.RjEpoch;          %#ok<PROP,NASGU>
            epochTags = obj.epochTags;      %#ok<PROP,NASGU>
            epochLast = obj.epochLast;      %#ok<PROP,NASGU>
            reference = obj.reference;      %#ok<PROP,NASGU>
            epochData = obj.epochData;      %#ok<PROP,NASGU>
            Wp = obj.Wp;                    %#ok<PROP,NASGU>
            DE = obj.DE;                    %#ok<PROP,NASGU>
            DatumCas = obj.DatumCas;        %#ok<PROP,NASGU>
            save(filename,'d','tabs','tabs_orig','fs','header','sce','PsyDataP','epochtime','baseline','CH_H','els',...
                    'plotES','RjCh','RjEpoch','epochTags','epochLast','reference','epochData','Wp','DE','DatumCas','-v7.3');  
            disp(['ulozeno do ' filename]); 
        end
        function obj = Load(obj,filename)
            % nacte veskere promenne tridy ze souboru
            assert(exist(filename,'file')==2, 'soubor s daty neexistuje, nejde o data tridy CHilbert?');
            vars = whos('-file',filename) ;
            assert(ismember('d', {vars.name}), 'soubor neobsahuje promennou d, nejde o data tridy CHilbert?'); 
            load(filename,'d','tabs','tabs_orig','fs','header','sce','epochtime','els','plotES','RjCh','RjEpoch','epochTags','epochLast','reference');            
            obj.d = d;                      %#ok<CPROP,PROP>
            obj.tabs = tabs;                %#ok<CPROP,PROP>
            obj.tabs_orig = tabs_orig;      %#ok<CPROP,PROP>
            obj.fs = fs;                    %#ok<CPROP,PROP>           
            obj.mults = ones(1,size(d,2));       %#ok<CPROP,PROP>
            obj.header = header;            %#ok<CPROP,PROP>
            obj.samples = sce(1); obj.channels=sce(2); obj.epochs = sce(3); 
            vars = whos('-file',filename);
            if ismember('PsyDataP', {vars.name})
                load(filename,'PsyDataP'); obj.PsyData = CPsyData(PsyDataP);%  %vytvorim objekt psydata ze struktury
            else
                load(filename,'PsyData');  obj.PsyData = PsyData ; %#ok<CPROP,PROP> %  %drive ulozeny objekt, nez jsem zavedl ukladani struct
            end
            if obj.epochs > 1
                if ismember('epochData', {vars.name}), load(filename,'epochData');  obj.epochData = epochData;   end   %#ok<CPROP,PROP> 
                if ismember('baseline',  {vars.name}), load(filename,'baseline');   obj.baseline = baseline;   end   %#ok<CPROP,PROP>     
                load(filename,'epochtime');                
                obj.epochtime = epochtime;      %#ok<CPROP,PROP>               
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
                obj.CH = CH; %#ok<CPROP,PROP> %  %drive ulozeny objekt, nez jsem zavedl ukladani struct
            end 
            
            if ismember('Wp', {vars.name})
                load(filename,'Wp');      obj.Wp = Wp; %#ok<CPROP,PROP>
            else
                obj.Wp = struct;
            end
            if ismember('DE', {vars.name}) %1.9.2016
                load(filename,'DE');      obj.DE = DE; %#ok<CPROP,PROP>
            else
                obj.Wp = struct;
            end
            if ismember('DatumCas', {vars.name}) %7.4.2017
                load(filename,'DatumCas');      obj.DatumCas = DatumCas; %#ok<CPROP,PROP>
            else
                obj.DatumCas = {};
            end
            obj.els = els;                  %#ok<CPROP,PROP>
            obj.plotES = plotES;            %#ok<CPROP,PROP>
            %obj.plotH = plotH;              %#ok<CPROP,PROP>
            obj.RjCh = RjCh;                %#ok<CPROP,PROP>     
            obj.RjEpoch = RjEpoch;          %#ok<CPROP,PROP>
            if exist('epochTags','var'),  obj.epochTags = epochTags;   else obj.epochTags = []; end         %#ok<CPROP,PROP>     
            if exist('epochLast','var'),  obj.epochLast = epochLast;   else obj.epochLast = []; end         %#ok<CPROP,PROP>
            if exist('reference','var'),  obj.reference = reference;   else obj.reference = 'original'; end         %#ok<CPROP,PROP> %14.6.2016
           
            obj.filename = filename;
            disp(['nacten soubor ' filename]); 
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
               case 'numpad4' %predchozi oznacena epocha
                   s = obj.plotES(2);
                   prevTag = obj.epochTags(obj.epochTags < s);
                   prevDel = obj.RjEpoch(obj.RjEpoch < s);
                   if numel(prevTag) > 0 || numel(prevDel)>0
                     obj.PlotElectrode(obj.plotES(1),max([prevTag prevDel]),obj.plotES(3),obj.plotES(4));
                   end                   
               case 'numpad6' %dalsi oznacena epocha
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
               case 'rightarrow' 
                   obj.PlotResponseCh( min( [obj.plotRCh.ch + 1 , obj.channels]));
               case 'pagedown' 
                   obj.PlotResponseCh( min( [obj.plotRCh.ch + 10 , obj.channels]));
               case 'leftarrow'
                   obj.PlotResponseCh( max( [obj.plotRCh.ch - 1 , 1]));
               case 'pageup'
                   obj.PlotResponseCh( max( [obj.plotRCh.ch - 10 , 1]));
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
               case {'divide','slash'} %lomeno na numericke klavesnici
                   obj.plotRCh.ylim = obj.plotRCh.range; %spocitalo se pri volani PlotResponseCh
                   obj.PlotResponseCh( obj.plotRCh.ch); %prekreslim grafy
               case 'space' %zobrazi i prumerne krivky
                   if isa(obj,'CHilbert'), obj.PlotResponseFreq(obj.plotRCh.ch,obj.Wp.kats); end %vykreslim vsechna frekvencni pasma
                   obj.PlotEpochs(obj.plotRCh.ch,obj.Wp.kats); %vykreslim prumery freq u vsech epoch
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
               case 'delete' %Del na numericke klavesnici
                   obj.plotEp.ylim = [];
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
            if isfield(obj.plotRCh,'ylim') && numel(obj.plotRCh.ylim)>=2 %pokud mam drive ulozene ylim
                    ylim( obj.plotRCh.ylim .*1.1); %udelam rozsah y o 10% vetsi
                    ymax = obj.plotRCh.ylim(2);
                    ymin = obj.plotRCh.ylim(1);
                else
                    ymax = 0; ymin = 0;
                    for katnum=kategories
                        katdata = obj.CategoryData(katnum); %epochy jedne kategorie
                        channels = 1:obj.channels; %#ok<PROP>
                        channels(ismember(channels, [obj.RjCh obj.CH.GetTriggerCh()]))=[]; %#ok<PROP> %vymazu rejectovana a triggerovane channels 
                        ymax = max([ ymax max(mean(katdata(:,channels,:),3))]); %#ok<PROP>
                        ymin = min([ ymin min(mean(katdata(:,channels,:),3))]); %#ok<PROP>
                    end  
                    ymin = ymin - 0.15*(ymax-ymin); %pridam patnact procent na napisy dole
                    %ylim( [ymin ymax].*1.1); %udelam rozsah y o 10% vetsi
                    obj.plotRCh.ylim = [ymin ymax];
            end
        end
    end
    
end

