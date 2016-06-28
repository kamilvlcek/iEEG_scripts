classdef CiEEGData < handle
    %CEEGDATA Trida na praci s datama ve formatu ISARG od Petra Jezdika
    %   Kamil Vlcek, FGU AVCR, since 2016 04
    
    properties (Access = public)
        d; %double nebo int matrix: time x channel, muze byt i time x channel x epoch
        tabs; 
        tabs_orig; %originalni tabs, ktere se zachovaji po epochaci. Downsamplovani se u nich dela
        fs; %vzorkovaci frekvence
        mults; %nepovinne
        header; %nepovinne
        
        samples; %pocet vzorku v zaznamu = rozmer v case
        channels; %pocet kanalu v zaznamu
        epochs;   %pocet epoch
        epochData; %cell array informaci o epochach; epochy v radcich, sloupce: kategorie, tab
        PsyData; %objekt ve formatu CPsyData (PPA, AEDist aj) podle prezentace KISARG
            %pole PsyData.P.data, sloupce, strings, interval, eegfile, pacientid
        epochtime; %delka eventu pre a po event v sekundach    
        CH; %objekt formatu CHHeader s Hammer headerem 
        els; %cisla poslednich kanalu v kazde elektrode
        plotES; % current electrode, second of plot/epoch, range of y values, time range, allels, rangey all els
        plotH;  % handle to plot
        plotRCh = struct; %stavove udaje o grafu PlotResponseCh
        RjCh; %seznam cisel rejectovanych kanalu
        RjEpoch; %seznam vyrazenych epoch
        epochTags; %seznam oznacenych epoch
        epochLast; %nejvyssi navstivena epocha
        filename;
        reference; %slovni popis reference original, avg perHeadbox, perElectrode, Bipolar
        yrange = [10 10 50 50]; %minimum y, krok y0, hranice y1, krok y1, viz funkce - a + v hybejPlot
        Wp = {}; %pole signifikanci pro jednotlive kanaly vuci baseline, vysledek  ResponseSearch        
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
                assert(size(mults,1)<= 1 || size(mults,1)==size(d,2),'d and mults have to have same number of channels');                
               
                obj.tabs = tabs;
                obj.tabs_orig = tabs;
                obj.fs = fs;
                if exist('mults','var') && ~isempty(mults)
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
            end
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
            [~, obj.els] = obj.CH.ChannelGroups();            
        end
         
        function obj = RejectChannels(obj,RjCh)
            %ulozi cisla vyrazenych kanalu - kvuli pocitani bipolarni reference 
            obj.RjCh = RjCh;
        end
        
        function obj = RejectEpochs(obj,RjEpoch)
            %ulozi cisla vyrazenych epoch - kvuli prevodu mezi touto tridou a CHilbert
            obj.RjEpoch = RjEpoch;
        end
        
        function ExtractEpochs(obj, psy,epochtime)
            % epochuje data v poli d, pridava do objektu:
            % cell array epochData, double(2) epochtime v sekundach, struct psy na tvorbu PsyData
            % upravuje obj.mults, samples channels epochs
            if obj.epochs > 1
                disp('already epoched data');
                return;
            end
            obj.PsyData = CPsyData(psy); %vytvorim objekt CPsyData
            obj.epochtime = epochtime; %v sekundach cas pred a po udalosti  , prvni cislo je zaporne druhe kladne
            iepochtime = round(epochtime.*obj.fs); %v poctu vzorku cas pred a po udalosti
            ts_podnety = obj.PsyData.TimeStimuli(); %timestampy vsech podnetu
            de = zeros(iepochtime(2)-iepochtime(1), size(obj.d,2), size(ts_podnety,1)); %nova epochovana data time x channel x epoch            
            tabs = zeros(iepochtime(2)-iepochtime(1),size(ts_podnety,1)); %#ok<PROP> %udelam epochovane tabs
            obj.epochData = cell(size(ts_podnety,1),3); % sloupce kategorie, cislo kategorie, timestamp
            for epoch = 1:size(ts_podnety,1) %pro vsechny eventy
                izacatek = find(obj.tabs<=ts_podnety(epoch), 1, 'last' ); %najdu index podnetu podle jeho timestampu
                    %kvuli downsamplovani Hilberta, kdy se mi muze ztratit presny cas zacatku
                    %epochy, beru posledni nizsi tabs nez je cas zacatku epochy
                [Kstring Knum] = obj.PsyData.Category(epoch);    %jmeno a cislo kategorie
                obj.epochData(epoch,:)= {Kstring Knum obj.tabs(izacatek)}; %zacatek epochy beru z tabs aby sedel na tabs pri downsamplovani
                for ch = 1:obj.channels %pro vsechny kanaly                    
                    baseline = mean(obj.d(izacatek+iepochtime(1):izacatek-1,ch)); %baseline toho jednoho kanalu, jedne epochy
                    de(:,ch,epoch) = obj.d( izacatek+iepochtime(1) : izacatek+iepochtime(2)-1,ch) - baseline; 
                    tabs(:,epoch) = obj.tabs(izacatek+iepochtime(1) : izacatek+iepochtime(2)-1); %#ok<PROP>
                end
            end
            obj.d = de; %puvodni neepochovana budou epochovana            
            obj.tabs = tabs; %#ok<PROP>
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
        end
        
        function [d]= CategoryData(obj, katnum)
            %vraci epochy ve kterych podnet byl kategorie/podminky katnum
            assert(obj.epochs > 1,'data not yet epoched'); %vyhodi chybu pokud data nejsou epochovana
            iEpochy = cell2mat(obj.epochData(:,2))==katnum ; %seznam epoch v ramci kategorie ve sloupci
            iEp=obj.GetEpochsExclude();            
            d = obj.d(:,:,iEpochy & iEp); %epochy z kategorie, ktere nejsou excludovane
        end      
        
        function ChangeReference(obj,ref)            
            assert(any(ref=='heb'),'neznama reference, mozne hodnoty h e b');
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
            
            if obj.epochs <= 1
                filtData = obj.d(:,H.selCh_H) * filterMatrix;
                assert(size(filtData,1) == size(obj.d,1)); %musi zustat stejna delka zaznamu  
                obj.d=filtData;                
            else
                dd = zeros(obj.samples*obj.epochs,numel(H.selCh_H));
                for ch = 1:numel(H.selCh_H) %predelam matici 3D na 2D
                    dd(:,ch) = reshape(obj.d(:,ch,:),obj.samples*obj.epochs,1);
                end                
                filtData = dd(:,H.selCh_H) * filterMatrix;
                assert(size(filtData,1) == size(dd,1)); %musi zustat stejna delka zaznamu  
                obj.d = zeros(obj.samples,size(filtData,2),obj.epochs); %nove pole dat s re-referencovanymi daty
                for ch=1:size(filtData,2)
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
        end
        
        function [iEp,epochsEx]=GetEpochsExclude(obj)
            %vraci iEp=index epoch k vyhodnoceni - bez chyb, treningu a rucniho vyrazeni
            %epochsEx=seznam vsech epoch s 1 u tech k vyrazeni
            chyby = obj.PsyData.GetErrorTrials();
            epochsEx = [chyby , zeros(size(chyby,1),1) ]; %pridam dalsi prazdny sloupec
            epochsEx(obj.RjEpoch,4)=1; %rucne vyrazene epochy podle EEG
            iEp = all(epochsEx==0,2); %index epoch k pouziti
        end
            
        function ResponseSearch(obj,timewindow,kats)
            %projede vsechny kanaly a hleda signif rozdil proti periode pred podnetem
            %timewindow - pokud dve hodnoty - porovnava prumernou hodnotu mezi nimi
            % -- pokud jedna hodnota, je to sirka klouzaveho okna - maximalni p z teto delky
            assert(obj.epochs > 1,'only for epoched data');
            iepochtime = round(obj.epochtime.*obj.fs); %v poctu vzorku cas pred a po udalosti, pred je zaporne cislo
            itimewindow = round(timewindow.*obj.fs); %
            iEp = obj.GetEpochsExclude(); %ziska seznam epoch k vyhodnoceni
            baseline = mean(obj.d(1: -iepochtime(1),:,iEp),1);  % 1 cas x kanaly x epochy - prumer za cas pred podnete,, pro vsechny epochy a jeden kanal
            if numel(itimewindow) == 2  %chci prumernou hodnotu d od do
                response = mean(obj.d( (itimewindow(1) : itimewindow(2)) - iepochtime(1) , : , iEp ),1); %prumer v case                
            else
                response = obj.d(- iepochtime(1):end , :, iEp ); %hodnoty po podnetu                  
            end
            Wp = CStat.Wilcox2D(response,baseline,1); %#ok<PROP>
            if numel(itimewindow) == 1 %chci maximalni hodnotu p z casoveho okna
                Wp = CStat.Klouzaveokno(Wp,itimewindow(1),'max',1); %#ok<PROP>
                obj.Wp.D2 = Wp; %#ok<PROP> %pole 2D signifikanci si ulozim kvuli kresleni - cas x channels
                obj.Wp.D2params = timewindow;
                obj.Wp.D2iEp = iEp; %index zpracovanych epoch pro zpetnou kontrolu
            else
                obj.Wp.D1 = Wp; %#ok<PROP>%pole 1D signifikanci - jedna hodnota pro kazdy kanal
                obj.Wp.D1params = timewindow;
                obj.Wp.D1iEp = iEp; %index zpracovanych epoch pro zpetnou kontrolu
            end
            if exist('kats','var') && numel(kats)>1  && numel(timewindow)==1                                      
                responsekat = cell(numel(kats),1); %response zvlast pro kazdou kat
                for k = 1:numel(kats)
                    katdata = obj.CategoryData(k-1); %epochy jedne kategorie, uz jsou vyrazeny vyrazene epochy, Kategorie se pocitaji od 0
                    responsekat{k,1} = katdata(- iepochtime(1):end,:,:); %jen cas po podnetu; 
                end
                WpKat = cell(numel(kats));
                for k = 1:numel(kats) %budu statisticky porovnavat kazdou kat s kazdou, bez ohledu na poradi
                    for j = k+1:numel(kats)
                        Wp = CStat.Wilcox2D(responsekat{k}, responsekat{j},1); %#ok<PROP>
                        WpKat{k,j} = Wp; %CStat.Klouzaveokno(Wp,itimewindow(1),'max',1); %#ok<PROP>
                    end
                end
                obj.Wp.WpKat = WpKat;
            end
            
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
                line([obj.els(j)+0.5 obj.els(j)+0.5],[1 size(CC,1)],'color','black');
                line([1 size(CC,1)],[obj.els(j)+0.5 obj.els(j)+0.5],'color','black');
            end  
            colorbar;
        end
       
        function PlotCategory(obj,katnum,channel)
            %vykresli vsechny a prumernou odpoved na kategorii podnetu
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
        
        function PlotElectrode(obj,e,s,range,time)
            %vykresli data (2 sekundy ) z jedne elektrody e od vteriny zaznamu s
            %osa y je v rozmezi [-r +r]
            %zatim jen neepochovana data
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
            if  allels==1  %chci zobrazit vsechny elektrody
                if mod(e,2) 
                        elmin = 1; 
                        elmax = obj.els( floor(numel(obj.els)/2));
                        els = obj.els(1 : find(obj.els==elmax)); %#ok<PROP> 
                        els(2,1) = 1;   %#ok<PROP>  %do els davam hranice elektrod k postupnemu vykresleni - v prvni radce jsou konce kazde elektrody                       
                else
                        elmin = obj.els( floor(numel(obj.els)/2))+1; 
                        elmax =  max(obj.els);
                        els = obj.els( find(obj.els > elmin, 1 ) : size(obj.els,2)); %#ok<PROP>
                        els(2,1) = elmin; %#ok<PROP>
                end  
                els(2,2:end) = els(1,1:end-1)+1; %#ok<PROP> %doplnim dolni radku - zacatky kazde elektrody
            else
                if e==1, elmin = 1; else elmin = obj.els(e-1)+1; end %index prvni elektrody kterou vykreslit
                elmax = obj.els(e);            % index posledni elektrody kterou vykreslit
                els = [elmax; elmin]; %#ok<PROP>
            end
            
            % -------- ziskam data k vykresleni do promenne dd -----------------
            if obj.epochs <= 1 %pokud data jeste nejsou epochovana
                iD = [ (s-1)*obj.fs + 1,  (s-1)*obj.fs + obj.fs*time]; %indexy eeg, od kdy do kdy vykreslit
                dd = obj.d( iD(1) : iD(2),elmin: elmax)' ;   %data k plotovani - prehodim poradi, prvni jsou kanaly
                t = linspace(iD(1)/obj.fs, iD(2)/obj.fs, iD(2)-iD(1)+1); %casova osa            
            else %pokud data uz jsou epochovana 
                time_n = time*obj.fs; %kolik vzorku v case chci zobrazit
                dd = zeros(elmax-elmin+1,time_n);
                time_nsum = 0; %kolik vzorku v case uz mam v poli dd
                ss = s-1; %cislo kreslene epochy
                while time_nsum < time_n %pridavam hodnoty do pole dd z nekolika epoch, dokud nenaplnim pozadovanou delku
                    ss = ss+1;
                    if time_n - time_nsum >= obj.samples 
                        if ss <= obj.epochs
                            dd(:,1+time_nsum : time_nsum+obj.samples) = squeeze(obj.d(:, elmin: elmax,ss))';  %data k plotovani - prehodim poradi, prvni jsou kanaly
                        end
                        time_nsum = time_nsum + obj.samples;
                    else %pokud mi nezbyva cela epocha do konce pozadovaneho casoveho rozsahu
                        if ss <= obj.epochs
                            dd(:,1+time_nsum : time_n ) = squeeze(obj.d(1:time_n - time_nsum, elmin: elmax,ss))';
                        end
                        time_nsum = time_n ;
                    end
                    
                end
                t = linspace(obj.epochtime(1), time+obj.epochtime(1), time_n); %casova osa pres nekolik epoch
            end
            % -------- KRESLIM -----------------
            %kod viz navod zde https://uk.mathworks.com/matlabcentral/newsreader/view_thread/294163            
            mi = repmat(-range,[size(dd,1) 1]); % rozsahu osy y mi:ma
            ma = repmat(+range,[size(dd,1) 1]);                        
            shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
            shift = repmat(shift,1,size(dd,2));
            colors = [ 'b' 'k'];
            c = 0;
            for el = els              %#ok<PROP>
                rozsahel = (el(2):el(1))-els(2,1)+1;  %#ok<PROP>
                rozsahel1 = setdiff(rozsahel,obj.RjCh); %nerejectovane kanaly                
                plot(t, bsxfun(@minus,shift(end,:),shift( rozsahel1,:)) + dd( rozsahel1,:) ,colors(c+1) );  
                %lepsi bude je nezobrazovat
                %rozsahel0 = intersect(rozsahel,obj.RjCh); %rejectovane kanaly                
                %if numel(rozsahel0)>0 
                    %hrj = plot(t, bsxfun(@minus,shift(end,:),shift( rozsahel0,:)) + dd( rozsahel0,:) ,'color',[.8 .8 .8] );  
                    %uistack(hrj,'bottom');
                %end
                hold on;
                c = 1-c;
            end
            hold off;
            set(gca,'ytick',shift(:,1),'yticklabel',elmax:-1:elmin); %znacky a popisky osy y
            grid on;
            
            ylim([min(min(shift))-range max(max(shift))+range]); %rozsah osy y
            ylabel(['Electrode ' num2str(e) '/' num2str(numel(obj.els)) ]);
            xlabel(['Seconds of ' num2str( round(obj.samples*obj.epochs/obj.fs)) ]);
            text(t(1),-shift(2,1),[ 'resolution +/-' num2str(range) 'mV']);         
            xlim([t(1) t(end)]);
            
            % -------- ulozim  handle na obrazek a nastaveni grafu -----------------
            methodhandle = @obj.hybejPlot;
            set(obj.plotH,'KeyPressFcn',methodhandle); 
            
            obj.plotES = [e s range time allels ]; %ulozim hodnoty pro pohyb klavesami
            obj.epochLast = max([s obj.epochLast]); %oznaceni nejvyssi navstivene epochy
            
            for j = 1:size(shift,1)
                yshift = shift(end,1)-shift(j,1);
                text(t(end),yshift,[ ' ' obj.CH.H.channels(1,elmin+j-1).neurologyLabel ',' obj.CH.H.channels(1,elmin+j-1).ass_brainAtlas]);
                text(t(1)-size(dd,2)/obj.fs/10,yshift,[ ' ' obj.CH.H.channels(1,elmin+j-1).name]);
                if find(obj.RjCh==elmin-1+j) %oznacim vyrazene kanaly
                    text(t(1),yshift+50,' REJECTED');
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
                if allels==1, ty = -shift(4,1); else ty = -shift(2,1); end
                text(t(end)-((t(end)-t(1))/10),ty,[ 'excluded ' num2str(numel(obj.RjEpoch))]); 
            end
            
        end
        
        function PlotResponses(obj)
            %vykresli uspesnost odpovedi spolu s chybami a vyrazenymi epochami
            obj.PsyData.PlotResponses();
            figure(obj.PsyData.fhR); %kreslim dal do stejneho obrazku
            [~,rt,kategorie] = obj.PsyData.GetResponses();            
            plot(obj.RjEpoch,rt(obj.RjEpoch),'*','MarkerSize',10,'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]); %vykreslim vyrazene epochy
            plot(obj.RjEpoch,kategorie(obj.RjEpoch),'*r','MarkerSize',5); %vykreslim vyrazene epochy
        end
        
        function PlotResponseCh(obj,ch)
            %vykresli odpovedi pro jednotlivy kanal
            assert(obj.epochs > 1,'only for epoched data');
            if ~exist('ch','var')
                if isfield(obj.plotRCh,'ch')
                    ch = obj.plotRCh.ch;
                else
                    ch = 1;
                end
            end
            T = linspace(obj.epochtime(1),obj.epochtime(2),size(obj.d,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            if isfield(obj.plotRCh,'fh')
                figure(obj.plotRCh.fh); %pouziju uz vytvoreny graf
                clf(obj.plotRCh.fh); %graf vycistim
            else
                obj.plotRCh.fh = figure('Name','W plot channel');
            end
            
            iEp=obj.GetEpochsExclude();  
            M = mean(obj.d(:,ch,iEp),3);                
            E = std(obj.d(:,ch,iEp),[],3)/sqrt(size(obj.d,3)); %std err of mean
            errorbar(T,M,E,'.','color',[.8 .8 .8]); %nejdriv vykreslim errorbars aby byly vzadu
            hold on;
            plot(T,M,'LineWidth',2);  %prumerna odpoved              
            xlim(obj.epochtime);
            if isfield(obj.plotRCh,'ylim')
                ylim( obj.plotRCh.ylim);
                ymax = obj.plotRCh.ylim(2);
            else
                ymax = 0; ymin = 0;
                for katnum=obj.PsyData.Categories()
                    katdata = obj.CategoryData(katnum); %epochy jedne kategorie
                    ymax = max([ ymax max(mean(katdata(:,:,:),3))]);
                    ymin = min([ ymin min(mean(katdata(:,:,:),3))]);
                end           
                ylim( [ymin ymax]);
                obj.plotRCh.ylim = [ymin ymax];
            end
            if isfield(obj.Wp,'D2') %krivka p hodnot z W testu
                Tr = linspace(0,obj.epochtime(2),size(obj.Wp.D2,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
                plot(Tr,obj.Wp.D2(:,ch),'b:');  %carkovana modra cara oznacuje signifikanci prumeru
                iWp = obj.Wp.D2(:,ch) <= 0.05;
                plot(Tr(iWp),obj.Wp.D2(iWp,ch),'m.'); %fialove jsou p < 0.05
                iWp = obj.Wp.D2(:,ch) <= 0.01;
                plot(Tr(iWp),obj.Wp.D2(iWp,ch),'r.'); %cervene jsou p < 0.01
            end
            colorskat = 'kgr';
            for katnum = obj.PsyData.Categories()
                katdata = obj.CategoryData(katnum); %epochy jedne kategorie
                M = mean(katdata(:,ch,:),3);
                %E = std(katdata(:,ch,:),[],3)/sqrt(size(katdata,3)); %std err of mean
                plot(T,M,'LineWidth',1,'Color',colorskat(katnum+1));  %prumerna odpoved  
                if isfield(obj.Wp,'WpKat')
                    for l = katnum+1:numel(obj.PsyData.Categories())-1 %katnum jde od nuly
                        if katnum==0, color=colorskat(l+1); else color = colorskat(1); end %green a red jsou proti kategorii 0, cerna je kat 1 vs kat 2
                        plot(Tr,obj.Wp.WpKat{katnum+1,l+1}(:,ch),[color ':']); %carkovana cara oznacuje signifikanci kategorie vuci jine kategorii
                        iWp = obj.Wp.WpKat{katnum+1,l+1}(:,ch)  <= 0.1; 
                        plot(Tr(iWp),ones(1,sum(iWp))*-0.2, [ color '.']); %
                        iWp = obj.Wp.WpKat{katnum+1,l+1}(:,ch)  <= 0.05;               
                        plot(Tr(iWp),ones(1,sum(iWp))*-0.2, [ color '*']); %
                    end
                end
            end
            title(['channel ' num2str(ch)]);
            text(-0.1,ymax*.95,[ obj.CH.H.channels(1,ch).name ' : ' obj.CH.H.channels(1,ch).neurologyLabel ',' obj.CH.H.channels(1,ch).ass_brainAtlas]);
            
            methodhandle = @obj.hybejPlotCh;
            set(obj.plotRCh.fh,'KeyPressFcn',methodhandle); 
            obj.plotRCh.ch = ch;
            
        end        
            
        function PlotResponseP(obj)
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
        
             
        
        
        %% SAVE AND LOAD FILE
        function Save(obj,filename)   
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
            save(filename,'d','tabs','tabs_orig','fs','header','sce','PsyDataP','epochtime','CH_H','els','plotES','RjCh','RjEpoch','epochTags','epochLast','reference','epochData','Wp','-v7.3');            
        end
        function obj = Load(obj,filename)
            % nacte veskere promenne tridy ze souboru
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
                if ismember('epochData', {vars.name}),           load(filename,'epochData');    obj.epochData = epochData;   end   %#ok<CPROP,PROP>             
                load(filename,'epochtime');                
                obj.epochtime = epochtime;      %#ok<CPROP,PROP>               
            else
                obj.epochtime = [];
                obj.epochData = [];
            end
            if ismember('CH_H', {vars.name})
                load(filename,'CH_H');      obj.CH = CHHeader(CH_H);
            else
                load(filename,'CH');
                obj.CH = CH; %#ok<CPROP,PROP> %  %drive ulozeny objekt, nez jsem zavedl ukladani struct
            end  
            if ismember('Wp', {vars.name})
                load(filename,'Wp');      obj.Wp = Wp;
            else
                obj.Wp = struct;
            end
            obj.els = els;                  %#ok<CPROP,PROP>
            obj.plotES = plotES;            %#ok<CPROP,PROP>
            %obj.plotH = plotH;              %#ok<CPROP,PROP>
            obj.RjCh = RjCh;                %#ok<CPROP,PROP>     
            obj.RjEpoch = RjEpoch;          %#ok<CPROP,PROP>
            if exist('epochTags','var'),  obj.epochTags = epochTags;   else obj.epochTags = []; end         %#ok<CPROP,PROP>     
            if exist('epochLast','var'),  obj.epochLast = epochLast;   else obj.epochLast = []; end         %#ok<CPROP,PROP>
            if exist('reference','var'),  obj.reference = reference;   else obj.reference = 'original'; end         %#ok<CPROP,PROP> %14.6.2016
           
            obj.filename = filename;          %
        end
    end
    %% privatni metody
    methods  (Access = private)
        function hybejPlot(obj,~,eventDat)           
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
                   if(obj.plotES(1))<numel(obj.els) %pokud je cislo elektrody ne maximalni
                        obj.PlotElectrode(obj.plotES(1)+1,obj.plotES(2),obj.plotES(3),obj.plotES(4));
                   end                   
               case 'uparrow'
                   if(obj.plotES(1))>1 %pokud je cislo elektrody vetsi nez 1
                        obj.PlotElectrode(obj.plotES(1)-1,obj.plotES(2),obj.plotES(3),obj.plotES(4));
                   end
               case 'add'     %signal mensi - vetsi rozliseni 
                   if obj.plotES(3)>=obj.yrange(3), pricist = obj.yrange(4);
                   else pricist = obj.yrange(2);                   
                   end
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3)+pricist,obj.plotES(4));
                   
               case 'subtract' %signal vetsi - mensi rozliseni   
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
        function hybejPlotCh(obj,~,eventDat)  
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
               case 'space' %zobrazi i prumerne krivky
                   obj.PlotResponseFreq(obj.plotRCh.ch);
                   figure(obj.plotRCh.fh); %dam puvodni obrazek dopredu
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
        
    end
    
end

