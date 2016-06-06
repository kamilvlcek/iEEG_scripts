classdef CiEEGData < handle
    %CEEGDATA Trida na praci s datama ve formatu ISARG od Petra Jezdika
    %   Kamil Vlcek, 2016-04
    
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
        plotES; % current electrode, second of plot/epoch, range of y values, time range
        plotH;  % handle to plot
        RjCh; %seznam cisel rejectovanych kanalu
        RjEpoch; %seznam vyrazenych epoch
    end
    
    methods (Access = public)
        function obj = CiEEGData(d,tabs,fs,mults,header)
            %konstruktor, parametry d,tabs,fs[,mults,header]
            if ischar(d) && ~exist('tabs','var') %pokud je prvni parametr retezec, tak ho beru jako nazev souboru, ktery nactu
                assert(exist(d, 'file') == 2,['Soubor' d 'neexistuje']);
                obj.Load(d);
            else
                obj.d = d;
                [obj.samples,obj.channels, obj.epochs] = obj.DSize();
                obj.tabs = tabs;
                obj.tabs_orig = tabs;
                obj.fs = fs;
                if exist('mults','var') && ~isempty(mults)
                    obj.mults = mults;
                else
                    obj.mults = ones(1,size(d,2)); %defaultove jednicky pro kazdy kanal
                end
                if exist('header','var')
                    obj.header = header;
                else
                    obj.header = [];
                end
                obj.plotES = [1 1 150 5 0]; %nastavim defaultni hodnoty grafy
            end
        end
        
        function [samples, channels, epochs] = DSize(obj)
            % vraci velikosti pole d - samples, channels, epochs
            samples = size(obj.d,1);
            channels = size(obj.d,2);
            epochs = size(obj.d,3);
        end
        
        function dd = DData(obj,ch,epoch)
            % vraci jeden kanal a jednu epochu ze zaznamu. Implementovano kvuli nasobeni mults
            if ~exist('epoch','var')
                epoch = 1;            
            end
            dd = double(obj.d(:,ch,epoch)) .* obj.mults(ch);                           
        end
        
        function obj = GetHHeader(obj,H)
            %nacte header z promenne H - 25.5.2016
            obj.CH = CHHeader(H);
            [~, obj.els] = obj.CH.ChannelGroups();            
        end
        
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
        end
        
        function obj = RejectChannels(obj,RjCh)
            %ulozi cisla vyrazenych kanalu - kvuli pocitani bipolarni reference 
            obj.RjCh = RjCh;
        end
        
        function obj = RejectEpochs(obj,RjEpoch)
            %ulozi cisla vyrazenych epoch - kvuli prevodu mezi touto tridou a CHilbert
            obj.RjEpoch = RjEpoch;
        end
        
        function [ranges]=PlotElectrode(obj,e,s,range,time)
            %vykresli data (2 sekundy ) z jedne elektrody e od vteriny zaznamu s
            %osa y je v rozmezi [-r +r]
            %zatim jen neepochovana data
            if ~exist('e','var'), e= obj.plotES(1); end %cislo elektrody
            if ~exist('s','var'), s= obj.plotES(2); end %cislo epochy nebo vteriny zaznamu
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
            
            if  allels==1  %chci zobrazit vsechny elektrody
                if mod(e,2) 
                        elmin = 1; elmax = obj.els( floor(numel(obj.els)/2));
                else
                        elmin = obj.els( floor(numel(obj.els)/2))+1; elmax =  max(obj.els);
                end
                
            else
                if e==1, elmin = 1; else elmin = obj.els(e-1)+1; end %index prvni elektrody kterou vykreslit
                elmax = obj.els(e);            % index posledni elektrody kterou vykreslit
            end
           
            if obj.epochs <= 1 %pokud data jeste nejsou epochovana
                iD = [ (s-1)*obj.fs + 1,  (s-1)*obj.fs + obj.fs*time]; %indexy eeg, od kdy do kdy vykreslit
                dd = obj.d( iD(1) : iD(2),elmin: elmax)';   %data k plotovani - prehodim poradi, prvni jsou kanaly
                t = linspace(iD(1)/obj.fs, iD(2)/obj.fs, iD(2)-iD(1)+1); %casova osa            
            else %pokud data uz jsou epochovana   
                dd = squeeze(obj.d(:, elmin: elmax,s))';  %data k plotovani - prehodim poradi, prvni jsou kanaly
                t = linspace(obj.epochtime(1), obj.epochtime(2), size(obj.d,1)); %casova osa    
            end
            %kod viz navod zde https://uk.mathworks.com/matlabcentral/newsreader/view_thread/294163
            if exist('range','var') && ~isempty(e)
                mi = repmat(-range,[size(dd,1) 1]);
                ma = repmat(+range,[size(dd,1) 1]);
            else
                mi = min(dd,[],2);          
                ma = max(dd,[],2);
                e = [];
            end
            
            shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
            shift = repmat(shift,1,size(dd,2));            
            plot(t,dd+shift);
            set(gca,'ytick',shift(:,1),'yticklabel',elmin:elmax);
            grid on;
            ylim([min(min(shift))-range max(max(shift))+range]); 
            ylabel(['Electrode ' num2str(e) '/' num2str(numel(obj.els)) ]);
            xlabel(['Seconds of ' num2str( round(obj.samples/obj.fs)) ]);
            text(t(1),-shift(2,1),[ 'resolution +/-' num2str(range) 'mV']); 
            
            methodhandle = @obj.hybejPlot;
            set(obj.plotH,'KeyPressFcn',methodhandle); 
            ranges = [mi ma];
            obj.plotES = [e s range time allels]; %ulozim hodnoty pro pohyb klavesami
            for j = 1:size(shift,1)
                text(t(end),shift(j,1),[ ' ' obj.CH.H.channels(1,elmin+j-1).neurologyLabel]);
                text(t(1)-size(dd,2)/obj.fs/10,shift(j,1),[ ' ' obj.CH.H.channels(1,elmin+j-1).name]);
                if find(obj.RjCh==elmin-1+j) %oznacim vyrazene kanaly
                    text(t(1),shift(j,1)+50,' REJECTED');
                end
            end  
            if obj.epochs > 1
                titul = ['Epoch ' num2str(s) '/' num2str(obj.epochs)];
                if find(obj.RjEpoch==s)
                    titul = [titul ' - EXCLUDED'];                    
                    line([t(1) t(end)],[shift(1,1) shift(end,1)],'Color','r');
                end
                title(titul);
            end
            
        end
        
        
        function ExtractEpochs(obj, PsyData,epochtime)
            % epochuje data v poli d, pridava do objektu:
            % cell array epochData, double(2) epochtime v sekundach, tridu PsyData 
            % upravuje obj.mults, samples channels epochs
            if obj.epochs > 1
                disp('already epoched data');
                return;
            end
            obj.PsyData = PsyData; %objekt CPsyData
            obj.epochtime = epochtime; %v sekundach cas pred a po udalosti  , prvni cislo je zaporne druhe kladne
            iepochtime = round(epochtime.*obj.fs); %v poctu vzorku cas pred a po udalosti
            ts_podnety = PsyData.TimeStimuli(); %timestampy vsech podnetu
            de = zeros(iepochtime(2)-iepochtime(1), size(obj.d,2), size(ts_podnety,1)); %nova epochovana data time x channel x epoch            
            tabs = zeros(iepochtime(2)-iepochtime(1),size(ts_podnety,1)); %#ok<PROP> %udelam epochovane tabs
            obj.epochData = cell(size(ts_podnety,1),3); % sloupce kategorie, cislo kategorie, timestamp
            for epoch = 1:size(ts_podnety,1) %pro vsechny eventy
                izacatek = find(obj.tabs<=ts_podnety(epoch), 1, 'last' ); %najdu index podnetu podle jeho timestampu
                    %kvuli downsamplovani Hilberta, kdy se mi muze ztratit presny cas zacatku
                    %epochy, beru posledni nizsi tabs nez je cas zacatku epochy
                [Kstring Knum] = PsyData.Category(epoch);    %jmeno a cislo kategorie
                obj.epochData(epoch,:)= {Kstring Knum obj.tabs(izacatek)}; %zacatek epochy beru z tabs aby sedel na tabs pri downsamplovani
                for ch = 1:obj.channels %pro vsechny kanaly                    
                    baseline = mean(obj.d(izacatek+iepochtime(1):izacatek-1));
                    de(:,ch,epoch) = double(obj.d( izacatek+iepochtime(1) : izacatek+iepochtime(2)-1,ch)).* obj.mults(ch) - baseline; 
                    tabs(:,epoch) = obj.tabs(izacatek+iepochtime(1) : izacatek+iepochtime(2)-1); %#ok<PROP>
                end
            end
            obj.d = de; %puvodni neepochovana budou epochovana
            obj.mults = ones(1,size(obj.d,2)); %nove pole uz je double defaultove jednicky pro kazdy kanal
            obj.tabs = tabs; %#ok<PROP>
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
        end
        
        function [d]= CategoryData(obj, katnum)
            %vraci epochy ve kterych podnet byl kategorie/podminky katnum
            assert(obj.epochs > 1,'data not yet epoched'); %vyhodi chybu pokud data nejsou epochovana
            iEpochy = cell2mat(obj.epochData(:,2))==katnum ; %seznam epoch v ramci kategorie ve sloupci
            d = obj.d(:,:,iEpochy);
        end
        
        function PlotCategory(obj,katnum,channel)
            %vykresli vsechny a prumernou odpoved na kategorii podnetu
            d1=obj.CategoryData(katnum); %epochy jedne kategorie
            d1m = mean(d1,3); %prumerne EEG z jedne kategorie
            T = (0 : 1/obj.fs : (size(obj.d,1)-1)/obj.fs) + obj.epochtime(1); %cas zacatku a konce epochy
            E = 1:obj.epochs; %vystupni parametr
            h1 = figure('Name','Mean Epoch'); %#ok<NASGU> %prumerna odpoved na kategorii
            plot(T,d1m(:,channel));
            xlabel('Time [s]'); 
            title(obj.PsyData.CategoryName(katnum));
            h2 = figure('Name','All Epochs');  %#ok<NASGU> % vsechny epochy v barevnem obrazku
            imagesc(T,E,squeeze(d1(:,channel,:))');
            colorbar;
            xlabel('Time [s]');
            ylabel('Epochs');
            title(obj.PsyData.CategoryName(katnum));
            
        end
        
        function Save(obj,filename)   
            %ulozi veskere promenne tridy do souboru
            d = obj.d;                      %#ok<PROP,NASGU>            
            tabs = obj.tabs;                %#ok<PROP,NASGU>
            tabs_orig = obj.tabs_orig;      %#ok<PROP,NASGU>
            fs = obj.fs;                    %#ok<PROP,NASGU>
            mults = obj.mults;              %#ok<PROP,NASGU>
            header = obj.header;            %#ok<PROP,NASGU>
            sce = [obj.samples obj.channels obj.epochs]; %#ok<NASGU>
            PsyData = obj.PsyData;          %#ok<PROP,NASGU>
            epochtime = obj.epochtime;      %#ok<PROP,NASGU>
            CH=obj.CH;                      %#ok<PROP,NASGU>
            els = obj.els;                  %#ok<PROP,NASGU>
            plotES = obj.plotES;            %#ok<PROP,NASGU>
            plotH = obj.plotH;              %#ok<PROP,NASGU>
            RjCh = obj.RjCh;                %#ok<PROP,NASGU>
            RjEpoch = obj.RjEpoch;          %#ok<PROP,NASGU>
            save(filename,'d','tabs','tabs_orig','fs','mults','header','sce','PsyData','epochtime','CH','els','plotES','plotH','RjCh','RjEpoch','-v7.3');            
        end
        function obj = Load(obj,filename)
            % nacte veskere promenne tridy ze souboru
            load(filename,'d','tabs','tabs_orig','fs','mults','header','sce','PsyData','epochtime','CH','els','plotES','plotH','RjCh','RjEpoch');            
            obj.d = d;                      %#ok<CPROP,PROP>
            obj.tabs = tabs;                %#ok<CPROP,PROP>
            obj.tabs_orig = tabs_orig;      %#ok<CPROP,PROP>
            obj.fs = fs;                    %#ok<CPROP,PROP>           
            obj.mults = mults;              %#ok<CPROP,PROP>
            obj.header = header;            %#ok<CPROP,PROP>
            obj.samples = sce(1); obj.channels=sce(2); obj.epochs = sce(3); 
            obj.PsyData = PsyData;          %#ok<CPROP,PROP>
            obj.epochtime = epochtime;      %#ok<CPROP,PROP>
            obj.CH = CH;                    %#ok<CPROP,PROP>
            obj.els = els;                  %#ok<CPROP,PROP>
            obj.plotES = plotES;            %#ok<CPROP,PROP>
            obj.plotH = plotH;              %#ok<CPROP,PROP>
            obj.RjCh = RjCh;                %#ok<CPROP,PROP>     
            obj.RjEpoch = RjEpoch;          %#ok<CPROP,PROP>
        end
    end
    methods  (Access = private)
        function hybejPlot(obj,~,eventDat)           
           switch eventDat.Key
               case 'rightarrow' 
                   if obj.epochs == 1
                       rightval = obj.plotES(2)+obj.plotES(4);
                       maxval = round(obj.samples/obj.fs);
                   else
                       rightval = obj.plotES(2);
                       maxval = obj.epochs; %pocet epoch
                   end
                   if( rightval < maxval)   %pokud je cislo vteriny vpravo mensi nez celkova delka                    
                        obj.PlotElectrode(obj.plotES(1),obj.plotES(2)+1,obj.plotES(3),obj.plotES(4));
                   end
               case 'leftarrow'
                   if(obj.plotES(2))>1 %pokud je cislo vteriny vetsi nez 1
                        obj.PlotElectrode(obj.plotES(1),obj.plotES(2)-1,obj.plotES(3),obj.plotES(4));
                   end
               case 'pagedown'
                   if obj.epochs == 1
                       rightval = obj.plotES(2)+obj.plotES(4);
                       maxval = round(obj.samples/obj.fs)-10;
                   else
                       rightval = obj.plotES(2);
                       maxval = obj.epochs-10; %pocet epoch
                   end
                   if( rightval < maxval)   %pokud je cislo vteriny vpravo mensi nez celkova delka
                        obj.PlotElectrode(obj.plotES(1),obj.plotES(2)+10,obj.plotES(3),obj.plotES(4));
                   end
               case 'pageup'
                   if(obj.plotES(2))>10 %pokud je cislo vteriny vetsi nez 1
                        obj.PlotElectrode(obj.plotES(1),obj.plotES(2)-10,obj.plotES(3),obj.plotES(4));
                   end
               case 'home'     % na zacatek zaznamu              
                        obj.PlotElectrode(obj.plotES(1),1,obj.plotES(3),obj.plotES(4));                   
               case 'uparrow'
                   if(obj.plotES(1))<numel(obj.els) %pokud je cislo elektrody ne maximalni
                        obj.PlotElectrode(obj.plotES(1)+1,obj.plotES(2),obj.plotES(3),obj.plotES(4));
                   end                   
               case 'downarrow'
                   if(obj.plotES(1))>1 %pokud je cislo elektrody vetsi nez 1
                        obj.PlotElectrode(obj.plotES(1)-1,obj.plotES(2),obj.plotES(3),obj.plotES(4));
                   end
               case 'add'     %signal mensi - vetsi rozliseni 
                   if obj.plotES(3)>=50, pricist = 50;
                   else pricist = 10;                   
                   end
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3)+pricist,obj.plotES(4));
                   
               case 'subtract' %signal vetsi - mensi rozliseni   
                   if obj.plotES(3)>50, odecist = 50;
                   elseif obj.plotES(3)>10 odecist = 10;
                   else odecist = 0;
                   end
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3)-odecist,obj.plotES(4));
               
               case 'space' %epoch exclusion
                   s = obj.plotES(2);
                   if find(obj.RjEpoch== s)                      
                        obj.RjEpoch = obj.RjEpoch(obj.RjEpoch~=s); %vymazu hodnoty s
                   else
                        obj.RjEpoch = [obj.RjEpoch  obj.plotES(2)]; %pridam hodnotu s
                   end   
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3),obj.plotES(4));
                   
               case 'return'  %prehazuje mezi zobrazeni jednotlivych elektrod a cele poloviny elektrod                   
                   obj.plotES(5) = 1-obj.plotES(5);
                   obj.PlotElectrode(obj.plotES(1),obj.plotES(2),obj.plotES(3),obj.plotES(4));
               otherwise
                   disp(['You just pressed: ' eventDat.Key]);                      
           end
        end
        
    end
    
end

