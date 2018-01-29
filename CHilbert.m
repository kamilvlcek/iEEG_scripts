classdef CHilbert < CiEEGData
    %HILBERT.CLASS sbirka funkci na analyzu pomoci hilbert trasform
    %   Kamil Vlcek, FGU AVCR, since 2016 04
    properties (Constant = true)
        decimatefactor = 8; %o kolik decimuju hiblertovu obalku oproti puvodi sampling rate; 2 je dostatecne konzervativni, 8 hodne setri pamet
    end
    properties (Access = public)
        HFreq; %hilberova obalka pro kazde frekvenci pasmo - time x channel x freq (x kategorie)
        HFreqEpochs; %Hf bez priemerovania cez epochy - time x channel x frequency x epoch
        Hf; %frekvencni pasma pro ktere jsou pocitany obalky - okraje pasem, pocet je tedy vetsi o 1 nez pocet spocitanych pasem
        Hfmean; %stredni hodnoty pasem  - pocet = pocet spocitanych pasem
        hfilename; %jmeno souboru CHilbert  
        plotF = struct; %udaje o stavu plotu PlotResponseFreq
        plotEpochs = struct; %udaje o stavu plotu PlotMovingEpochs - Nada
    end
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS 
        function obj = CHilbert(d,tabs,fs,mults,header)            
            if ~exist('header','var'), header = []; end %nejakou hodnotu dat musim
            if ~exist('mults','var'),  mults = []; end %nejakou hodnotu dat musim
            if ~exist('d','var') %konstruktor uplne bez parametru - kvuli CHilbertMulti
                d = []; tabs = []; fs = [];
            elseif ischar(d) && ~exist('tabs','var') %pokud je prvni parametr retezec, tak ho beru jako nazev souboru, ktery nactu
                tabs=[]; fs = [];
                % volani Load z CiEEGData mi zavola Load z CHilbert, takze d=filename predelavat nemusim
            end            
            obj@CiEEGData(d,tabs,fs,mults,header); %volani konstruktoru nemuze byt v if bloku 
            try                
                if ~isempty(obj.Hf)
                    disp(['Frequency bands: ' num2str(numel(obj.Hf)) ': ' num2str(obj.Hf(1)) ':' num2str(obj.Hf(2)-obj.Hf(1)) ':' num2str(obj.Hf(end)) ' Hz' ]);
                else
                    disp('no Frequency bands');
                end
            catch exception %#ok<NASGU>
                disp('no Frequency bands');
            end
            
        end        
        
        function obj = PasmoFrekvence(obj,freq,channels,prekryv,decimatefactor)
            %EEG2HILBERT prevede vsechny kanaly na prumer hilbertovych obalek
            %   podle puvodni funkce EEG2Hilbert
            %   pouziva data d z parentu a take fs
            %   freq    seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi
            if ~exist('channels','var') || isempty(channels), channels = 1:obj.channels; end
            if ~exist('prekryv','var'), prekryv = 0; end %kolik se maji prekryvat sousedni frekvencni pasma, napriklad 0.5
            if ~exist('decimatefactor','var') || isempty(decimatefactor), decimatefactor = obj.decimatefactor; end; %volitelny parametr decimatefactor 
            samples = ceil(obj.samples/decimatefactor); 
            disp(['vytvarim pole ' num2str(samples) 'x' num2str(obj.channels) 'x' num2str(numel(freq)-1) ... 
                '=' num2str(samples*obj.channels*(numel(freq)-1)*8/1024/1024) ' MBytes']); %zpravu abych vedel, v jakych velikostech se pohybuju
            obj.HFreq = zeros(samples,obj.channels,numel(freq)-1); %inicializace pole   
            tic; %zadnu meric cas
            fprintf('kanal ze %i: ', max(channels) );
            for ch = channels %jednotlive elektrody
                %fprintf('channel %i: Hz ',ch);
                           
                fprintf('%i,',ch);
                for fno = 1:numel(freq)-1 %seznam frekvenci
                    loF = freq(fno) -prekryv*(freq(fno+1)-freq(fno)); 
                    hiF = freq(fno+1)-0.1 +prekryv*(freq(fno+1)-freq(fno));  %napr 50 - 59.9
                    hh = obj.hilbertJirka(obj.d(:,ch),loF,hiF,obj.fs); %cista hilbertova obalka, tohle i skript hodne zrychli
                    hh = decimate(hh,decimatefactor); % mensi sampling rate                    
                    obj.HFreq(:,ch,fno) = (hh./mean(hh)); %podil prumeru = prumerna hodnota
                    %fprintf('%i Hz, ',loF);
                end
                %fprintf('\n'); %tisk znova na stejnou radku
            end
            obj.d = squeeze(mean(obj.HFreq,3)); %11.5.2016 - prepisu puvodni data prumerem
            obj.fs = obj.fs/decimatefactor;
            obj.tabs = downsample(obj.tabs,decimatefactor);
            obj.tabs_orig = downsample(obj.tabs_orig,decimatefactor); %potrebuju zdecimovat i druhy tabs. Orig znamena jen ze nepodleha epochovani
            obj.Hf = freq;
            obj.Hfmean = (freq(1:end-1) + freq(2:end)) ./ 2;            
            obj.mults = ones(1,size(obj.d,2)); %nove pole uz je double defaultove jednicky pro kazdy kanal
            obj.yrange = [1 1 5 5]; %zmenim rozliseni osy y v grafu
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
            fprintf('\n'); %ukoncim radku
            toc; %ukoncim mereni casu a vypisu
            obj.DatumCas.HilbertComputed = datestr(now);
            disp(['vytvoreno ' num2str(numel(obj.Hfmean)) ' frekvencnich pasem']); 
        end
        
        function obj = ExtractEpochs(obj, PsyData,epochtime,baseline)
            % rozdeli hilbertovu obalku podle epoch
            % i u ni odecte baseline pred podnetem
            if ~exist('baseline','var') || isempty(baseline), baseline = [epochtime(1) 0]; end %defaultni baseline je do 0 sec
            ExtractEpochs@CiEEGData(obj,PsyData, epochtime,baseline); %to mi zepochuje prumernou obalku za frekvencni pasma v poli d
            if(numel(obj.HFreq)>0)
                %ted epochace vsech frekvencnich pasem zvlast, hlavne kvuli obrazkum
                %prumer za kazdou kategorii, statistiku z toho delat nechci
                 iepochtime = round(epochtime(1:2).*obj.fs); %v poctu vzorku cas pred a po udalosti, prvni cislo je zaporne druhe kladne             
                 ibaseline =  round(baseline.*obj.fs); %v poctu vzorku cas pred a po udalosti
                 kategorie = cell2mat(obj.PsyData.P.strings.podminka(:,2)); %cisla karegorii ve sloupcich
                 Hfreq2 = zeros(iepochtime(2)-iepochtime(1), size(obj.d,2), numel(obj.Hfmean),size(kategorie,1)); %nova epochovana data time x channel x freq x kategorie=podminka
                 obj.HFreqEpochs = zeros(iepochtime(2)-iepochtime(1),size(obj.HFreq,2),size(obj.HFreq,3),obj.epochs); % time x channel x frequency x epoch
                 %cyklus po kategoriich ne po epochach
                 for katnum = kategorie' %potrebuji to v radcich
                     Epochy = find(cell2mat(obj.epochData(:,2))==katnum); %seznam epoch v ramci kategorie ve sloupci 
                     for epoch = Epochy' %potrebuji to v radcich
                         izacatek = find(obj.tabs_orig==obj.epochData{epoch,3}); %najdu index podnetu, podle jeho timestampu. v tretim sloupci epochData jsou timestampy
                         for ch=1:obj.channels
                            baseline_mean = mean(obj.HFreq(izacatek + ibaseline(1) : izacatek+ibaseline(2)-1,ch,:),1); %baseline pro vsechny frekvencni pasma dohromady
                            epoch_data = bsxfun(@minus,obj.HFreq(izacatek + iepochtime(1) : izacatek+iepochtime(2)-1,ch,:) , baseline_mean); %odecteni baseline pro aktualni epochu a kanal
                            obj.HFreqEpochs(:,ch,:,epoch) = epoch_data;
                            Hfreq2(:,ch,:,katnum+1) = Hfreq2(:,ch,:,katnum+1) + epoch_data; %soucet power pro kategorii, pres prislusne epochy
                            %tady se mi to mozna odecetlo blbe? KOntrola
                         end
                     end
                     Hfreq2(:,:,:,katnum+1) = Hfreq2(:,:,:,katnum+1)./numel(Epochy); %prumer pred epochy - soucet podelim prumerem
                 end             
                 obj.HFreq = Hfreq2;
            end
        end
        
        function obj = Decimate(obj,podil,rtrim)
            %zmensi frekvencni data na nizsi vzorkovaci frekvenci
            if ~exist('rtrim','var') || isempty(rtrim), rtrim = []; end 
            Decimate@CiEEGData(obj,podil,rtrim);
            if obj.decimatefactor == 1 %zatim pouzivam pouze, pokud nejsou data uz decimovana vuci CiEEGdata
                fprintf('channels to decimate HFreq (z %i):',numel(obj.channels));
                HFreq = zeros(ceil(size(obj.HFreq,1)/podil) , size(obj.HFreq,2), size(obj.HFreq,3),size(obj.HFreq,4));  %#ok<PROPLC>
                HFreqEpochs = zeros(ceil(size(obj.HFreq,1)/podil) , size(obj.HFreq,2), size(obj.HFreq,3),obj.epochs); %#ok<PROPLC>
                for ch = 1:obj.channels                    
                    fprintf('%i, ',ch);
                    for f = 1:size(obj.HFreq,3) %pocet frekvencnich pasem
                       for kat = 1:size(obj.HFreq,4) %pocet kategorii podnetu     
                            HFreq(:,ch,f,kat) = decimate(obj.HFreq(:,ch,f,kat),podil); %#ok<PROPLC>                            
                       end
                       for ep = 1:obj.epochs %frekvencni data se vsemi epochami
                            HFreqEpochs(:,ch,f,ep) = decimate(obj.HFreqEpochs(:,ch,f,ep),podil); %#ok<PROPLC>    
                       end
                    end                    
                end
                obj.HFreq = HFreq; %#ok<PROPLC>
                obj.HFreqEpochs = HFreqEpochs; %#ok<PROPLC>
                fprintf('... done\n');
                if exist('rtrim','var') && ~isempty(rtrim)
                    obj.HFreq = obj.HFreq(1:rtrim,:,:,:);                    
                end
            end
        end
        
        function obj = PlotResponseFreq(obj,ch,kategories)
            %uchovani stavu grafu, abych ho mohl obnovit a ne kreslit novy
            if ~exist('ch','var')
                if isfield(obj.plotF,'ch'), ch = obj.plotF.ch;
                else ch = 1; obj.plotF.ch = ch; end
            else
                obj.plotF.ch = ch;
            end
            
            if ~exist('kategories','var')
                if isfield(obj.Wp, 'kats')
                    kategories = obj.Wp.kats; 
                elseif isfield(obj.plotF,'kategories')
                    kategories = obj.plotF.kategories;
                else
                    kategories = obj.PsyData.Categories(); 
                    obj.plotF.kategories = kategories;
                end
            else
                obj.plotF.kategories = kategories;
            end
            
            if isfield(obj.plotF,'fh') && ishandle(obj.plotF.fh)
                figure(obj.plotF.fh); %pouziju uz vytvoreny graf
                %clf(obj.plotF.fh); %graf vycistim
            else
                obj.plotF.fh = figure('Name','ResponseFreq','Position', [20, 500, 1200, 300]);
                colormap jet; %aby to bylo jasne u vsech verzi matlabu - i 2016
            end   
            
            maxy = 0;
            miny = 0;
            for k = 1:numel(kategories)
                subplot(1,numel(kategories),k);
                T = obj.epochtime(1):0.1:obj.epochtime(2);
                F =  obj.Hfmean;
                D = squeeze(obj.HFreq(:,ch,:,kategories(k)+1));
                imagesc(T,F, D');
                maxy = max([maxy max(max( D ))]);
                miny = min([miny min(min( D ))]);                
                axis xy;               
                xlabel('time [s]');                
            end
            if isfield(obj.plotF,'ylim') && numel(obj.plotF.ylim)>=2 %nactu nebo ulozim hodnoty y
                miny = obj.plotF.ylim(1); maxy = obj.plotF.ylim(2);
            else
                obj.plotF.ylim = [miny maxy];
            end
            for k = 1:numel(kategories)
                subplot(1,numel(kategories),k);
                caxis([miny,maxy]);
                title( obj.PsyData.CategoryName(kategories(k)));
                if k == 1, ylabel(['channel ' num2str(ch) ' - freq [Hz]']); end
                if k == numel(kategories), colorbar('Position',[0.92 0.1 0.02 0.82]); end
            end
            
            methodhandle = @obj.hybejPlotF;
            set(obj.plotF.fh,'KeyPressFcn',methodhandle);             
        end 
        
        %% PLOT FUNCTIONS
        function PlotAllEpochs(obj, icondition, channel,zlimits)
            % PlotAllEpochs(obj, icondition, channel,ylimits) - cislo podminky, kanal, rozsah z osy
            % plots all available time x frequency epoch maps for given channel and condition (aedist: 1 = ego, 2 = allo, 0 = red)
            % Nada since 2018/01
            assert(~isempty(obj.HFreqEpochs),'soubor s frekvencnimi daty pro epochy neexistuje');
            if(~exist('zlimits','var')), zlimits =[-1 1]; end
            condition = obj.PsyData.CategoryName(icondition);  %zjisti jmeno kategorie z jejiho cisla
            figure('Name',[condition, ' - channel ', num2str(channel)], 'NumberTitle', 'off');
            time = linspace(obj.epochtime(1), obj.epochtime(2), size(obj.HFreqEpochs,1));
            
            epochs = find(obj.PsyData.P.data(:, 7) == icondition);
            chyby = obj.PsyData.GetErrorTrials();
            correct = 0;
            correct_epochs = zeros(90,90);
            subplot_x = ceil(sqrt(length(epochs)/15)*3); % dynamicke rozmery subplotu pre rozny pocet epoch (pre screen cca 3:5)
            subplot_y = ceil(sqrt(length(epochs)/15)*5); % vychadzam z 3x * 5x = pocet_epoch
            for i = 1:length(epochs)
                
               subplot(subplot_x, subplot_y, i);
               imagesc(squeeze(obj.HFreqEpochs(:,channel,:,epochs(i)))', 'XData', time, 'YData', obj.Hf);
               colormap parula; %aby to bylo jasne u vsech verzi matlabu - i 2010 - vypada v takovem mnozstvi lip nez jet
               caxis(zlimits);
               set(gca,'YDir','normal');
               title([num2str(obj.PsyData.P.data(epochs(i),5)), '-Epoch (', num2str(epochs(i)), ')'],'FontSize',8);
               
               % plot rejected line
               if obj.PlotRejected(channel, time, epochs(i), sum(chyby(epochs(i),:))) %jestli jde o vyrazenou epochu
                  correct = correct +1; % ukladam spravne epochy pre vykreslenie mapy mean frequency powers pre vsetky epochy 
                  correct_epochs(correct,:) = mean(squeeze(obj.HFreqEpochs(:,channel,:,epochs(i)))',1);
                  title(sprintf('%d - epoch %d (%d) ', obj.PsyData.P.data(epochs(i),5), correct, epochs(i)),'FontSize',8) ;
               end
               
               % plot response time 
               response_time = obj.PsyData.P.data(epochs(i), 4);
               hold on; 
               plot([response_time response_time], [obj.Hf(1)-10 obj.Hf(end)+10],'black','LineWidth',2);               
               hold on;
           end
            
           hold off;
           mtit([condition, ' - channel ', num2str(channel)],'fontsize',20,'color','black');
           % plot mean frequencies of all epochs
           figure('Name','Mean frequencies for all epochs');
           imagesc(correct_epochs); %kreslime jen nevyrazene epochy
           title([condition, ' - channel ', num2str(channel)],'fontsize',20,'color','black');
           colormap parula; %aby to bylo jasne u vsech verzi matlabu - i 2016
           colorbar;
           caxis(zlimits);
        end
        
        function PlotMovingEpochs(obj, channels)
            %pro zadane channels kresli graf cas x frekvence pro kazdou epochu zvlast
            %sipkami se da prochazet pres epochy a kanaly - private function MovePlotEpochs
            %Nada since 2018/01
            assert(~isempty(obj.HFreqEpochs),'soubor s frekvencnimi daty pro epochy neexistuje');
            obj.plotEpochs.f = figure('Name','All Epochs','Position', [20, 100, 1000, 600]);
            obj.plotEpochs.channels = channels;
            set(obj.plotEpochs.f, 'KeyPressFcn', @obj.MovePlotEpochs);
            
            obj.plotEpochs.iChannel = 1; % initiate channel index
            obj.plotEpochs.iEpoch = 1; % initiate epoch index
            obj.plotEpochs.T = linspace(obj.epochtime(1), obj.epochtime(2), size(obj.HFreqEpochs,1)); % time
            obj.plotEpochs.rejectedEpochs = obj.PsyData.GetErrorTrials(); % get rejected epoch trials
            
            % calculate zlimits for all channels
            obj.plotEpochs.zlimits = zeros(channels,2);
            for ch = 1:size(channels)
                obj.plotEpochs.zlimits(ch, :) = obj.getZlimits(obj.plotEpochs.channels(ch));
            end
            
            obj.plotEpochData();
        end    
         
        function plotEpochData(obj)     
            %vykresli vlavo time x frequency graf pre danu epochu 
            %vpravo average power cez vsetky frekvencie danej epochy 
            %pouziva sa v PlotMovingEpochs
            
            subplot(1,2,1) % subplot time x frequency power for given epoch
            imagesc(squeeze(obj.HFreqEpochs(:,obj.plotEpochs.channels(obj.plotEpochs.iChannel),:,obj.plotEpochs.iEpoch))', 'XData', obj.plotEpochs.T, 'YData', obj.Hf);
            caxis(obj.plotEpochs.zlimits(obj.plotEpochs.iChannel,:));
            colormap parula; %aby to bylo jasne u vsech verzi matlabu - i 2016
            set(gca,'YDir','normal');
            colorbar;
            
            hold on; % plot rejected line 
            obj.PlotRejected(obj.plotEpochs.channels(obj.plotEpochs.iChannel), obj.plotEpochs.T, obj.plotEpochs.iEpoch, sum(obj.plotEpochs.rejectedEpochs(obj.plotEpochs.iEpoch,:)));
            response_time = obj.PsyData.P.data(obj.plotEpochs.iEpoch, 4);
            
            hold on; % plot response time 
            plot([response_time response_time], [obj.Hf(1)-10 obj.Hf(end)+10],'black','LineWidth',4);
            hold on;
            
            subplot(1,2,2) % subplot mean power across all frequencies
            plot(obj.plotEpochs.T, obj.d(:,obj.plotEpochs.channels(obj.plotEpochs.iChannel), obj.plotEpochs.iEpoch)');
            
            hold on; % plot response time 
            plot([response_time response_time], obj.plotEpochs.zlimits(obj.plotEpochs.iChannel,:), 'black', 'LineWidth', 4);
            ylim(obj.plotEpochs.zlimits(obj.plotEpochs.iChannel,:)); % y axis = zlimits (default/specified by user)
            xlim([obj.plotEpochs.T(1) obj.plotEpochs.T(end)]); % x axis = time
            
            title(sprintf('%s - channel %d epoch %d', obj.epochData{obj.plotEpochs.iEpoch,1}, obj.plotEpochs.channels(obj.plotEpochs.iChannel), obj.plotEpochs.iEpoch), 'FontSize', 12);
            hold off;
        end
        
        %% SAVE AND LOAD FILE
        %dve funkce na ulozeni a nacteni vypocitane Hilbertovy obalky, protoze to trva hrozne dlouho
        %uklada se vcetne dat parenta CiEEGData
        %trida se musi jmenovat jinak nez v parentovi, protoze jinak se vola tato overloaded function, i z parenta kdyz to nechci
        function obj = Save(obj,filename)
            if ~exist('filename','var')
                filename = obj.hfilename;
                assert( ~isempty(filename), 'no filename given or saved before');
            else
                obj.hfilename = filename;
            end            
            Save@CiEEGData(obj,CHilbert.filenameE(filename));  %ulozim do prvniho souboru data z nadrazene tridy          
            if ~isempty(obj.HFreq)                
                HFreq = obj.HFreq;   %#ok<PROPLC,NASGU>
                Hf = obj.Hf;         %#ok<PROPLC,NASGU> 
                Hfmean = obj.Hfmean; %#ok<PROPLC,NASGU> 
                HFreqEpochs = obj.HFreqEpochs; %#ok<PROPLC,NASGU> %time x channel x frequency x epoch
                yrange = obj.yrange; %#ok<NASGU>
                save(CHilbert.filenameH(filename),'HFreq','Hf','Hfmean','HFreqEpochs','yrange','-v7.3'); %do druheho souboru data z teto tridy
            end
        end
        
        %pokud je treti parametr 1, nenacitaji se data z nadrazene tridy
        function obj = Load(obj,filename,onlyself)            
            if ~exist('onlyself','var') || onlyself == 0
                assert(exist(CHilbert.filenameE(filename),'file')==2, ['soubor s daty neexistuje, mozna se jedna o data tridy CiEEGData?:' char(10) CHilbert.filenameE(filename)]);    
                Load@CiEEGData(obj,CHilbert.filenameE(filename));  
            end
            if exist(CHilbert.filenameH(filename),'file')                
                load(CHilbert.filenameH(filename),'HFreq','Hf','yrange');
                obj.HFreq = HFreq;        %#ok<CPROPLC>
                obj.Hf = Hf;               %#ok<CPROPLC>                 
                obj.yrange = yrange;
                vars = whos('-file',filename);
                if ismember('Hfmean', {vars.name}) %7.4.2017
                    load(filename,'Hfmean');      obj.Hfmean = Hfmean; %#ok<CPROPLC>
                else
                    obj.Hfmean = (obj.Hf(1:end-1) + obj.Hf(2:end)) ./ 2;
                end
                if ismember('HFreqEpochs', {vars.name}) %7.4.2017
                    load(filename,'HFreqEpochs');      obj.HFreqEpochs = HFreqEpochs; %#ok<CPROPLC>
                else
                    obj.HFreqEpochs = [];
                end
            else
                warning(['soubor s frekvencnimi pasmy neexistuje ' CHilbert.filenameH(filename)]);
            end
            obj.hfilename = filename; 
        end 
        function obj = ExtractData(obj,chns,label)
            %ExtractData(obj,chns,filename)
            %vytvori data z vyberu elektrod, pro sdruzeni elektrod pres vsechny pacienty. 
            %pole d, tabs, RjEpochCh a header H
            %jen epochovana data, bipolarni reference
            assert(obj.epochs > 1,'nejsou epochovana data');
            assert(strcmp(obj.reference,'Bipolar'),'neni bipolarni reference');
            d = obj.d(:,chns,:); %#ok<NASGU> %vsechny casy a epochy, vyber kanalu
            tabs = obj.tabs; %#ok<NASGU> %to je spolecne pro vsechny kanaly; time x epochs
            tabs_orig = obj.tabs_orig; %#ok<NASGU> 
            fs = obj.fs; %#ok<NASGU> 
            P = obj.PsyData.P; %#ok<NASGU> %psychopy data
            epochtime = obj.epochtime; %#ok<NASGU> %abych vedel kde je podnet
            baseline = obj.baseline; %#ok<NASGU> 
            RjEpochCh = obj.RjEpochCh(chns,:); %#ok<NASGU> %kanaly vs epochy
            epochData = obj.epochData; %#ok<NASGU> %identita jednotlivych epoch. Musi byt stejna pres pacienty
            DatumCas = obj.DatumCas;
            DatumCas.Extracted = datestr(now);          
            H = obj.CH.H;
            H = rmfield(H,'electrodes'); %smazu nepotrebna pole
            H = rmfield(H,'selCh_H');
            H = rmfield(H,'triggerCH');
            H.channels = H.channels(chns); %vyfiltruju kanaly
            for ch = 1:numel(H.channels)
                H.channels(ch).name = [H.subjName ' ' H.channels(ch).name]; %pridam ke jmenu kanalu jmeno subjektu
            end
            Hf = obj.Hf; %#ok<PROPLC> 
            Hfmean = obj.Hfmean;  %#ok<PROPLC> 
            if isempty(Hfmean), Hfmean = (Hf(1:end-1) + Hf(2:end)) ./ 2; end %#ok<PROPLC,NASGU>             
            HFreq = obj.HFreq(:,chns,:,:); %#ok<PROPLC,NASGU>  %time x channel x freq (x kategorie)
            
            [filepath,fname,~] = fileparts(obj.filename);
            podtrzitko = strfind(fname,'_'); %chci zrusit cast za poslednim podtrzitkem
            filename =[filepath filesep fname(1:podtrzitko(end)-1) ' ' label '_Extract']; 
            save(filename,'d','tabs','tabs_orig','fs','P','epochtime','baseline','RjEpochCh','epochData','DatumCas','H','Hf','Hfmean','HFreq','-v7.3'); 
        end
        
        function CB = ExtractBrainPlotData(obj,chns)
            %vytvori data, pro CBrainPlot.PlotBrain3D
            CB = struct;            
            if ~exist('chns','var'), chns = []; end %pokud neni definovane, je prazdne a pak vytvarim jen d            
            CB.intervals = iff(isempty(chns),[1],[0 1]); %budu mit dve pole hodnoty, vybrane kanaly a vsechny kanaly s vybranymi vyznacenyma
            CB.katstr = iff(isempty(chns),{'all'},{'selected','all'});
            
            CB.VALS = cell(1,iff(isempty(chns),1,2));
            CB.NAMES = cell(1,iff(isempty(chns),1,2));                        
            CB.MNI = cell(1,iff(isempty(chns),1,2));
            CB.LABELS = cell(1,iff(isempty(chns),1,2)); %sem budu ukladata neurologyLabel od Martina Tomaska  
            CB.EPI = cell(1,iff(isempty(chns),1,2)); %pridam jeste udaje o epilepticke aktivite, ktera pak muzu pouzit v zobrazeni mozku            
            
            CB.NAMES{1}= cell(obj.channels,1);            
            CB.MNI{1} = struct('MNI_x',{},'MNI_y',{},'MNI_z',{});            
            CB.VALS{1} = zeros(obj.channels,1); %tam pak doplnim 1 na mista vybranych kanalu         
            CB.LABELS{1}= cell(obj.channels,1);
            CB.EPI{1} = struct('seizureOnset',{},'interictalOften',{},'rejected',{});
            if ~isempty(chns)
                CB.NAMES{2}= cell(numel(chns),1);
                CB.VALS{2} = ones(numel(chns),1);
                CB.MNI{2} = struct('MNI_x',{},'MNI_y',{},'MNI_z',{});                        
                CB.LABELS{2}= cell(numel(chns),1);
                CB.EPI{2} = struct('seizureOnset',{},'interictalOften',{},'rejected',{});
            end
      
            for ch = 1:obj.channels
                CB.NAMES{1}{ch} = obj.CH.H.channels(ch).name;
                CB.LABELS{1}{ch} = obj.CH.H.channels(ch).neurologyLabel;
                CB.MNI{1}(ch).MNI_x = obj.CH.H.channels(ch).MNI_x;
                CB.MNI{1}(ch).MNI_y = obj.CH.H.channels(ch).MNI_y;
                CB.MNI{1}(ch).MNI_z = obj.CH.H.channels(ch).MNI_z;
                if isfield(obj.CH.H.channels,'seizureOnset')
                    CB.EPI{1}(ch).seizureOnset = obj.CH.H.channels(ch).seizureOnset;
                    CB.EPI{1}(ch).interictalOften = obj.CH.H.channels(ch).interictalOften;
                    CB.EPI{1}(ch).rejected = obj.CH.H.channels(ch).rejected;
                end
            end            
                  
            for ch = 1:numel(chns)
                CB.VALS{1}(chns(ch)) = 1;
                CB.NAMES{2}{ch} = obj.CH.H.channels(chns(ch)).name;     
                CB.LABELS{2}{ch} = obj.CH.H.channels(chns(ch)).neurologyLabel; 
                CB.MNI{2}(ch).MNI_x = obj.CH.H.channels(chns(ch)).MNI_x;
                CB.MNI{2}(ch).MNI_y = obj.CH.H.channels(chns(ch)).MNI_y;
                CB.MNI{2}(ch).MNI_z = obj.CH.H.channels(chns(ch)).MNI_z;
                if isfield(obj.CH.H.channels,'seizureOnset')
                    CB.EPI{2}(ch).seizureOnset = obj.CH.H.channels(chns(ch)).seizureOnset;
                    CB.EPI{2}(ch).interictalOften = obj.CH.H.channels(chns(ch)).interictalOften;
                    CB.EPI{2}(ch).rejected = obj.CH.H.channels(chns(ch)).rejected;
                end
            end
        end
    end 
        
    %  --------- privatni metody ----------------------
    methods (Static,Access = private)
        function [ freqPow ] = hilbertJirka(rawData, loF, hiF, srate )                
            %HILBERTJIRKA hilbertJirka( rawData, loF, hiF, srate )
            %   vrati Power vybraneho frekvencniho pasma

            % Hilberova obalka podle Jirky Hammera - mail 1.8.2014
            %loF = pasmo; hiF = pasmo + 10; srate = 1000;
            %rawData = yy;

            %1) filrovani v gamma pasmu: loF = 60, hiF = 100, srate = sampling rate
            Wn = [loF, hiF]/(srate/2); % normalized bandpass frequencies
            n = 4; % butterworth order
            [b,a] = butter(n, Wn); % returns polynoms of Butterw. filter
            filtData = filtfilt(b, a, rawData);

             %2) analyticky signal pomoci Hilbertovy transformace:
            tmp = hilbert(filtData);

             %3) gamma amplitude, power
            %gammaAmp = abs(tmp);
            freqPow = abs(tmp).^2; %power je druha mocnicna, 
            % viz https://en.wikipedia.org/wiki/Spectral_density 
            %-  "power" is simply reckoned in terms of the square of the signal,
        end
    end
    methods (Static,Access = public)
        function filename2 = filenameE(filename)
            %vraci jmeno souboru s daty tridy CiEEGData
           filename=strrep(filename,'_CHilb',''); %odstranim pripony vytvorene pri save
           filename=strrep(filename,'_CiEEG','');
           [pathstr,fname,ext] = fileparts(filename);  
           if strcmp(ext,'.mat')==false || numel(ext)<1
               fname = [fname ext]; %pokud pripona neni mat, pridam ji na konec jmena a vytvorim priponu mat
               ext = '.mat';
           end           
           filename2 = fullfile(pathstr,[fname '_CiEEG' ext]);
        end
        function filename2 = filenameH(filename)
             %vraci jmeno souboru s daty teto tridy
           filename=strrep(filename,'_CHilb',''); %odstranim pripony vytvorene pri save
           filename=strrep(filename,'_CiEEG','');
           [pathstr,fname,ext] = fileparts(filename);
           if strcmp(ext,'.mat')==false || numel(ext)<1
               fname = [fname ext]; %pokud pripona neni mat, pridam ji na konec jmena a vytvorim priponu mat
               ext = '.mat';
           end            
           filename2 = fullfile(pathstr,[fname '_CHilb' ext]);
        end
    end
    methods  (Access = private)
        function obj = hybejPlotF(obj,~,eventDat)  
           %reaguje na udalosti v grafu PlotResponseCh
           switch eventDat.Key
               case 'rightarrow' 
                   obj.PlotResponseFreq( min( [obj.plotF.ch + 1 , obj.channels]));
               case 'pagedown' 
                   obj.PlotResponseFreq( min( [obj.plotF.ch + 10 , obj.channels]));
               case 'leftarrow'
                   obj.PlotResponseFreq( max( [obj.plotF.ch - 1 , 1]));
               case 'pageup'
                   obj.PlotResponseFreq( max( [obj.plotF.ch - 10 , 1]));
               case 'space' %zobrazi i prumerne krivky
                   obj.PlotResponseCh(obj.plotF.ch);
                   figure(obj.plotF.fh); %dam puvodni obrazek dopredu
               case {'multiply','8'} %hvezdicka na numericke klavesnici
                   %dialog na vlozeni minima a maxima osy y
                   answ = inputdlg('Enter ymax and min:','Yaxis limits', [1 50],{num2str(obj.plotF.ylim)});
                   if numel(answ)>0  %odpoved je vzdy cell 1x1 - pri cancel je to cell 0x0
                       if isempty(answ{1}) || any(answ{1}=='*') %pokud vlozim hvezdicku nebo nic, chci znovy spocitat max a min
                           obj.plotF.ylim = [];
                       else %jinak predpokladam dve hodnoty
                           data = str2num(answ{:});  %#ok<ST2NM>
                           if numel(data)>= 2 %pokud nejsou dve hodnoty, nedelam nic
                             obj.plotF.ylim = [data(1) data(2)];
                           end
                       end
                   end
                   obj.PlotResponseFreq( obj.plotF.ch); %prekreslim grafy
                case 'delete' %Del na numericke klavesnici
                   obj.plotF.ylim = [];
                   obj.PlotResponseFreq( obj.plotEp.ch); %prekreslim grafy
           end
        end
        
        function obj = MovePlotEpochs(obj,~,eventDat)
            %zpracovava stlaceni klavesy pro graf PlotMovingEpochs
            switch eventDat.Key
                case 'rightarrow' % +1 epoch
                    obj.plotEpochs.iEpoch = min([obj.plotEpochs.iEpoch + 1, size(obj.HFreqEpochs,4)]);
                case 'leftarrow'  % -1 epoch
                    obj.plotEpochs.iEpoch = max([obj.plotEpochs.iEpoch - 1, 1]);
                case 'uparrow'    % -1 channel
                    obj.plotEpochs.iChannel = max([obj.plotEpochs.iChannel - 1, 1]);
                case 'downarrow'  % +1 channel
                    obj.plotEpochs.iChannel = min([obj.plotEpochs.iChannel + 1, length(obj.plotEpochs.channels)]);
                case 'numpad6' % nasledujuca s rovnakou condition
                    obj.plotEpochs.iEpoch = min([obj.getNextCondition(1), size(obj.HFreqEpochs,4)]);
                case 'numpad4' % predchadzajuca s rovnakou condition
                    obj.plotEpochs.iEpoch = max([obj.getLastCondition(1), 1]);
                case 'pageup' % nasledujuca condition
                    obj.plotEpochs.iEpoch = min([obj.getNextCondition(0), size(obj.HFreqEpochs,4)]);
                case 'pagedown' % predchadzajuca condition
                    obj.plotEpochs.iEpoch = max([obj.getLastCondition(0), 1]);
                case {'multiply','8'} %hvezdicka na numericke klavesnici
                   %dialog na vlozeni minima a maxima osy y
                   answ = inputdlg('Enter ymax and min:','Yaxis limits', [1 50], {num2str(obj.plotEpochs.zlimits(obj.plotEpochs.iChannel,:))});
                   if numel(answ)>0  %odpoved je vzdy cell 1x1 - pri cancel je to cell 0x0
                       if isempty(answ{1}) || any(answ{1}=='*') %pokud vlozim hvezdicku nebo nic, chci znovy spocitat max a min
                           obj.plotEpochs.zlimits(obj.plotEpochs.iChannel,:) = obj.getZlimits(obj.plotEpochs.channels(obj.plotEpochs.iChannel));
                       else %jinak predpokladam dve hodnoty
                           data = str2num(answ{:});  %#ok<ST2NM>
                           if numel(data)>= 2 %pokud nejsou dve hodnoty, nedelam nic
                             obj.plotEpochs.zlimits(obj.plotEpochs.iChannel,:) = [data(1) data(2)];
                           end
                       end
                   end
                   obj.plotEpochData(); %prekreslim grafy
                otherwise  
                   display(['key pressed: ' eventDat.Key]); %vypise stlacenou klavesu
            end
            obj.plotEpochData();
        end
        
        function correct = PlotRejected(obj, channel, time, epoch, sum_chyby)
           %vraci true, pokud epocha neni oznacena jako vyrazena
           %pouziva se v PlotAllEpochs
           correct = false;           
           if obj.PsyData.P.data(epoch, 6) % mark trening answers with red
               hold on; plot([time(1) time(end)], [obj.Hf(end) obj.Hf(1)],'red','LineWidth',6)
           
           elseif ~obj.PsyData.P.data(epoch, 3) || sum_chyby > 0 % mark wrong answers with black
               hold on; plot([time(1) time(end)], [obj.Hf(end) obj.Hf(1)],'black','LineWidth',6)
           
           elseif ismember(epoch, obj.RjEpoch) % mark rejected epochs with green
               hold on; plot([time(1) time(end)], [obj.Hf(end) obj.Hf(1)],'green','LineWidth',6)
           
           elseif obj.RjEpochCh(channel, epoch) % mark rejected channel's epochs with blue
               hold on; plot([time(1) time(end)], [obj.Hf(end) obj.Hf(1)],'blue','LineWidth',6)
           else
              correct = true;
           end
        end
        
        function zlimits = getZlimits(obj, ch)
            %vypocita minimalnu a maximalnu power pre dany channel cez
            %vsetky epochy a frekvencie
            %pouziva sa v PlotMovingEpochs
            ymin = min(min(obj.d(:, ch, :)));
            ymax = max(max(obj.d(:, ch, :)));
            zlimits = [ymin ymax];
        end
        
        function last = getLastCondition(obj, same)
            %najde index poslednej najblizsej epochy 
            %rovnakej(same=1)/rozdielnej(same=0) kategorie
            %pouziva sa v PlotMovingEpochs pri numpad4/pagedown
            if same
                last = find(obj.PsyData.P.data(1:(obj.plotEpochs.iEpoch-1),7) == obj.epochData{obj.plotEpochs.iEpoch,2}, 1, 'last');
            else 
                last = find(obj.PsyData.P.data(1:(obj.plotEpochs.iEpoch-1),7) ~= obj.epochData{obj.plotEpochs.iEpoch,2}, 1, 'last');
            end  
            if isempty(last)
                last = obj.plotEpochs.iEpoch;
            end
        end
        
        function next = getNextCondition(obj, same)
            %najde index najblizsej epochy 
            %rovnakej(same=1)/rozdielnej(same=0) kategorie
            %pouziva sa v PlotMovingEpochs pri numpad6/pageup
            if same
                next = obj.plotEpochs.iEpoch + find(obj.PsyData.P.data((obj.plotEpochs.iEpoch+1):end,7) == obj.epochData{obj.plotEpochs.iEpoch,2}, 1);
            else 
                next = obj.plotEpochs.iEpoch + find(obj.PsyData.P.data((obj.plotEpochs.iEpoch+1):end,7) ~= obj.epochData{obj.plotEpochs.iEpoch,2}, 1);
            end
            if isempty(next)
                next = obj.plotEpochs.iEpoch;
            end
        end
       
    end
end

