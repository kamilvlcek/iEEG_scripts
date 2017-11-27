classdef CHilbert < CiEEGData
    %HILBERT.CLASS sbirka funkci na analyzu pomoci hilbert trasform
    %   Kamil Vlcek, FGU AVCR, since 2016 04
    properties (Constant = true)
        decimatefactor = 8; %o kolik decimuju hiblertovu obalku oproti puvodi sampling rate; 2 je dostatecne konzervativni, 8 hodne setri pamet
    end
    properties (Access = public)
        HFreq; %hilberova obalka pro kazde frekvenci pasmo - time x channel x freq (x kategorie)
        Hf; %frekvencni pasma pro ktere jsou pocitany obalky - okraje pasem, pocet je tedy vetsi o 1 nez pocet spocitanych pasem
        Hfmean; %stredni hodnoty pasem  - pocet = pocet spocitanych pasem
        hfilename; %jmeno souboru CHilbert  
        plotF = struct; %udaje o stavu plotu PlotResponseFreq
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
        
        function obj = PasmoFrekvence(obj,freq,channels,prekryv)
            %EEG2HILBERT prevede vsechny kanaly na prumer hilbertovych obalek
            %   podle puvodni funkce EEG2Hilbert
            %   pouziva data d z parentu a take fs
            %   freq    seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi
            if ~exist('channels','var') || isempty(channels), channels = 1:obj.channels; end
            if ~exist('prekryv','var'), prekryv = 0; end %kolik se maji prekryvat sousedni frekvencni pasma, napriklad 0.5
            
            obj.HFreq = zeros(ceil(obj.samples/obj.decimatefactor),obj.channels,numel(freq)-1); %inicializace pole   
            tic; %zadnu meric cas
            fprintf('kanal ze %i: ', max(channels) );
            for ch = channels %jednotlive elektrody
                %fprintf('channel %i: Hz ',ch);
                           
                fprintf('%i,',ch);
                for fno = 1:numel(freq)-1 %seznam frekvenci
                    loF = freq(fno) -prekryv*(freq(fno+1)-freq(fno)); 
                    hiF = freq(fno+1)-0.1 +prekryv*(freq(fno+1)-freq(fno));  %napr 50 - 59.9
                    hh = obj.hilbertJirka(obj.d(:,ch),loF,hiF,obj.fs); %cista hilbertova obalka, tohle i skript hodne zrychli
                    hh = decimate(hh,obj.decimatefactor); % mensi sampling rate                    
                    obj.HFreq(:,ch,fno) = (hh./mean(hh)); %podil prumeru = prumerna hodnota
                    %fprintf('%i Hz, ',loF);
                end
                %fprintf('\n'); %tisk znova na stejnou radku
            end
            obj.d = squeeze(mean(obj.HFreq,3)); %11.5.2016 - prepisu puvodni data prumerem
            obj.fs = obj.fs/obj.decimatefactor;
            obj.tabs = downsample(obj.tabs,obj.decimatefactor);
            obj.tabs_orig = downsample(obj.tabs_orig,obj.decimatefactor); %potrebuju zdecimovat i druhy tabs. Orig znamena jen ze nepodleha epochovani
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
                 %cyklus po kategoriich ne po epochach
                 for katnum = kategorie' %potrebuji to v radcich
                     Epochy = find(cell2mat(obj.epochData(:,2))==katnum); %seznam epoch v ramci kategorie ve sloupci 
                     for epoch = Epochy' %potrebuji to v radcich
                         izacatek = find(obj.tabs_orig==obj.epochData{epoch,3}); %najdu index podnetu, podle jeho timestampu. v tretim sloupci epochData jsou timestampy
                         for ch=1:obj.channels
                            baseline_mean = mean(obj.HFreq(izacatek + ibaseline(1) : izacatek+ibaseline(2)-1,ch,:),1); %baseline pro vsechny frekvencni pasma dohromady
                            Hfreq2(:,ch,:,katnum+1) = Hfreq2(:,ch,:,katnum+1) +  bsxfun(@minus,obj.HFreq(izacatek + iepochtime(1) : izacatek+iepochtime(2)-1,ch,:) , baseline_mean); 
                            %tady se mi to mozna odecetlo blbe? KOntrola
                         end
                     end
                     Hfreq2(:,:,:,katnum+1) = Hfreq2(:,:,:,katnum+1)./numel(Epochy);
                 end             
                 obj.HFreq = Hfreq2;
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
                HFreq = obj.HFreq;  %#ok<PROP,NASGU>
                Hf = obj.Hf;         %#ok<PROP,NASGU>  
                Hfmean = obj.Hfmean; %#ok<PROP,NASGU>  
                yrange = obj.yrange; %#ok<NASGU> 
                save(CHilbert.filenameH(filename),'HFreq','Hf','Hfmean','yrange','-v7.3'); %do druheho souboru data z teto tridy
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
                obj.HFreq = HFreq;  %#ok<CPROP,PROP>            
                obj.Hf = Hf;                %#ok<CPROP,PROP>
                obj.yrange = yrange;
                vars = whos('-file',filename);
                if ismember('Hfmean', {vars.name}) %7.4.2017
                    load(filename,'Hfmean');      obj.Hfmean = Hfmean; %#ok<CPROP,PROP>
                else
                    obj.Hfmean = (obj.Hf(1:end-1) + obj.Hf(2:end)) ./ 2;
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
            d = obj.d(:,chns,:); %vsechny casy a epochy, vyber kanalu
            tabs = obj.tabs; %to je spolecne pro vsechny kanaly; time x epochs
            tabs_orig = obj.tabs_orig;
            fs = obj.fs;
            P = obj.PsyData.P; %psychopy data
            epochtime = obj.epochtime; %abych vedel kde je podnet
            RjEpochCh = obj.RjEpochCh(chns,:); %kanaly vs epochy
            epochData = obj.epochData; %identita jednotlivych epoch. Musi byt stejna pres pacienty
            DatumCas = obj.DatumCas;
            DatumCas.Extracted = datestr(now);
            H = obj.CH.H;
            H = rmfield(H,'electrodes'); %smazu nepotrebna pole
            H = rmfield(H,'selCh_H');
            H = rmfield(H,'triggerCH');
            H.channels = H.channels(chns); %vyfiltruju kanaly
            [filepath,fname,~] = fileparts(obj.filename);
            podtrzitko = strfind(fname,'_'); %chci zrusit cast za poslednim podtrzitkem
            filename =[filepath filesep fname(1:podtrzitko(end)-1) ' ' label '_Extract']; 
            save(filename,'d','tabs','tabs_orig','fs','P','epochtime','RjEpochCh','epochData','DatumCas','H','-v7.3'); 
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
           if numel(ext)<1, ext = '.mat'; end
           filename2 = fullfile(pathstr,[fname '_CiEEG' ext]);
        end
        function filename2 = filenameH(filename)
             %vraci jmeno souboru s daty teto tridy
           filename=strrep(filename,'_CHilb',''); %odstranim pripony vytvorene pri save
           filename=strrep(filename,'_CiEEG','');
           [pathstr,fname,ext] = fileparts(filename); 
           if numel(ext)<1, ext = '.mat'; end
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
    end
end

