classdef CHilbert < CiEEGData
    %HILBERT.CLASS sbirka funkci na analyzu pomoci hilbert trasform
    %   Kamil Vlcek, FGU AVCR, since 2016 04
    properties (Constant = true)
        decimatefactor = 8; %o kolik decimuju hiblertovu obalku oproti puvodi sampling rate; 2 je dostatecne konzervativni, 8 hodne setri pamet
    end
    properties (Access = public)
        HFreq; %hilberova obalka pro kazde frekvenci pasmo - time x channel x freq (x kategorie)
        Hf; %frekvencni pasma pro ktere jsou pocitany obalky
        hfilename; %jmeno souboru CHilbert        
    end
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS 
        function obj = CHilbert(d,tabs,fs,mults,header)
            if ~exist('header','var'), header = []; end %nejakou hodnotu dat musim
            if ~exist('mults','var'),  mults = []; end %nejakou hodnotu dat musim
            if ischar(d) && ~exist('tabs','var') %pokud je prvni parametr retezec, tak ho beru jako nazev souboru, ktery nactu
                tabs=[]; fs = [];
                filename = d;
                d = CHilbert.filenameE(d); %predelam filename na filename pro CiEEGData                
            end
            obj@CiEEGData(d,tabs,fs,mults,header); %volani konstruktoru nemuze byt v if bloku ! Sakra
            if exist('filename','var')
                obj.LoadH(filename,1); %zkusim nacist jeste data z teto tridy
            end
            
        end
        
        function obj = PasmoFrekvence(obj,freq,channels)
            %EEG2HILBERT prevede vsechny kanaly na prumer hilbertovych obalek
            %   podle puvodni funkce EEG2Hilbert
            %   pouziva data d z parentu a take fs
            %   freq    seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi
            if ~exist('channels','var'), channels = 1:obj.channels; end
            obj.HFreq = zeros(ceil(obj.samples/obj.decimatefactor),obj.channels,numel(freq)-1); %inicializace pole   
            fprintf('kanal ze %i: ', obj.channels);
            for ch = channels %jednotlive elektrody
                %fprintf('channel %i: Hz ',ch);
                           
                fprintf('%i,',ch);
                for fno = 1:numel(freq)-1 %seznam frekvenci
                    loF = freq(fno);
                    hiF = freq(fno+1)-1; 
                    hh = obj.hilbertJirka(obj.DData(ch),loF,hiF,obj.fs); %cista hilbertova obalka, tohle i skript hodne zrychli
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
            obj.mults = ones(1,size(obj.d,2)); %nove pole uz je double defaultove jednicky pro kazdy kanal
            obj.yrange = [1 1 5 5]; %zmenim rozliseni osy y v grafu
            fprintf('\n'); %ukoncim radku
        end
        
        function ExtractEpochs(obj, PsyData,epochtime)
            % rozdeli hilbertovu obalku podle epoch
            % i u ni odecte baseline pred podnetem
            
            ExtractEpochs@CiEEGData(obj,PsyData, epochtime); %to mi zepochuje prumernou obalku za frekvencni pasma v poli d
            if(numel(obj.HFreq)>0)
                %ted epochace vsech frekvencnich pasem zvlast, hlavne kvuli obrazkum
                %prumer za kazdou kategorii, statistiku z toho delat nechci
                 iepochtime = round(epochtime.*obj.fs); %v poctu vzorku cas pred a po udalosti, prvni cislo je zaporne druhe kladne             
                 kategorie = cell2mat(obj.PsyData.P.strings.podminka(:,2)); %cisla karegorii ve sloupcich
                 Hfreq2 = zeros(iepochtime(2)-iepochtime(1), size(obj.d,2), numel(obj.Hf)-1,size(kategorie,1)); %nova epochovana data time x channel x freq x kategorie=podminka
                 %cyklus po kategoriich ne po epochach
                 for katnum = kategorie' %potrebuji to v radcich
                     Epochy = find(cell2mat(obj.epochData(:,2))==katnum); %seznam epoch v ramci kategorie ve sloupci 
                     for epoch = Epochy' %potrebuji to v radcich
                         izacatek = find(obj.tabs_orig==obj.epochData{epoch,3}); %najdu index podnetu, podle jeho timestampu. v tretim sloupci epochData jsou timestampy
                         for ch=1:obj.channels
                            baseline = mean(obj.HFreq(izacatek + iepochtime(1) : izacatek-1,ch,:),1); %baseline pro vsechny frekvencni pasma dohromady
                            Hfreq2(:,ch,:,katnum+1) = Hfreq2(:,ch,:,katnum+1) +  bsxfun(@minus,obj.HFreq(izacatek + iepochtime(1) : izacatek+iepochtime(2)-1,ch,:) , baseline); 
                            %tady se mi to mozna odecetlo blbe? KOntrola
                         end
                     end
                     Hfreq2(:,:,:,katnum+1) = Hfreq2(:,:,:,katnum+1)./numel(Epochy);
                 end             
                 obj.HFreq = Hfreq2;
            end
        end
        
        %% SAVE AND LOAD FILE
        %dve funkce na ulozeni a nacteni vypocitane Hilbertovy obalky, protoze to trva hrozne dlouho
        %uklada se vcetne dat parenta CiEEGData
        %trida se musi jmenovat jinak nez v parentovi, protoze jinak se vola tato overloaded function, i z parenta kdyz to nechci
        function SaveH(obj,filename)
            if ~exist('filename','var')
                filename = obj.hfilename;
                assert( ~isempty(filename), 'no filename given or saved before');
            else
                obj.hfilename = filename;
            end            
            obj.Save(CHilbert.filenameE(filename));  %ulozim do prvniho souboru data z nadrazene tridy          
            if ~isempty(obj.HFreq)                
                HFreq = obj.HFreq;  %#ok<PROP,NASGU>
                Hf = obj.Hf;         %#ok<PROP,NASGU>           
                yrange = obj.yrange; %#ok<NASGU> 
                save(CHilbert.filenameH(filename),'HFreq','Hf','yrange','-v7.3'); %do druheho souboru data z teto tridy
            end
        end
        
        %pokud je treti parametr 1, nenacitaji se data z nadrazene tridy
        function obj = LoadH(obj,filename,onlyself)            
            if ~exist('onlyself','var') || onlyself == 0
                obj.Load(obj,CHilbert.filenameE(filename));            
            end
            if exist(CHilbert.filenameH(filename),'file')                
                load(CHilbert.filenameH(filename),'HFreq','Hf','yrange');
                obj.HFreq = HFreq;  %#ok<CPROP,PROP>            
                obj.Hf = Hf;                %#ok<CPROP,PROP>
                obj.yrange = yrange;       
            end
            obj.hfilename = filename; 
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
           [pathstr,fname,ext] = fileparts(filename); 
           filename2 = fullfile(pathstr,[fname '_CiEEG' ext]);
        end
        function filename2 = filenameH(filename)
             %vraci jmeno souboru s daty teto tridy
           [pathstr,fname,ext] = fileparts(filename); 
           filename2 = fullfile(pathstr,[fname '_CHilb' ext]);
        end
    end
    
end

