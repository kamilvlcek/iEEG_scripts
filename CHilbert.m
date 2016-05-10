classdef CHilbert < CiEEGData
    %HILBERT.CLASS sbirka funkci na analyzu pomoci hilbert trasform
    %   Kamil Vlcek, 2016-04
    properties (Constant = true)
        decimatefactor = 8; %o kolik decimuju hiblertovu obalku oproti puvodi sampling rate; 2 je dostatecne konzervativni, 8 hodne setri pamet
    end
    properties (Access = public)
        HFreq; %hilberova obalka pro kazde frekvenci pasmo - time x channel x freq
        HMean; %prumerna hilbertova obalka za frekvence
        Hfs; %sampling frekvence hilberovy obalky
        Htabs; %decimovane originaly tabs
        Hf; %frekvencni pasma pro ktere jsou pocitany obalky
    end
    methods (Access = public)
        function obj = CHilbert(d,tabs,fs,mults,header)
            if ~exist('header','var'), header = []; end %nejakou hodnotu dat musim
            if ~exist('mults','var'),  mults = []; end %nejakou hodnotu dat musim
            obj@CiEEGData(d,tabs,fs,mults,header);
        end
        
        function obj = PasmoFrekvence(obj,freq)
            %EEG2HILBERT prevede vsechny kanaly na prumer hilbertovych obalek
            %   podle puvodni funkce EEG2Hilbert
            %   pouziva data d z parentu a take fs
            %   freq    seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi
            
            obj.HFreq = zeros(ceil(obj.samples/obj.decimatefactor),obj.channels,numel(freq)-1); %inicializace pole   
            fprintf('kanal ze %i: ', obj.channels);
            for ch = 1:obj.channels %jednotlive elektrody
                %fprintf('channel %i: Hz ',ch);
                           
                fprintf('%i,',ch);
                for fno = 1:numel(freq)-1 %seznam frekvenci
                    loF = freq(fno);
                    hiF = freq(fno+1)-1; 
                    hh = obj.hilbertJirka(obj.DData(ch),loF,hiF,obj.fs); %cista hilbertova obalka, tohle i skript hodne zrychli
                    hh = decimate(hh,obj.decimatefactor); %dvakrat mensi sampling rate                    
                    obj.HFreq(:,ch,fno) = (hh./mean(hh)).* 100; %podil prumeru 100 = prumerna hodnota
                    %fprintf('%i Hz, ',loF);
                end
                %fprintf('\n'); %tisk znova na stejnou radku
            end
            obj.HMean = squeeze(mean(obj.HFreq,3));
            obj.Hfs = obj.fs/obj.decimatefactor;
            obj.Htabs = downsample(obj.tabs,obj.decimatefactor);
            obj.Hf = freq;
            fprintf('\n'); %ukoncim radku
        end
        
        %dve funkce na ulozeni a nacteni vypocitane Hilbertovy obalky, protoze to trva hrozne dlouho
        %ukladaji pouze vysledky funkce PasmoFrekvence
        function SaveH(obj,filename)
            HFreq = obj.HFreq; %#ok<PROP,NASGU>
            HMean = obj.HMean; %#ok<PROP,NASGU>
            Hfs = obj.Hfs;     %#ok<PROP,NASGU>
            Htabs = obj.Htabs; %#ok<PROP,NASGU>
            Hf = obj.Hf;       %#ok<PROP,NASGU>
            save(filename,'HFreq','HMean','Hfs','Htabs','Hf','-v7.3');
        end
        function obj = LoadH(obj,filename)
            load(filename,'HFreq','HMean','Hfs','Htabs','Hf');
            obj.HFreq = HFreq; %#ok<CPROP,PROP>
            obj.HMean = HMean; %#ok<CPROP,PROP>
            obj.Hfs = Hfs;     %#ok<CPROP,PROP>
            obj.Htabs = Htabs; %#ok<CPROP,PROP>
            obj.Hf = Hf;        %#ok<CPROP,PROP>
        end
        
        function ExtractEpochs(obj, PsyData,epochtime)
            % rozdeli hilbertovu obalku podle epoch
            % i u ni odecte baseline pred podnetem
             ExtractEpochs@CiEEGData(obj,PsyData, epochtime);
             iepochtime = round(epochtime.*obj.Hfs); %v poctu vzorku cas pred a po udalosti, prvni cislo je zaporne druhe kladne
             Hmean2 = zeros(iepochtime(2)-iepochtime(1), size(obj.Hmean,2), size(obj.epochData,1)); %nova epochovana data time x channel x epoch  
             for epoch = 1:size(obj.epochData,1) 
                 izacatek = find(obj.Htabs==obj.epochData{epoch,2}); %najdu index podnetu, podle jeho timestampu. v druhem sloupci epochData jsou timestampy
                 for ch=1:obj.channels
                     baseline = mean(obj.Hmean(izacatek + iepochtime(1) : izacatek-1,ch)); %prumerna hodnota v case pred podnetem
                     Hmean2(:,ch,epoch)=  obj.Hmean(izacatek + iepochtime(1) : izacatek+iepochtime(2)-1,ch) - baseline; %i u Hilbertovy obalky odecitam baseline
                 end
             end
             obj.Hmean = Hmean2; %nahradim puvodni prumernou hilbertovu obalku tou epochovanou             
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
    
end

