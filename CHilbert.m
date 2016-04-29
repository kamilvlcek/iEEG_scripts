classdef CHilbert < CiEEGData
    %HILBERT.CLASS sbirka funkci na analyzu pomoci hilbert trasform
    %   Kamil Vlcek, 2016-04
    
    properties (Access = public)
        HFreq; %hilberova obalka pro kazde frekvenci pasmo - time x channel x freq
        HMean; %prumerna hilbertova obalka za frekvence
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
            
            obj.HFreq = zeros(obj.samples,obj.channels,numel(freq)-1); %inicializace pole   
            fprintf('kanal ze %i: ', obj.channels);
            for ch = 1:obj.channels %jednotlive elektrody
                %fprintf('channel %i: Hz ',ch);
                           
                fprintf('%i,',ch);
                for fno = 1:numel(freq)-1 %seznam frekvenci
                    loF = freq(fno);
                    hiF = freq(fno+1)-1; 
                    hh = obj.hilbertJirka(obj.DData(ch),loF,hiF,obj.fs); %cista hilbertova obalka 
                    obj.HFreq(:,ch,fno) = (hh./mean(hh)).* 100; %podil prumeru 100 = prumerna hodnota
                    %fprintf('%i Hz, ',loF);
                end
                %fprintf('\n'); %tisk znova na stejnou radku
            end
            obj.HMean = squeeze(mean(obj.HFreq,3));
            fprintf('\n'); %ukoncim radku
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
            freqPow = abs(tmp).^2;
        end
    end
    
end

