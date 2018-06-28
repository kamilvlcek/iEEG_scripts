classdef  CMorlet < CHilbert
    properties (Access = public)
        freal; % original filtered data
    end    
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS 
        function obj = CMorlet(d,tabs,fs,mults,header)
            if ~exist('header','var'), header = []; end %nejakou hodnotu dat musim
            if ~exist('mults','var'),  mults = []; end %nejakou hodnotu dat musim
            if ischar(d) && ~exist('tabs','var') %pokud je prvni parametr retezec, tak ho beru jako nazev souboru, ktery nactu
                tabs=[]; fs = [];
                % volani Load z CiEEGData mi zavola Load z CHilbert, takze d=filename predelavat nemusim
            end
            obj@CHilbert(d,tabs,fs,mults,header); %volani konstruktoru nemuze byt v if bloku 
            disp('vytvoren objekt CMorlet');             
        end
        
        function obj = PasmoFrekvence(obj,freq,channels,~,decimatefactor)
            %   podle figure 1311 z knihy Mike X Cohen 
            %   pouziva data d z parentu a take fs
            %   freq    seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi
            if ~exist('channels','var') || isempty(channels) , channels = 1:obj.channels; end
            %if ~exist('prekryv','var'), prekryv = 0; end %kolik se maji prekryvat sousedni frekvencni pasma, napriklad 0.5
            if ~exist('decimatefactor','var') || isempty(decimatefactor) , decimatefactor = obj.decimatefactor; end; %volitelny parametr decimatefactor 
            HFreq = zeros(ceil(obj.samples/decimatefactor),obj.channels,numel(freq)); %inicializace pole - power
            fphase = zeros(ceil(obj.samples/decimatefactor),obj.channels,numel(freq)); % inicializace pole - faze
            ffreal = zeros(ceil(obj.samples/decimatefactor),obj.channels,numel(freq)); % inicializace pole - filtrovana puvodni data
            timer = tic; %zacnu merit cas
            fprintf('kanal ze %i: ', numel(channels) );
            
            if freq(1) < 10
                time = -2:1/obj.fs:2; %pro nizke frekvence pouziju delsi wavelet - +- 2s - % 40 cyklu wavelet pro 10Hz
            else
                time = -1:1/obj.fs:1; %pro vyssi frekvence staci wv +- 1s %alespon 20 cyklu waveletu pro 10Hz. Neni to zbytecne moc ?
            end
            s    =  5./(2*pi*freq);   %gaussian width (or stdev) - 5=pet cyklu waveletu - pro kazdou frekvenci zvlas
            n_wavelet            = length(time); %delka waveletu ve vzorcich - 1025 pro 512 a 1s
            n_data               = obj.samples; %data nejsou epochovana
            n_convolution        = n_wavelet+n_data-1; %pocet vzorku vysledne konvoluce - vic bodu nez maji data, takze budu pridavat 0 na konce
            n_conv_pow2          = pow2(nextpow2(n_convolution)); %prevede pocet vzorku konvoluce na nejblizsi vyssi mocninu 2 
            half_of_wavelet_size = (n_wavelet-1)/2;
            
            for ch = channels %jednotlive elektrody                                           
                fprintf('%i,',ch);
                if sum(obj.d(:,ch))==0, continue; end %pro vyrazene kanaly jsou hodnoty 0 pri jine nez bipol ref. Z tech pak vznikne nan, pri tomhle cyklu, coz vadi dal                                  
                eegfft = fft(obj.d(:,ch)',n_conv_pow2); %FFT of eeg data, potrebuju to dat do radku aby stejne jak-o wavelet                
                for fno = 1:numel(freq) %seznam frekvenci                    
                    wavelet = fft( sqrt(1/(s(fno)*sqrt(pi))) * exp(1i*2*pi*freq(fno).*time) .* exp(-time.^2./(2*(s(fno)^2))) , n_conv_pow2 );
                    % fft ( (A=frequency band-specific scaling factor) * complex sin * gaussian )
                    
                    eegconv = ifft(wavelet.*eegfft); %convoluce pomoci convolution theorem
                    eegconv = eegconv(1:n_convolution);
                    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);    
                    fpower = eegconv .* conj(eegconv); % nasobeni conj je pry 2x rychlejsi  - abs(eegconv).^2; %z komplexniho cisla vypocitam power
                    if decimatefactor > 1
                        fpower = decimate(fpower,decimatefactor); % mensi sampling rate (moving average vubec nepomohl)
                    end
                    HFreq(:,ch,fno) = fpower; % povodna normalizacia (fpower./mean(fpower)) premiestnena do funkcie CHilbert.Normalize                   
                    %fprintf('%i Hz, ',loF);
                    fphase0 = angle(eegconv); %#ok<PROPLC> %faze frekvence 
                    if decimatefactor > 1                   
                        fphase0 = decimate(fphase0,decimatefactor);    %#ok<PROPLC>                        
                    end
                    fphase(:,ch,fno) = fphase0;
                    freal0 = real(eegconv);
                    if decimatefactor > 1
                        freal0 = decimate(freal0,decimatefactor);    %#ok<PROPLC>                        
                    end
                    ffreal(:,ch,fno) = freal0;
                end
                %fprintf('\n'); %tisk znova na stejnou radku
            end
            obj.HFreq = HFreq; %#ok<PROPLC>
            obj.fphase = fphase;%#ok<PROPLC>
            obj.freal = ffreal;
            toc(timer); %ukoncim mereni casu a vypisu,
            obj.d = squeeze(mean(obj.HFreq,3)); %11.5.2016 - prepisu puvodni data prumerem pres frekvence
            obj.fs = obj.fs/decimatefactor;
            obj.tabs = downsample(obj.tabs,decimatefactor);
            obj.tabs_orig = downsample(obj.tabs_orig,decimatefactor); %potrebuju zdecimovat i druhy tabs. Orig znamena jen ze nepodleha epochovani
            obj.Hf = freq;
            obj.Hfmean = freq; %wavelety nepocitaji silu rozmezi frekvenci,ale primo zadanych frekvenci
            obj.mults = ones(1,size(obj.d,2)); %nove pole uz je double - defaultove jednicky pro kazdy kanal
            obj.yrange = [1 1 5 5]; %zmenim rozliseni osy y v grafu
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
            fprintf('\n'); %ukoncim radku
            toc(timer); %ukoncim mereni casu a vypisu
            disp(['vytvoreno ' num2str(numel(obj.Hf)) ' frekvencnich pasem v case ' num2str(toc(timer)) 's']); 
        end
        function PasmoFrekvenceCVUT(obj,freq,channels)
            %   pouzivam VlnkovaTransformacia od Bortela            
            %   pouziva data d z parentu a take fs
            %   freq    seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi
            if ~exist('channels','var'), channels = 1:obj.channels; end
            obj.HFreq = zeros(ceil(obj.samples/obj.decimatefactor), obj.channels,numel(freq)); %inicializace pole   cas x kanal x freq
            tic; %zacnu merit cas
            %fprintf('kanal ze %i: ', numel(channels) ); 
            
            scales = 5./(2*pi*(freq)); %prevedu frekvence na scales - to jsem okopiroval z Mike X Cohena a funguje to
              % = gaussian width (or stdev) - 5=pet cyklu waveletu - pro kazdou frekvenci zvlas
            [S,f,t] = VlnkovaTransformacia(obj.d(:,channels),scales,obj.fs); %#ok<ASGLU> %rozmery S: frekvence x delka x kanaly
            S = permute(S,[ 2 3 1]); % predelam freq 1 x cas 2 x kanal 3  na cas x kanal x freq 
            fprintf('... computing mean'); %ukoncim radku
            for ch = 1:numel(channels) %S ma v sobe jen jeden kanal, ale Hfreq vsechny kanaly
                for s = 1:numel(scales)
                    eegconv = decimate(squeeze(abs(S(:,ch,s))),obj.decimatefactor); % mensi sampling rate  
                    obj.HFreq(:,channels(ch),s) = (eegconv./mean(eegconv)); %podil prumeru = prumerna hodnota
                end
            end            
            obj.d = squeeze(mean(obj.HFreq,3)); %11.5.2016 - prepisu puvodni data prumerem
            obj.fs = obj.fs/obj.decimatefactor;
            obj.tabs = downsample(obj.tabs,obj.decimatefactor);
            obj.tabs_orig = downsample(obj.tabs_orig,obj.decimatefactor); %potrebuju zdecimovat i druhy tabs. Orig znamena jen ze nepodleha epochovani
            obj.Hf = freq;
            obj.Hfmean = freq; %wavelety nepocitaji silu rozmezi frekvenci,ale primo zadanych frekvenci
            obj.mults = ones(1,size(obj.d,2)); %nove pole uz je double - defaultove jednicky pro kazdy kanal
            obj.yrange = [1 1 5 5]; %zmenim rozliseni osy y v grafu
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
            fprintf('\n'); %ukoncim radku
            toc; %ukoncim mereni casu a vypisu
            disp(['vytvoreno ' num2str(numel(obj.Hf)) ' frekvencnich pasem']);      
        end
    end
 end