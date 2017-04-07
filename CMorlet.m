classdef  CMorlet < CHilbert
    
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
        
        function obj = PasmoFrekvence(obj,freq,channels)
            %   podle figure 1311 z knihy Mike X Cohen 
            %   pouziva data d z parentu a take fs
            %   freq    seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi
            if ~exist('channels','var'), channels = 1:obj.channels; end
            obj.HFreq = zeros(ceil(obj.samples/obj.decimatefactor),obj.channels,numel(freq)); %inicializace pole   
            tic; %zacnu merit cas
            fprintf('kanal ze %i: ', numel(channels) );
            
            if freq(1) < 10 
                time = -2:1/obj.fs:2; %pro nizke frekvence pouziju delsi wavelet - +- 2s
            else
                time = -1:1/obj.fs:1; %pro vyssi frekvence staci wv +- 1s
            end
            s    =  5./(2*pi*freq);   %gaussian width (or stdev) - 5=pet cyklu waveletu - pro kazdou frekvenci zvlas
            n_wavelet            = length(time); %delka waveletu ve vzorcich - 1025 pro 512 a 1s
            n_data               = obj.samples; %data nejsou epochovana
            n_convolution        = n_wavelet+n_data-1; %pocet vzorku vysledne konvoluce - vic bodu nez maji data, takze budu pridavat 0 na konce
            n_conv_pow2          = pow2(nextpow2(n_convolution)); %prevede pocet vzorku konvoluce na nejblizsi vyssi mocninu 2 
            half_of_wavelet_size = (n_wavelet-1)/2;
            
            for ch = channels %jednotlive elektrody                                           
                fprintf('%i,',ch);
                eegfft = fft(obj.d(:,ch)',n_conv_pow2); %FFT of eeg data, potrebuju to dat do radku aby stejne jako wavelet
                for fno = 1:numel(freq) %seznam frekvenci                    
                    wavelet = fft( sqrt(1/(s(fno)*sqrt(pi))) * exp(1i*2*pi*freq(fno).*time) .* exp(-time.^2./(2*(s(fno)^2))) , n_conv_pow2 );
                    % (A=frequency band-specific scaling factor) * complex sin * gaussian
                    
                    eegconv = ifft(wavelet.*eegfft); %convoluce pomoci convolution theorem
                    eegconv = eegconv(1:n_convolution);
                    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);    
                    eegconv = abs(eegconv).^2; %z komplexniho cisla vypocitam power
                    eegconv = decimate(eegconv,obj.decimatefactor); % mensi sampling rate                    
                    obj.HFreq(:,ch,fno) = (eegconv./mean(eegconv)); %podil prumeru = prumerna hodnota
                    %fprintf('%i Hz, ',loF);
                end
                %fprintf('\n'); %tisk znova na stejnou radku
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