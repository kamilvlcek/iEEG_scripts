classdef CStat
    %CSTAT Class for custom statistical functions
    %   Kamil Vlcek, FGU AVCR, since 2016 06
    
    properties        
    end
    
    methods (Static,Access = public)
        
        function W = Wilcox2D(A,B,print,fdr,msg,RjEpChA,RjEpChB)
            %Wilcox2D(A,B,print,fdr,msg,RjEpChA,RjEpChB)
            %srovna dve 3D matice proti sobe, ohledne hodnot v poslednim rozmeru
            %A musi mit oba prvni rozmery > rozmery B, 
            %B muze mit jeden nebo oba prvni rozmer = 1 - pak se porovnava se vsemi hodnotami v A
            %pokud fdr=1 nebo prazne, provadi fdr korekci
            %pokud print = 0 nebo prazne, netiskne nic
            if ~exist('print','var'), print = 0; end
            if ~exist('fdr','var') || isempty(fdr), fdr = 1; end %min striktni je default           
            if ~exist('msg','var'), msg = ''; end
            if ~exist('RjEpChA','var'), RjEpChA = false(size(A,2),size(A,3)); end
            if ~exist('RjEpChB','var'), RjEpChB = false(size(B,2),size(B,3)); end
            W = zeros(size(A,1),size(A,2));
           
            if print, fprintf(['Wilcox Test 2D - ' msg ': ']); end
            for j = 1:size(A,1) % napr cas
                if print && mod(j,50)==0, fprintf('%d ', j); end %tisknu jen cele padesatky
                for k = 1:size(A,2) %napr kanaly                  
                   aa = squeeze (A(j,k,~RjEpChA(k,:))); %jen nevyrazene epochy 
                   bb = squeeze (B( min(j,size(B,1)) , min(k,size(B,2)) , ~RjEpChB(k,:) )); %jen nevyrazene epochy
                   if numel(aa) >= 2 && numel(bb) >= 2 
                      W(j,k) = ranksum(aa,bb); % Statistics and Machine Learning Toolbox
                   else
                      W(j,k) = 1; %pokud jen jedna hodnota, nelze delat statistika
                   end
                end
            end
            if fdr
                if fdr == 2, method = 'dep'; else method = 'pdep'; end
                [~, ~, adj_p]=fdr_bh(W,0.05,method,'no'); %dep je striktnejsi nez pdep
                W = adj_p; %prepisu puvodni hodnoty korigovanymi podle FDR
            end
            if print, fprintf('%d .. done\n',j); end
        end
        function p_corrected = PermStat(A,B,print,msg,RjEpChA,RjEpChB) %#ok<INUSD>
            %P = PermStat(A,B,print,msg,RjEpChA,RjEpChB)            
            % A a B jsou cas x channels x epochs
            % RjEpChA a B jsou channels x epochs
            % funkce od Radka Bortela pro permutacni statistiku - 29.1.2018
            if ~exist('print','var'), print = 0; end
            if ~exist('msg','var'), msg = ''; end
            if ~exist('RjEpChA','var'), RjEpChA = false(size(A,2),size(A,3)); end %#ok<NASGU>
            if ~exist('RjEpChB','var'), RjEpChB = false(size(B,2),size(B,3)); end %#ok<NASGU>                          
            %vyrazovani epoch pro jednotlive kanaly se zatim nepouziva, funkce to neumoznuje
            if (size(A,3)==0 || size(B,3)==0 ) 
                warning('Permutation Test 2D - zadne epochy');                
            end
          
            if print, fprintf(['Permutation Test 2D - ' msg ': ']); end
            timer = tic; %zadnu merit cas
            p_corrected=CmpPerm(A,B,30000,2000); %pocty permutaci, cisla doporucena od Radka
            cas = toc(timer); %ukoncim mereni casu a vypisu
            if print, fprintf(' .. done in %.1f min \n',cas/60); end            
        end
        function [W] = Klouzaveokno(W,oknosirka, funkce,dimension)
            % oknosirka je v poctu bodu, funkce muze byt min, max, mean
            % 27.4.2015 - vynato z wilcoxmap, 21.6.2016 presunuto do CStat      
                        
            if oknosirka >= 1 
                if dimension == 1 %pokud chci pocitat klouzave okno v prvnim rozmeru, transponuju na zacatku i na konci
                    W = W';
                end                
                pulsirka = ceil(oknosirka/2); %ktery sloupec se povazuje za pulku okna - tam se hodnota ulozi
                W2 = zeros(size(W,1),size(W,2)); %musim udelat kopii, jinak si prepisuju hodnoty ze kterych pak pocitam
                for sloupec = 1:size(W,2) 
                    iW = max([1 sloupec-pulsirka+1]) : min([size(W,2) sloupec-pulsirka+oknosirka]); 
                    switch funkce
                        case 'min'
                            W2(:,sloupec)=min(W(:,iW),[],2);
                        case 'max'
                            W2(:,sloupec)=max(W(:,iW),[],2);
                        case 'mean'
                            W2(:,sloupec)=mean(W(:,iW),2);
                        otherwise
                            W2(:,sloupec)=W(:,sloupec); %pokud jina funkce, nemenim matici
                    end 
                end
                if dimension == 1 
                    W = W2';
                else
                    W = W2;
                end
            end
        end
        
        function N = round(n,dig)
            %zakrouhuje na dany pocet desetinnych mist
            N = round(n * 10^dig) / 10 ^ dig;
        end
        function [frequencies,fft_d_abs]=Fourier(dd,fs,method)
            %[frequencies,fft_d_abs]=Fourier(dd,fs) 
            % spocita Fourierovu fransformaci dat z jednoho kanalu 
            %predpoklada data s casem v prvnim rozmeru, ostatni rozmery nejsou nebo 1
            % 4.5.2017 kopie z CiEEGdata.Fourier
            if ~exist('method','var'), method = 'pwelch'; end 
            switch method   
                
                case 'fft'
            %1. varianta podle MikeXCohen
            frequencies       = linspace(0,fs/2,length(dd)/2+1); %maximalni frekvence je fs/2, ale frekvencni rozliseni je N/2+1
            fft_d    = fft(dd)/length(dd); %fast fourier transform
           
            %vezmu jen tolik frekvenci, kolik je v frequencies - realnych frekvenci. ostatni jsou imaginarni frekvence
            fft_d_abs = abs(fft_d(1:length(frequencies)))*2; %dvema nasobim kvuli tem imaginarnim frekvencim. Viz MikeCohenP.            
            %prvni frekvence je DC
            fft_d_abs = 10*log10(fft_d_abs);
            
                case 'periodogram'
            %2. varianta z Matlab Answers https://uk.mathworks.com/matlabcentral/answers/114942-how-to-calculate-and-plot-power-spectral-density-of-a-given-signal
            %ale je to v podstate stejne jako fft
            [pxx,frequencies] = periodogram(dd,[],length(dd),fs);  %tohle je jen abs(fft())^2/n
            fft_d_abs = 10*log10(pxx);
            
                case 'pwelch'
            %3. pouziti pwelch, viz https://stackoverflow.com/questions/27079289/on-the-use-and-understanding-of-pwelch-in-matlab            
            M = 8*fs; %round(length(dd)/20); 
            [pxx, frequencies] = pwelch(dd,M,round(M/2),fs*20,fs); %2.-4. parametr od Radka Bortela 29.06.2017
            fft_d_abs = 10*log10(pxx);
            end
        end
            
        function [filter_result] = FIR(freq,dd,fs,firtype)
            %[filter_result] = FIR(freq,dd,fs) 
            % provede FIR filter podle Cohen Ch 14.  
            % 4.5.2017   
            if ~exist('firtype','var'), firtype = 'fir1'; end
            assert(strcmp(firtype,'fir1') || strcmp(firtype,'firls'), ['neznamy typ kernelu ' firtype]);
            
            nyquist =fs/2;                      
                        
            if strcmp(firtype,'fir1')
                %FIR1 - jednodussy kernel doporuceny v eeglab
                if numel(freq)==1 %highpass
                    filter_order       = round(3*(fs/freq(1))); %trojnasobek dolni frekvence - delka filter kernelu 
                    filterweights = fir1(filter_order,freq./nyquist, 'high'); %vytvorim filter kernel 
                elseif freq(1)==0 %lowpass
                    filter_order       = round(3*(fs/freq(2))); %trojnasobek dolni frekvence - delka filter kernelu 
                    filterweights = fir1(filter_order,freq(2)./nyquist); %vytvorim filter kernel 
                else
                    filter_order       = round(3*(fs/freq(1))); %trojnasobek dolni frekvence - delka filter kernelu 
                    filterweights = fir1(filter_order,freq./nyquist); %vytvorim filter kernel
                end
                %kdyz u fir1 dam jen jednu frekvenci, automaticky pocita lowpass - viz eeglab:eegfilt
            else
                %FIRLS
                %freqspread    = (freq(2)-freq(1))/2; % Hz +/- the center frequency
                %center_freq   = freq(1) + freqspread;   
                assert(numel(freq)>1 && freq(1)>0,'firls muze byt jen bandpass');
                filter_order       = round(3*(fs/freq(1))); %trojnasobek dolni frekvence - delka filter kernelu 
                transwid      = .10; % transition zone withth
                ffrequencies  = [ 0 (1-transwid)*(freq(1)) (freq(1)) (freq(2)) (1+transwid)*(freq(2)) nyquist ]/nyquist;
                ffrequencies( ffrequencies>1 ) = 1; %pokud jsem nastavil jako pasmo=nyquist, zlobilo by to jinak
                ffrequencies( ffrequencies<0 ) = 0; %pokud jsem nastavil jako pasmo=0, zlobilo by to jinak
                idealresponse = [ 0 0 1 1 0 0 ];
                filterweights = firls(filter_order,ffrequencies,idealresponse); %vytvorim filter kernel                          
            end
            filter_result = filtfilt(filterweights,1,dd);
            
        end
        
        function [AUC] = ROCKrivka(epochData,dd,katnames,kresli)
            %11/2018 - spocita a vykresli roc krivku
            %epochData - cellarray z CiEEGData
            %dd - eegdata pro jeden channel a jeden sample
            %sample - cislo sample
            %katnames - jmena dvou kategorii, positive a negative
            if ~exist('kresli','var'), kresli = 0; end
            labels = epochData(:,1);
            scores = dd;
            posclass = katnames{1}; %napr Scene
            if numel(katnames) > 1
                negclass = katnames{2}; %napr Face
            else
                negclass = 'all';
            end
            
            %labels = cell array/int array nazvu kategorii - co radek to event, napr Scene, Object, Face
            %scores = array se scores pro jednotlive - co radek to event - power - prumer nebo urcity cas
            %posclass = string - Positive class label - nazev vyhodnocovane kategorie
            %options - 'NegClass','versicolor' - negative class, default je all
            
            [X,Y,~,AUC] = perfcurve(labels,scores,posclass,'NegClass',negclass);            
            if kresli
                figure('Name',['ROC curve for ' posclass 'x' negclass ]);
                subplot(1,2,1);
                plot(X,Y);     
                subplot(1,2,2);
                posdata = dd(not(cellfun('isempty',strfind(epochData(:,1),posclass))));
                negdata = dd(not(cellfun('isempty',strfind(epochData(:,1),negclass))));            
                plot(rand(size(posdata))*0.5+1,posdata,'or');
                hold on;
                plot(rand(size(negdata))*0.5+2,negdata,'ob');
            end
            
        end
        function [AUC] = ROCAnalysis(E,ch,time,katnames)
            %time je cas, ve kterem chci spocitat ROC ve vterinach
            sample = round((time - E.epochtime(1))*E.fs);
            
            %musim vyradit spatne epochy
            katnums = E.PsyData.CategoryNum(katnames);
            [~,~,~,iEpP] = E.CategoryData(katnums(1),[],[],ch); %novy parametr iEp se seznamem vsech validnich epoch pro tento kanal
            [~,~,~,iEpN] = E.CategoryData(katnums(2),[],[],ch);                                                           
            
            if numel(sample) == 1 %chci udela ROC krivku jen pro jeden bod v case
                AUC = CStat.ROCKrivka(E.epochData(iEpP | iEpN,:),squeeze(E.d(sample,ch,iEpP | iEpN)),katnames,1);
            else %udalam ROC krivku pro kazd
                AUC = zeros(diff(sample)+1,1);
                for s = sample(1):sample(2)
                    AUC(s-sample(1)+1) = CStat.ROCKrivka(E.epochData(iEpP | iEpN,:),squeeze(E.d(s,ch,iEpP | iEpN)),katnames,0);                   
                end
                figure('Name', ['AUC waveform for '  katnames{1} 'x' katnames{2} ]);
                
                subplot(2,1,1); %prvni plot je AUC krivka
                X = linspace(time(1),time(2),numel(AUC));
                plot(X,AUC,'.-');
                title('AUC');
                
                subplot(2,1,2); %do eruheho plot vykreslim prubeh power obou kategorii a jejich rozdil
                title('power');
                dataP = squeeze(E.d(sample(1):sample(2),ch,iEpP));
                dataN = squeeze(E.d(sample(1):sample(2),ch,iEpN));
                hold on;
                %prvni kategorie
                MP = mean(dataP,2);
                Er= std(dataP,[],2)/sqrt(size(dataP,2)); %std err of mean                                             
                ciplot(MP+Er, MP-Er, X,[1 .8 .8]); %funguje dobre pri kopii do corelu, ulozim handle na barevny pas
                plot(X,MP,'.-r');  
                
                %druha kategorie
                MN = mean(dataN,2);
                Er = std(dataN,[],2)/sqrt(size(dataN,2)); %std err of mean
                ciplot(MN+Er, MN-Er, X,[.8 .8 1]); %funguje dobre pri kopii do corelu, ulozim handle na barevny pas
                plot(X,MN,'.-b');                              
                                
                plot(X,MP-MN,'o-g'); %zelene bude rozdil
                
                
            end
        end
    end
    
end

