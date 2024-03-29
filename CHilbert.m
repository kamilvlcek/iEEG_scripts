classdef CHilbert < CiEEGData
    %HILBERT.CLASS sbirka funkci na analyzu pomoci hilbert trasform
    %   Kamil Vlcek, FGU AVCR, since 2016 04
    properties (Constant = true)
        decimatefactor = 8; %o kolik decimuju hiblertovu obalku oproti puvodi sampling rate; 2 je dostatecne konzervativni, 8 hodne setri pamet
    end
    properties (Access = public)
        HFreq; %hilberova obalka pro kazde frekvenci pasmo - time x channel x freq (x kategorie)
        HFreq_ChEpochs; %true for channels/epochs, which are included in Hfreq, channels x epochs
        HFreqEpochs; %Hf bez priemerovania cez epochy - time x channel x frequency x epoch        
        fphaseEpochs; % epochovane fazy fphase ~HFreqEpochs
        Hf; %frekvencni pasma pro ktere jsou pocitany obalky - okraje pasem, pocet je tedy vetsi o 1 nez pocet spocitanych pasem
        Hfmean; %stredni hodnoty pasem  - pocet = pocet spocitanych pasem
        hfilename; %jmeno souboru CHilbert  
        plotF = struct; %udaje o stavu plotu PlotResponseFreq
        fphase; %faze vsech zpracovavanych frekvenci - premiestnene z CMorlet pre vykreslenie a porovnanie faz z MW a Hilberta do buducna    %time x channel x frequency    
        frealEpochs; % epochovane filtrovane eeg
        normalization; %typ normalizace
    end
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS 
        function obj = CHilbert(d,tabs,fs,mults,header)            
            if ~exist('header','var'), header = []; end %nejakou hodnotu dat musim
            if ~exist('mults','var'),  mults = []; end %nejakou hodnotu dat musim
            if ~exist('d','var') || isempty(d) %konstruktor uplne bez parametru - kvuli CHilbertMulti
                d = []; tabs = []; fs = [];
            elseif ischar(d) && ~exist('fs','var') %pokud je prvni parametr retezec, tak ho beru jako nazev souboru, ktery nactu
                fs = []; 
                if ~exist('tabs','var'), tabs=[]; end
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
            HFreq = zeros(samples,obj.channels,numel(freq)-1); %#ok<PROPLC> %inicializace pole   
            timer = tic; %zadnu merit cas
            fprintf('kanal ze %i: ', max(channels) );
            freq1 = freq(2:end); %kvuli parfor, aby bylo konzistentni indexovani
            fs = obj.fs; %jen kvuli parfor  
            fnocyclenum = numel(freq)-1; %kvuli parfor, aby jasny pocet cyklu
            for ch = channels %jednotlive elektrody 
                % outer parfor mi zatim nefunguje, hrozne pameti se nacetlo (d se asi ze ctyrnasobilo pro 4 workers) a pak to vyhodilo chybu:
                  
                %fprintf('channel %i: Hz ',ch);                         
                fprintf('%i,',ch); 
                d = obj.d(:,ch);
                if sum(d)==0, continue; end %pro vyrazene kanaly jsou hodnoty 0 pri jine nez bipol ref. Z tech pak vznikne nan, pri tomhle cyklu, coz vadi dal                                  
                parfor fno = 1:fnocyclenum %seznam frekvenci
                    loF = freq(fno) -prekryv*(freq1(fno)-freq(fno)); 
                    hiF = freq1(fno)-0.1 +prekryv*(freq1(fno)-freq(fno));  %napr 50 - 59.9
                    hh = CHilbert.hilbertJirka(d,loF,hiF,fs); %cista hilbertova obalka, tohle i skript hodne zrychli
                    hh = decimate(hh,decimatefactor); % mensi sampling rate                    
                    HFreq(:,ch,fno) = hh;  %#ok<PROPLC> povodna normalizacia (hh./mean(hh)) premiestnena do funkcie Normalize
                    %fprintf('%i Hz, ',loF);
                end
                %fprintf('\n'); %tisk znova na stejnou radku
            end
            obj.HFreq = HFreq; %#ok<PROPLC>
            toc(timer); %ukoncim mereni casu a vypisu, skore inner parfor=136s, for=175s
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
            obj.DatumCas.HilbertComputed = datestr(now);
            disp(['vytvoreno ' num2str(numel(obj.Hfmean)) ' frekvencnich pasem v case ' num2str(toc(timer)) 's']); 
        end
        function obj = Normalize(obj, type, channels)
         %function for different normalization methods
         %to be used after PasmoFrekvence
         %TODO musi se osetrit, kdyz je prumer 0
         %according to MikeXCohen the normalization should be done using mean and std of baseline. 
         
            if ~exist('channels','var'), channels = 1:obj.channels; end
            if ~exist('type','var') || isempty(type), type='orig'; end %default normalization 
            fprintf('Normalize %s: ', type);
            switch type
                case 'orig' %(fpower./mean(fpower)) - povodna normalizacia
                    fprintf('channel:');
                    for ch = channels
                        fprintf('%i,', ch );
                        for f = 1:size(obj.HFreq,3) %over frequency bands
                            obj.HFreq(:,ch,f) = obj.HFreq(:,ch,f)./mean(obj.HFreq(:,ch,f)); 
                        end                        
                    end
                case 'mean' %(fpower./mean(fpower))*100 ~ Bastin 2012   
                    for ch = channels
                        for f = 1:size(obj.HFreq,3)
                            obj.HFreq(:,ch,f) = (obj.HFreq(:,ch,f)/mean(obj.HFreq(:,ch,f)))*100;
                        end
                    end
                case 'z' %z-transform (fpower-mean(fwpower))/std(fpower)
                    for ch = channels
                        for f = 1:size(obj.HFreq,3)
                            obj.HFreq(:,ch,f) = (obj.HFreq(:,ch,f)-mean(obj.HFreq(:,ch,f)))./std(obj.HFreq(:,ch,f));
                        end
                    end
                case 'db' %db=10*log10(fpower/mean(fpower)) 
                    for ch = channels
                        for f = 1:size(obj.HFreq,3)
                            chdata = obj.HFreq(:,ch,f)./mean(obj.HFreq(:,ch,f));
                            obj.HFreq(:,ch,f) = 10*log10(chdata + 1 - min(chdata) ); % values >=1 for log10
                        end
                    end
                case 'valmax' %normalize all data to have the same average valmax of the response
                    assert(size(obj.d,3) > 1,'the data should be epoched');
                    kategories = flip(obj.Wp(obj.WpActive).kats);    %categories according to currently selected contrast                 
                    fprintf('category:');
                    for kk = 1:numel(kategories)                        
                       katnum = cellval(kategories,kk); % if kategories is cellarray, the katnum will also be cellarray %cellval(kategories,k); %funkce cellval funguje at je to cell array nebo neni
                       valmax = obj.ResponseTriggerTime(1, 1, {katnum}, channels); %valmax for all channels %here katnum should be in one column                                     
                       for k = 1:numel(katnum) %categories individually if they are combined
                            fprintf('%i,', cellval(katnum,k));
                            %for f = 1:size(obj.HFreq,3) %the frequency bands
                            %dimensions of HFreq are: samples, channels, frequency bands, kategories
                            obj.HFreq(:,:,:,cellval(katnum,k)+1) = obj.HFreq(:,:,:,cellval(katnum,k)+1)./valmax; %matlab matches 2. dimentions (channels) with valmax dimension by its size
                            %end                   
                            iEpochy =  ismember(cell2mat(obj.epochData(:,2)), cellval(katnum,k)) ;   %epochs with this stimulus category                           
                            %dimensions of d are: samples, channels, epochs
                            obj.d(:,:,iEpochy) = obj.d(:,:,iEpochy)./valmax; %M2016b and after: implicit expansion of arrays with compatible sizes.                            
                       end                       
                    end                                        
                otherwise
                    error(['unknown normalization: ' type]);
            end
            if size(obj.d,3) == 1 %only for non-epoched data
                obj.d = squeeze(mean(obj.HFreq,3)); %mean over freq bands 
            end
            obj.normalization = type; %pro zpetnou referenci
            fprintf('... finished\n');
        end
        function obj = ExtractEpochs(obj, PsyData,epochtime,baseline,freqepochs,filter)
            % rozdeli hilbertovu obalku podle epoch
            % i u ni odecte baseline pred podnetem
            % freqepochs - urcuje jestli ukladat freq pasma pro vsechny epochy do pole obj.HFreqEpochs
            % filter - cellarray: filter{1} - column number in obj.P.data, filter{2} - looked for values in this column
            if ~exist('baseline','var') || isempty(baseline), baseline = [epochtime(1) 0]; end %defaultni baseline je do 0 sec
            if ~exist('freqepochs','var') || isempty(freqepochs), freqepochs = 0; end %defaultne se neukladaji frekvencni pasma pro vsechny epochy
            if ~exist('filter','var') || isempty(filter), filter = []; end %default filter is empty  
            
            ExtractEpochs@CiEEGData(obj,PsyData, epochtime,baseline,filter); %to mi zepochuje prumernou obalku za frekvencni pasma v poli d
            fprintf('CHilbert.ExtractEpochs: category ' );
            if(numel(obj.HFreq)>0)
                %now epoching of all frequency band individually, mostly because of figures
                %average for each category, no statistic needed from this 
                 iepochtime = round(obj.epochtime(1:2).*obj.fs); %v poctu vzorku cas pred a po udalosti, prvni cislo je zaporne druhe kladne             
                 ibaseline =  round(obj.baseline.*obj.fs); %v poctu vzorku cas pred a po udalosti
                 kategorie = cell2mat(obj.PsyData.P.strings.podminka(:,2)); %cisla kayegorii ve sloupcich
                 iepochs = obj.PsyData.FilteredIn((1:obj.epochs)',filter); %to be processed epochs according to the filter,  array of 0/1 for each epoch
                 iepochs = iepochs & ~any(obj.PsyData.GetErrorTrials(),2); %exclude all epochs with errors
                 iepochs(obj.RjEpoch) = 0; %reject globally rejected epochs
                 Hfreq2 = zeros(iepochtime(2)-iepochtime(1), size(obj.d,2), numel(obj.Hfmean),size(kategorie,1)); %new epoched power data: time x channel x freq x kategorie=podminka
                 HFreq_ChEpochs = false(size(obj.RjEpochCh));  %#ok<PROPLC> %true for channels/epochs, which will be included in Hfreq, channels x epochs
                 if freqepochs %if we want to store all epochs for all frequency bands - for statistics, much larger data
                     %nan values - the data for rejected epochs will stay nan
                     obj.HFreqEpochs = nan(iepochtime(2)-iepochtime(1),size(obj.HFreq,2),size(obj.HFreq,3),sum(iepochs)); % time x channel x frequency x epoch
                     obj.fphaseEpochs = nan(iepochtime(2)-iepochtime(1),size(obj.HFreq,2),size(obj.HFreq,3),sum(iepochs)); % time x channel x frequency x epoch
                     obj.frealEpochs = nan(iepochtime(2)-iepochtime(1),size(obj.HFreq,2),size(obj.HFreq,3),sum(iepochs)); % time x channel x frequency x epoch
                 else %normal analysis                    
                     obj.HFreqEpochs = [];
                     obj.fphaseEpochs = [];
                     obj.frealEpochs = [];                  
                 end
                 %cycle over categories, not over epochs 
                 for katnum = kategorie' %categories in rows 
                     iepochskat = iepochs & (cell2mat(obj.epochData(:,2))==katnum); 
                     Epochs = find(iepochskat);  %epochs x 1: epochs numbers for this category
                     if (numel(Epochs)) > 0 % when using filter, for some categorie, we will have no epochs. 
                         fprintf('%i,', katnum);                      
                         nEpochsChKat = zeros(1,obj.channels); %number of epoch data stored in Hfreq2 for each channel for this category
                         for epoch = Epochs' %number of epochs in a column - convert into row
                             izacatek = find(obj.tabs_orig==obj.epochData{epoch,3}); %find sample number of stimulus for this epochs, using its timestamp. In the third columnt of epochData stimulus timestamps
                             for ch=1:obj.channels %over channels for this epoch
                                if obj.RjEpochCh(ch,epoch)==0 % 23.10.2023 - if this epoch for this channel is not rejected
                                    if ibaseline(1)==ibaseline(2) %baseline not used if there is the same start and end time
                                        baseline_mean = 0;
                                        epoch_data = obj.HFreq(izacatek + iepochtime(1) : izacatek+iepochtime(2)-1,ch,:); %timex1xfreq - epoch data with baseline substracted for the current epoch and channel
                                    else 
                                        baseline_mean = mean(obj.HFreq(izacatek + ibaseline(1) : izacatek+ibaseline(2)-1,ch,:),1); %1x1xfreq: baseline for each freq band  - mean over time                            
                                        %epoch_data = bsxfun(@minus,obj.HFreq(izacatek + iepochtime(1) : izacatek+iepochtime(2)-1,ch,:) , baseline_mean); %timex1xfreq - epoch data with baseline substracted for the current epoch and channel
                                        epoch_data = obj.HFreq(izacatek + iepochtime(1) : izacatek+iepochtime(2)-1,ch,:) - baseline_mean; %timex1xfreq - epoch data with baseline substracted for the current epoch and channel
                                            %requires implicit expansion of arrays with compatible size in Matlab 2016b and later, 
                                    end
                        
                                    Hfreq2(:,ch,:,katnum+1) = Hfreq2(:,ch,:,katnum+1) + epoch_data; %sum over all epochs of power for this channel and category, over all timesamples and frequecies
                                        %for this channel and katnum, this line is executed ones for each epoch
                                    nEpochsChKat(1,ch) = nEpochsChKat(1,ch)  + 1;
                                    HFreq_ChEpochs(ch,epoch)=1;  %#ok<PROPLC>
                                    if freqepochs  %if we want to store all epochs for all frequency bands - for statistics, much larger data
                                        obj.HFreqEpochs(:,ch,:,epoch) = epoch_data - baseline_mean; % timeXfreq - here data for each epoch, channel and freq separately are saved, independent of category
                                            %normalize power for all epochs
                                        if isprop(obj,'fphase') && ~isempty(obj.fphase) %if phase data were created in CMorlet.PasmoFrekvence
                                            obj.fphaseEpochs(:,ch,:,epoch) = squeeze(obj.fphase(izacatek + iepochtime(1) : izacatek + iepochtime(2)-1, ch, :)); %phase data for this epoch and channel (no baseline substracted
                                        end                                
                                        if isprop(obj,'freal') && ~isempty(obj.freal) %if this was created - the bandpass-filtered signal - projection on the real axis                                    
                                            obj.frealEpochs(:,ch,:,epoch) = squeeze(obj.freal(izacatek + iepochtime(1) : izacatek + iepochtime(2)-1, ch, :));
                                        end
                                    end    
                                end                                                        
                             end
                         end                     
                         Hfreq2(:,:,:,katnum+1) = Hfreq2(:,:,:,katnum+1)./nEpochsChKat(1,:); %average over all epochs in this category  - the stored sum divided by number 
                     end
                 end             
                 obj.HFreq = Hfreq2;
                 obj.HFreq_ChEpochs = HFreq_ChEpochs;  %#ok<PROPLC>
                 fprintf('\n');
            end            
        end
        function obj = NormalizeEpochs(obj, baseline)
            % normalizes the obj.HFreq property by substracting mean baseline activity in each channel and category for individual frequency bands
            % normalizes also HFreqEpochs if it exists 
            % needed after appending two CHilbert objects (joining two parts of epochs, e.g. in memact test in delayed epochs with jitter)  
            % Sofiia 2023/07 Memact 
            assert(obj.epochs > 1, 'data should be epoched');
            if ~exist('baseline','var') || isempty(baseline), baseline = [obj.epochtime(1) 0]; end % by default, baseline time - the whole period before stimulus
            
            NormalizeEpochs@CiEEGData(obj, baseline);
            
            iepochtime = round(obj.epochtime(1:2).*obj.fs); % epochtime in samples
            ibaseline =  round(baseline.*obj.fs); % baseline in samples 
            baseline_mean = mean(obj.HFreq(-iepochtime(1)+ibaseline(1)+1 : -iepochtime(1)+ibaseline(2), :, :, :),1); % 1 x channels x frequencies x categories, average baseline for each channel, freq and category
            HFreqnorm = obj.HFreq - baseline_mean; % time x ch x frequencies x categories - substract mean baseline from the entire epoch
                 %requires implicit expansion of arrays with compatible size in Matlab 2016b and later
            obj.HFreq = HFreqnorm; % replace by normalized data
            
            if ~isempty(obj.HFreqEpochs) % data for all epochs for all frequency bands
                baseline_meanEp = mean(obj.HFreqEpochs(-iepochtime(1)+ibaseline(1)+1 : -iepochtime(1)+ibaseline(2), :, :, :),1); % 1 x channel x frequency x epoch, average baseline
                HFreqEpochsnorm = obj.HFreqEpochs - baseline_meanEp; % time x channels x frequencies x epochs, substract mean baseline from each individual epoch
                obj.HFreqEpochs = HFreqEpochsnorm; % replace by normalized data
            end
        end
        function GetITPC(obj)
            assert(obj.epochs > 1, 'data musi byt epochovana');
            PN = CPlotsN(obj);
            fprintf('CalculateITPC ... ');
            [itpc, ~, itpc_pmean,n_epoch] = PN.CalculateITPC(obj.Wp(obj.WpActive).kats); %TODO nada pridat parametr kats, ktery umozni i cell array
            % d = time x channel x epoch
            % itpc_pmean = ch x condition x time - p values z mean itpc pres frekvence
            
            fprintf('CalculateITPCdiffs ... ');            
            [ditpc, ditpc_p, ~,diffs] = PN.CalculateITPCcontrasts(obj.Wp(obj.WpActive).kats,itpc,n_epoch);
            % ditpc_pmean = ch x diffs x time - p values z mean ditpc pres frekvence
            % ditpc_p - ch x contrast1 x contrast2 x time x fq - p values
            % citpc_pmean ch x contrast1 x contrast2 x time - p values z mean itpc pres frekvence
            
            itpc = permute(itpc,[3 1 4 2]); %time x channel x freq x epoch
            obj.HFreq = movmean(itpc,5,1); 
            obj.d = movmean(squeeze(mean(itpc,3)),5,1); %prumer pres frekvence            
            obj.PsyData.Cond2Epochs();            
            obj.RejectEpochs([],[]);
            obj.RjEpoch = [];
            obj.epochs = size(itpc,4); 
            katnum = obj.PsyData.Categories();
            epochData2 = cell(numel(katnum),3);
            for k = 1:numel(katnum)
               ikat = find([obj.epochData{:,2}]==katnum(k),1);
               epochData2(k,:) = obj.epochData(ikat,:);
            end
            obj.epochData = epochData2;
            %smazu statistiku, aby nebyla na obrazcich
            EEEGStat = CEEGStat(obj.d,obj.fs);
            baseline = EEEGStat.Baseline(obj.epochtime,obj.baseline);
            ibaseline = round(baseline.*obj.fs);       
            iepochtime = round(obj.epochtime(1:2).*obj.fs);     
            
            for k = 1:numel(obj.Wp(obj.WpActive).kats)
                stat = itpc_pmean(:,obj.Wp(obj.WpActive).kats(k)+1,abs(iepochtime(1)-ibaseline(2))+1 : end ); %vsechny hodnoty po konci baseline  
                obj.Wp(obj.WpActive).WpKatBaseline{k,1} = permute(squeeze(stat),[2 1]); %chci mit rozmer time x ch
                %TODO pridat FDR korekci
                %pridat korekci na podminky, nebo jinou kterou navrhuje Cohen? Korekce pro cas neni potreba, pocita se to podle hladiny, ktera je zavisla na poctu epoch. 
                
                %obj.Wp(obj.WpActive).WpKatBaseline{w} = ones(size(obj.Wp(obj.WpActive).WpKatBaseline{w}));
            end            
            
            for d = 1:size(diffs,1)
                stat = squeeze(ditpc_p(:,diffs(d,1),diffs(d,2),abs(iepochtime(1)-ibaseline(2))+1 : end,:));
                [~,iminp] = min(mean(stat,2),[],3);  %mean pres cas, min pres frekvence              
                stat2 = zeros(size(stat,1),size(stat,2));
                for ch = 1:size(stat,1)
                    stat2(ch,:) = stat(ch,:,iminp(ch));
                end
                obj.Wp(obj.WpActive).WpKat{diffs(d,1),diffs(d,2)} = permute(stat2,[2 1]);  %chci mit rozmer time x ch 
                %obj.Wp(obj.WpActive).WpKat{w} = ones(size(obj.Wp(obj.WpActive).WpKat{w}));    
                %TODO - kontrola - vracet frekvence s maximalni signif a do d ukladat ne prumer ale tu frekvenci
            end
            fprintf('done\n');
        end
        function obj = Decimate(obj,podil,rtrim)
            %zmensi frekvencni data na nizsi vzorkovaci frekvenci
            if ~exist('rtrim','var') || isempty(rtrim), rtrim = []; end 
            Decimate@CiEEGData(obj,podil,rtrim);
            if obj.decimatefactor == 1 %zatim pouzivam pouze, pokud nejsou data uz decimovana vuci CiEEGdata
                fprintf('channels to decimate HFreq (z %i):',numel(obj.channels));
                HFreq = zeros(ceil(size(obj.HFreq,1)/podil) , size(obj.HFreq,2), size(obj.HFreq,3),size(obj.HFreq,4));  %#ok<PROPLC>
                if ~isempty(obj.HFreqEpochs)
                    HFreqEpochs = zeros(ceil(size(obj.HFreq,1)/podil) , size(obj.HFreq,2), size(obj.HFreq,3),obj.epochs); %#ok<PROPLC>
                else
                    HFreqEpochs = []; %#ok<PROPLC>
                end
                for ch = 1:obj.channels                    
                    fprintf('%i, ',ch);
                    for f = 1:size(obj.HFreq,3) %pocet frekvencnich pasem
                       for kat = 1:size(obj.HFreq,4) %pocet kategorii podnetu     
                            HFreq(:,ch,f,kat) = decimate(obj.HFreq(:,ch,f,kat),podil); %#ok<PROPLC>                            
                       end
                       if ~isempty(obj.HFreqEpochs)
                       for ep = 1:obj.epochs %frekvencni data se vsemi epochami
                            HFreqEpochs(:,ch,f,ep) = decimate(obj.HFreqEpochs(:,ch,f,ep),podil); %#ok<PROPLC>    
                       end
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
            %plot of all frequency bands
            %which channel
            if ~exist('ch','var')
                if isfield(obj.plotF,'ch'), ch =  obj.CH.sortorder(obj.plotF.ch); %vytahnu cislo kanalu podle ulozeneho indexu
                else, obj.plotF.ch = 1; ch =  obj.CH.sortorder(1); end
            else
                obj.plotF.ch = ch; %tady bude ulozeny index sortorder, parametr ch urcuje index v sortorder
                ch =  obj.CH.sortorder(ch); %promenna ch uz urcuje skutecne cislo kanalu
            end
            %stimulus/response categories
            if ~exist('kategories','var')
                if isfield(obj.plotF,'kategories') 
                    kategories = obj.plotF.kategories;
                elseif isfield(obj.Wp(obj.WpActive), 'kats')
                    kategories = obj.Wp(obj.WpActive).kats;                 
                else
                    kategories = obj.PsyData.Categories(); 
                    obj.plotF.kategories = kategories;
                end
            else
                obj.plotF.kategories = kategories;
            end
            %how to plot
            if isfield(obj.plotF,'fh') && ishandle(obj.plotF.fh)
                figure(obj.plotF.fh); %pouziju uz vytvoreny graf
                %clf(obj.plotF.fh); %graf vycistim
            else
                obj.plotF.fh = figure('Name','ResponseFreq','Position', [20, 500, 1200, 300]);
                colormap jet; %aby to bylo jasne u vsech verzi matlabu - i 2016            
            end   
           
            T = linspace(obj.epochtime(1),obj.epochtime(2),size(obj.d,1)); %time scale - all samples
            
            maxy = 0; %scale of the power (z or y axis), over all categories
            miny = 0;            
            for k = 1:numel(kategories) %cycle over categories
                subplot(1,numel(kategories),k);                                 
                if iscell(kategories(k))                    
                    dd = zeros(size(obj.HFreq,1),size(obj.HFreq,3),numel(kategories{k}));
                    for ikat = 1:numel(kategories{k})
                        dd(:,:,ikat) = squeeze(obj.HFreq(:,ch,:,kategories{k}(ikat)+1));                        
                    end
                    D = mean(dd,3); %mean over combined categories
                else
                    D = squeeze(obj.HFreq(:,ch,:,kategories(k)+1));
                end
                imagesc(T,obj.Hfmean, D');%middles of the frequency bands
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
                caxis([miny,maxy]);     %colormap limits          
                title( obj.PsyData.CategoryName(cellval(kategories,k)), 'Interpreter', 'none');
                if k == 1
                    chstr = iff(isempty(obj.CH.sortedby),num2str(ch), [ num2str(ch) '(' obj.CH.sortedby  num2str(obj.plotF.ch) ')' ]);
                    if ~isempty(obj.HFreq_ChEpochs) && size(obj.HFreq_ChEpochs,1)>=ch                        
                        chstr = [ chstr ' (' num2str(sum(obj.HFreq_ChEpochs(ch,:))) ' epochs) ']; %#ok<AGROW>
                    end
                    ylabel(['channel ' chstr ' - freq [Hz]']); 
                    if isprop(obj,'plotRCh') && isfield(obj.plotRCh,'selCh') && any(obj.plotRCh.selCh(ch,:),2)==1        
                        klavesy = 'fghjkl'; %abych mohl vypsat primo nazvy klaves vedle hvezdicky podle selCh
                        text(0,obj.Hf(1),['*' klavesy(logical(obj.plotRCh.selCh(ch,:)))], 'FontSize', 15,'Color','red');
                    end
                end
                
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
                HFreq = obj.HFreq;   %#ok<PROPLC,NASGU>
                HFreq_ChEpochs = obj.HFreq_ChEpochs; %#ok<PROPLC,NASGU> %ch,epoch - bool - from which epochs HFreq was computed
                Hf = obj.Hf;         %#ok<PROPLC,NASGU> 
                Hfmean = obj.Hfmean; %#ok<PROPLC,NASGU> 
                HFreqEpochs = obj.HFreqEpochs; %#ok<PROPLC,NASGU> %time x channel x frequency x epoch                
                yrange = obj.yrange; %#ok<NASGU>
                fphase = obj.fphase; %#ok<NASGU,PROPLC> %15.5.2018
                fphaseEpochs = obj.fphaseEpochs; %#ok<NASGU,PROPLC> %15.5.2018
                frealEpochs = obj.frealEpochs;  %#ok<NASGU,PROPLC> %30.5.2018
                if isprop(obj,'freal'), freal = obj.freal; else, freal=[]; end   %#ok<NASGU>                 
                normalization = obj.normalization; %#ok<NASGU,PROPLC> 
                save(CHilbert.filenameH(filename),'HFreq','Hf','Hfmean','HFreqEpochs','freal','frealEpochs','yrange','fphase','fphaseEpochs','normalization','HFreq_ChEpochs','-v7.3'); %do druheho souboru data z teto tridy
            end
        end
        %pokud je treti parametr 1, nenacitaji se data z nadrazene tridy
        function obj = Load(obj,filename,loadall,onlyself)
            %parametr loadall se hodi pro FE data se vsemi ulozenymi epochami, ktere jsou giganticke
            if ~exist('loadall','var') || isempty(loadall) , loadall = 1; end
            if ~exist('onlyself','var') || onlyself == 0
                assert(exist(CHilbert.filenameE(filename),'file')==2, ['soubor s daty CHilbert neexistuje:' newline CHilbert.filenameE(filename) newline 'mozna se jedna o data tridy CiEEGData?']);    
                Load@CiEEGData(obj,CHilbert.filenameE(filename));  
            end
            if exist(CHilbert.filenameH(filename),'file')  
                filenameH = CHilbert.filenameH(filename);
                load(filenameH,'HFreq','Hf','yrange');
                obj.HFreq = HFreq;        %#ok<CPROPLC>                
                obj.Hf = Hf;               %#ok<CPROPLC>                 
                obj.yrange = yrange;
                vars = whos('-file',filenameH);
                if ismember('Hfmean', {vars.name}) %7.4.2017
                    load(filenameH,'Hfmean');      obj.Hfmean = Hfmean; %#ok<CPROPLC>
                else
                    obj.Hfmean = (obj.Hf(1:end-1) + obj.Hf(2:end)) ./ 2;
                end
                
                if ismember('HFreq_ChEpochs', {vars.name}) %7.4.2017
                    load(filenameH,'HFreq_ChEpochs');      obj.HFreq_ChEpochs = HFreq_ChEpochs; %#ok<CPROPLC>
                else
                    obj.HFreq_ChEpochs = [];
                end
                
                if ismember('fphase', {vars.name}) %15.5.2018
                    load(filenameH,'fphase');      obj.fphase = fphase; %#ok<CPROPLC>
                else
                    obj.fphase = [];
                end
                if ismember('normalization', {vars.name}) 
                    load(filenameH,'normalization');      obj.normalization = normalization; %#ok<CPROPLC>
                else
                    obj.normalization = [];
                end
                if loadall 
                    if ismember('HFreqEpochs', {vars.name}) %7.4.2017
                        load(filenameH,'HFreqEpochs');      obj.HFreqEpochs = HFreqEpochs; %#ok<CPROPLC>
                    else
                        obj.HFreqEpochs = [];
                    end
                    if ismember('fphaseEpochs', {vars.name}) %15.5.2018
                        load(filenameH,'fphaseEpochs');      obj.fphaseEpochs = fphaseEpochs; %#ok<CPROPLC>
                    else
                        obj.fphaseEpochs = [];
                    end
                    if ismember('frealEpochs', {vars.name}) %15.5.2018
                        load(filenameH,'frealEpochs');      obj.frealEpochs = frealEpochs; %#ok<CPROPLC>
                    else
                        obj.frealEpochs = [];
                    end
                else
                    obj.HFreqEpochs = [];
                    obj.fphaseEpochs = [];
                    obj.frealEpochs = [];
                end
                
                disp(['nacten soubor ' CHilbert.filenameH(filename)]); 
            else
                warning(['soubor s frekvencnimi pasmy neexistuje ' CHilbert.filenameH(filename)]);
            end
            obj.hfilename = filename; 
        end 
        function [filename,basefilename,skipped] = ExtractData(obj,PAC,label,overwrite)
            %vytvori data z vyberu elektrod, pro sdruzeni elektrod pres vsechny pacienty. 
            %pole d, tabs, RjEpochCh a header H
            %jen epochovana data, bipolarni reference
            %overwrite urcuje, jestli se maji prepisovat existujici soubory
            if ~exist('overwrite','var'), overwrite = 0; end %defaultne se nemaji prepisovat existujici soubory
            assert(obj.epochs > 1,'data have to be epoched');
            chns = [PAC.ch];
            [filepath,fname,ext] = CHilbert.matextension(obj.filename);
            podtrzitko = strfind(fname,'_'); %chci zrusit cast za poslednim podtrzitkem
            basefilename = [fname(1:podtrzitko(end)-1) ' ' label '_Extract' ext]; %jmeno bez cesty
            extracts_path = [filepath filesep 'Extracts']; %podadresar pro extrakty dat
            if exist(extracts_path,'dir') ~= 7, mkdir(extracts_path);   end
            filename =[extracts_path filesep basefilename]; %podadresar pro extrakty
            
            if exist(filename','file')~=2 || overwrite 
                %pokraduju jen pokud extrakt neexistuje nebo se ma prepsat
                %assert(strcmp(obj.reference,'Bipolar'),'neni bipolarni reference');
                d = obj.d(:,chns,:); %#ok<NASGU> %vsechny casy a epochy, vyber kanalu
                tabs = obj.tabs; %#ok<NASGU> %to je spolecne pro vsechny kanaly; time x epochs
                tabs_orig = obj.tabs_orig; %#ok<NASGU> 
                fs = obj.fs; %#ok<NASGU> 
                P = obj.PsyData.P;  %psychopy data
                if isempty(P.pacientid) || numel(P.pacientid)<=1 % u nekterych pacientu je pacientid praznde
                    P.pacientid = obj.CH.PacientTag();
                end
                epochtime = obj.epochtime; %#ok<NASGU> %abych vedel kde je podnet
                baseline = obj.baseline; %#ok<NASGU> 
                RjEpochCh = obj.RjEpochCh(chns,:); %#ok<NASGU> %kanaly vs epochy
                RjEpoch = obj.RjEpoch; %#ok<NASGU> %na to jsem zapomnel - 17.8.2018
                epochData = obj.epochData; %#ok<NASGU> %identita jednotlivych epoch. Musi byt stejna pres pacienty
                DatumCas = obj.DatumCas;
                DatumCas.Extracted = datestr(now);          
                H = obj.CH.H;
                if isfield(H,'electrodes') %Sofiia 23.6.2023
                    H = rmfield(H,'electrodes'); %smazu nepotrebna pole
                end
                H = rmfield(H,'selCh_H');
                H = rmfield(H,'triggerCH');
                H.channels = H.channels(chns); %vyfiltruju kanaly
                subjName = obj.PacientID(false);
                for ch = 1:numel(H.channels)
                    H.channels(ch).name = [subjName ' ' H.channels(ch).name]; %pridam ke jmenu kanalu jmeno subjektu
                end
                Hf = obj.Hf; %#ok<PROPLC> 
                Hfmean = obj.Hfmean;  %#ok<PROPLC> 
                if isempty(Hfmean), Hfmean = (Hf(1:end-1) + Hf(2:end)) ./ 2; end %#ok<PROPLC,NASGU>             
                HFreq = obj.HFreq(:,chns,:,:); %#ok<PROPLC,NASGU>  %time x channel x freq (x kategorie)             
                Wp = obj.Wp;  %#ok<NASGU>  %exportuju statistiku
                reference = obj.reference; %#ok<NASGU>  %exportuju referenci
                if isfield(PAC,'label')
                    chnlabels =  {PAC(:).label}; %#ok<NASGU> %optional field generated in CBrainPlot.MNIFind 
                else
                    chnlabels =  cell(size(PAC)); %#ok<NASGU>
                end
                save(filename,'d','tabs','tabs_orig','fs','P','epochtime','baseline','RjEpochCh','RjEpoch','epochData','DatumCas','H','Hf','Hfmean','HFreq','Wp','reference','chnlabels','label','-v7.3');                 
                %for now its impossible to import fphase and freal to CHilberMulti, as they are unepoched and of variable lenght
                %so there no reason to save them here                
                if isprop(obj,'HFreqEpochs') && ~isempty(obj.HFreqEpochs) %for epoched data with all epochs for all freq saved 
                    HFreqEpochs = obj.HFreqEpochs(:,chns,:,:); %#ok<NASGU,PROPLC> % time x channel x frequency x epoch
                    fphaseEpochs = obj.fphaseEpochs(:,chns,:,:); %#ok<NASGU,PROPLC> %time x channel x freq x epoch
                    frealEpochs = obj.frealEpochs(:,chns,:,:); %#ok<NASGU,PROPLC>
                    save(filename,'HFreqEpochs','fphaseEpochs','frealEpochs','-append'); %pridam k existujicimu souboru 
                end
                disp(['extract saved to "' basefilename '"']);
                skipped = 0;
            else
                disp(['extract already exists, skipped: "' basefilename '"']);
                skipped = 1;
            end
                
        end
        function BPD = ExtractBrainPlotData(obj,kategorie,signum,dofig)
            %vytvori Brain Plot Data, pro CBrainPlot.PlotBrain3D
            %kategorie, hodnoty a jejich signifikance vytahne na zaklade aktualne zvolene statistiky
            %kategorie jsou indexy v katsnames, ktere vzniknou spojenim kategorii vuc baseline, kategorii vuci sobe a vsech elektrod. Pro tri kategorie je to 1-7
            
            BPD = struct;                                   
            if ~exist('signum','var') || isempty(dofig) , signum = 0; end %jestli chci jen kat1>kat2 (1), nebo obracene (-1), nebo vsechny (0)
            if ~exist('dofig','var'), dofig = 0; end %jestli chci obrazek z IntervalyResp
            BPD.signum = signum; %jen abych mel info v datech
            if isprop(obj,'label') && ~isempty(obj.label)
                [~,intervaly,~] = CHilbertMulti.GetLabelInfo(obj.label); %zjisti interval z label, pokud label existuje
                if ~isnumeric(intervaly) %muze to byt retezec s nazvem oblasti
                   intervaly = [0.1 obj.epochtime(2)]; 
                end
            else
                intervaly = [0.1 obj.epochtime(2)]; 
            end  
            
            [prumery, MNI,names,intervaly,katsnames,neurologyLabels] = obj.IntervalyResp(intervaly,[],signum,dofig); %tohle funguje pro aktualni nastavenou statistiku
            %prumery obsahuji 0 tam, kde neni signif rozdil
            katsnames = union(katsnames,'AllEl','stable'); %pridam posledni kategorii all
            
            if ~exist('kategorie','var'), kategorie = 1:numel(katsnames); end %kategorie ze ktere chci ziskat hodnoty - odpovida kategoriim z IntervalyResp a pak v CBrainPLot
            %kategorie budou indexy v katsnames
            
            BPD.intervals = intervaly; %budu mit dve pole hodnoty, vybrane kanaly a vsechny kanaly s vybranymi vyznacenyma           
            BPD.katstr = katsnames(kategorie);
            BPD.testname = obj.PsyData.testname; %jmeno testu, aedist, menrot nebo ppa
            BPD.reference = obj.reference; %reference
            BPD.Hf = obj.Hf;%seznam frekvencnich pasem
            celltpl = cell(size(intervaly,1),numel(kategorie)); %size (interval,kat) - templat velikosti pro ostatni promenne
            %jen definice velikosti
            BPD.selCh = celltpl;            
            BPD.NAMES = celltpl;                      
            BPD.MNI = celltpl;
            BPD.VALS = celltpl;
            BPD.NLABELS = celltpl;%sem budu ukladata neurologyLabel od Martina Tomaska  
            BPD.EPI = celltpl; %pridam jeste udaje o epilepticke aktivite, ktera pak muzu pouzit v zobrazeni mozku            
            BPD.filename = basename(iff(isa(obj,'CHilbertMulti'),obj.mfilename,obj.hfilename)); %jmeno zdrojoveho souboru
            BPD.label = iff(isa(obj,'CHilbertMulti'),obj.label,'');  %pokud je tohle instance CHilbertMulti, vezmi z ni label
            %nejdriv udaje pro vsechny elektrody
            for int = 1:size(intervaly,1)
                for kat = 1:numel(kategorie)
                    %TODO 17.9.2018
                    if kat <= size(prumery,3)
                        ich = prumery(:,int, kat) ~= 0; % chci i zaporny rozdil ; aby tam neco bylo                     
                        BPD.VALS{int,kat} = prumery(ich,int,kategorie(kat));
                    else
                        ich = true(size(prumery,1),1);
                        BPD.VALS{int,kat} = zeros(size(prumery,1),1); %zero values for last catefory with all channels
                    end
                        
                    BPD.NAMES{int,kat}= names(ich);
                    BPD.MNI{int,kat} = MNI(ich);            

                    BPD.NLABELS{int,kat}= neurologyLabels(ich);
                    if isfield(obj.CH.H.channels,'seizureOnset')                        
                        BPD.EPI{int,kat} = struct('seizureOnset',{obj.CH.H.channels(ich).seizureOnset},'interictalOften',{obj.CH.H.channels(ich).interictalOften},'rejected',{obj.CH.H.channels(ich).interictalOften});
                    else
                        BPD.EPI{int,kat} = struct('seizureOnset',{nan(numel(ich),1)},'interictalOften',{nan(numel(ich),1)},'rejected',{nan(numel(ich),1)});
                    end
                    if isprop(obj,'plotRCh') && isfield(obj.plotRCh,'selCh') && sum(sum(obj.plotRCh.selCh))>0 ...
                            && ~strcmp(katsnames{kat},'AllEl') %if it is not all channels category
                        BPD.selCh{int,kat} = obj.plotRCh.selCh(ich,:); %novy format z 22.8.2018 - kanaly x marks
                    else
                        BPD.selCh{int,kat} = true(sum(ich),6);
                    end                                                                              
                end
            end            
        end
        function obj = RemoveChannels(obj,channels)       
            %smaze se souboru vybrane kanaly. Kvuli redukci velikost aj
            assert(obj.channels > numel(channels),['the number of existing channels (' num2str(obj.channels) ') is not higher than channels to remove (' num2str(numel(channels)) ')']);           
            keepch = setdiff(1:obj.channels,channels); %channels to keep            
            obj.channels = obj.channels - numel(channels);
            obj.d = obj.d(:,keepch,:);
            if isprop(obj,'mults') && ~isempty(obj.mults), obj.mults = obj.mults(:,keepch); end
            obj.HFreq = obj.HFreq(:,keepch,:,:);
            if isprop(obj,'HFreqEpochs') && ~isempty(obj.HFreqEpochs), obj.HFreqEpochs = obj.HFreqEpochs(:,keepch,:,:); end
            if isprop(obj,'fphaseEpochs') && ~isempty(obj.fphaseEpochs), obj.fphaseEpochs = obj.fphaseEpochs(:,keepch,:,:); end
            if isprop(obj,'frealEpochs') && ~isempty(obj.frealEpochs), obj.frealEpochs = obj.frealEpochs(:,keepch,:,:); end            
            if isprop(obj,'plotRCh') 
                if isfield(obj.plotRCh,'selCh')  && ~isempty(obj.plotRCh.selCh)
                    obj.plotRCh.selCh = obj.plotRCh.selCh(keepch,:);                         
                end
                if isfield(obj.plotRCh,'selChSave') && ~isempty(obj.plotRCh.selChSave)
                    for j = 1:numel(obj.plotRCh.selChSave)
                        obj.plotRCh.selChSave(j).selCh = obj.plotRCh.selChSave(j).selCh(keepch,:);
                    end
                end
            end            
            if isprop(obj,'RjEpochCh') && ~isempty(obj.RjEpochCh), obj.RjEpochCh = obj.RjEpochCh(keepch,:); end
            obj.CH.RemoveChannels(channels);
            obj.els = obj.CH.els; %ty uz se redukuji v CHHeader
            obj.RjCh = obj.CH.RjCh;      
            
            for j = 1:numel(obj.Wp)
                if isfield(obj.Wp(j),'D2')
                    obj.Wp(j).D2 = obj.Wp(j).D2(:,keepch);
                end
                obj.Wp(j).DiEpCh = obj.Wp(j).DiEpCh(keepch,:);
                for k = 1:numel(obj.Wp(j).WpKat)
                    if numel(obj.Wp(j).WpKat{k}) > 0
                        obj.Wp(j).WpKat{k} = obj.Wp(j).WpKat{k}(:,keepch);
                    end
                end
                for k = 1:numel(obj.Wp(j).WpKatBaseline)
                    if numel(obj.Wp(j).WpKatBaseline{k}) > 0
                        obj.Wp(j).WpKatBaseline{k} = obj.Wp(j).WpKatBaseline{k}(:,keepch);
                    end
                end
            end
            disp([ num2str(numel(channels)) ' channels removed']);
        end
        function AppendData(obj,E2, startTimeE2)
            % startTimeE2 - time in sec relative to the stimulus in E2, at which the part of the object E2 should be taken
            assert(size(obj.HFreq,2)==size(E2.HFreq,2) && size(obj.HFreq,3)==size(E2.HFreq,3) && size(obj.HFreq,4)==size(E2.HFreq,4),'same number of channels, frequencies and categories are required');
            assert(isequal(obj.Hf,E2.Hf),'same frequency bands required');
            assert(size(obj.HFreqEpochs,2)==size(E2.HFreqEpochs,2) && size(obj.HFreqEpochs,3)==size(E2.HFreqEpochs,3) && size(obj.HFreqEpochs,4)==size(E2.HFreqEpochs,4),'same number of channels, frequencies and epochs are required');
            assert(isequal(obj.normalization,E2.normalization),'same normalization required');
            if ~exist('startTimeE2','var') || isempty(startTimeE2), startTimeE2 = E2.epochtime(1); end % if no time is given, take the entire period of epoch in E2 including time before stimulus
            AppendData@CiEEGData(obj,E2, startTimeE2); 
            
            iepochtimeE2 = round(E2.epochtime(1:2).*obj.fs);
            istartTimeE2 = round(startTimeE2*obj.fs); % startTime of E2 in samples
            if istartTimeE2 == iepochtimeE2(1)
                iD = 1 : iepochtimeE2(2)-iepochtimeE2(1); % all time points of epoch in E2 including time before stimulus (e.g. in memact in file after delay, it'll take [-1.95, 2])
            else
                iD = abs(iepochtimeE2(1)) + istartTimeE2 : iepochtimeE2(2)-iepochtimeE2(1); % time points in E2 considering start time (might be before of after stimulus, e.g. [-.5, 2] or [.5, 2])
            end
%             if(iepochtimeE2(1) < 0)
%                 iD = abs(iepochtimeE2(1)): iepochtimeE2(2)-iepochtimeE2(1);   %import only part after stimulus of the E2 data
%             else
%                 iD = 1:size(E2.d,1); %if the epoch starts after 0, we want to start the actual epoch start
%             end
            
            obj.HFreq = cat(1,obj.HFreq,E2.HFreq(iD,:,:,:));
            if ~isempty(obj.HFreqEpochs)
                obj.HFreqEpochs = cat(1,obj.HFreqEpochs,E2.HFreqEpochs(iD,:,:,:));
                obj.fphaseEpochs = cat(1,obj.fphaseEpochs,E2.fphaseEpochs(iD,:,:,:));
                obj.frealEpochs = cat(1,obj.frealEpochs,E2.frealEpochs(iD,:,:,:));
            end
            if ~isempty(obj.fphase)
                 obj.fphase = cat(1,obj.fphase,E2.fphase(iD,:,:,:));
            end
        end
    end 
    methods (Static,Access = public)
        function filename2 = filenameE(filename)
            %vraci jmeno souboru s daty tridy CiEEGData
           filename=strrep(filename,'_CHilb',''); %odstranim pripony vytvorene pri save
           filename=strrep(filename,'_CiEEG','');
           filename=strrep(filename,'_CHMult','');
           [pathstr,fname,ext] = CiEEGData.matextension(filename);         
           filename2 = fullfile(pathstr,[fname '_CiEEG' ext]);
        end
        function filename2 = filenameH(filename)
             %vraci jmeno souboru s daty teto tridy
           filename=strrep(filename,'_CHilb',''); %odstranim pripony vytvorene pri save
           filename=strrep(filename,'_CiEEG','');
           filename=strrep(filename,'_CHMult','');
           [pathstr,fname,ext] = CiEEGData.matextension(filename);            
           filename2 = fullfile(pathstr,[fname '_CHilb' ext]);
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
            assert(sum(isnan(filtData))==0, ['hilbertJirka: spatne definovany filtr ' num2str(loF) '-' num2str(hiF) ' Hz']);
             %2) analyticky signal pomoci Hilbertovy transformace:
            tmp = hilbert(filtData);

             %3) gamma amplitude, power
            %gammaAmp = abs(tmp);
            freqPow = abs(tmp).^2; %power je druha mocnicna, 
            % viz https://en.wikipedia.org/wiki/Spectral_density 
            %-  "power" is simply reckoned in terms of the square of the signal,
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
                   obj.PlotEpochs(obj.plotRCh.ch,obj.Wp(obj.WpActive).kats); %vykreslim prumery freq u vsech epoch
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
                case {'divide','slash'} %lomeno na numericke klavesnici - automaticke meritko na ose z - power
                   obj.plotF.ylim = [];
                   obj.PlotResponseFreq( obj.plotF.ch); %prekreslim grafy
                case {'add' ,  'equal','s'}     % + oznaceni kanalu
                   obj.SelChannel(obj.plotF.ch);
                   obj.PlotResponseFreq( obj.plotF.ch); %prekreslim grafy
                case {'numpad6','d'}     % + oznaceni kanalu                   
                   chn2 = obj.plotRCh.selCh( find(obj.plotRCh.selCh>obj.plotF.ch,1) );
                   obj.PlotResponseFreq( iff(isempty(chn2),obj.plotF.ch,chn2) ); %prekreslim grafy
               case {'numpad4','a'}     % + oznaceni kanalu
                   chn2 = obj.plotRCh.selCh( find(obj.plotRCh.selCh<obj.plotF.ch,1,'last') );
                   obj.PlotResponseFreq( iff(isempty(chn2),obj.plotF.ch,chn2) ); %prekreslim grafy
           end
        end       
     end
end

