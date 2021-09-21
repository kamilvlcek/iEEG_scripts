classdef CPlotsN < handle
    %CPlotsN pokrocile funkce na zobrazeni frekvencnich dat, pouziti napr. PN = CPlotsN(E)
    %   Nada Bednarova 2018/04
    
    properties (Access = public)
        E;                      % CiEEGData object
        HOrigData;              % original untransformed EEG data - unepoched
        plotData = struct;      % contains figure and plot info
        plotISPC = struct;        
    end
    
    methods (Access = public)
        
        function obj = CPlotsN(E, HOrigData)
            obj.E = E;            
            if size(obj.E.d,3)>1 %for epoched data
                obj.HOrigData = []; %no original data - used only in plotUnepochedData
            elseif ~exist('HOrigData', 'var') || isempty(HOrigData) %only for non-epoched data
                obj.HOrigData = downsample(E.d,E.decimatefactor); %do funkce plotUnepochedData
            else
                %argument used for non-epoched data
                obj.HOrigData = HOrigData;           
            end
        end
        
        function PlotUnepoched(obj, type, psy, channels, frequencies, limits, timeDelay, startTime, plotGroup)
            % Plots powers and phases of transformed eeg signal for :
            %   * all selected frequencies for one channel
            %   * or one frequency for all channels if called from PlotElectrodeGroup
            % INPUT:
            % type          - 0 = power (HFreq), 1 = original filtered data (obj.freal)
            % psy           - psy structure, aka aedist for stimuli and response time
            % channels      - if called directly, only one channel
            % frequencies   - frequencies to plot (not indexes)
            % limits        - y axis limits, default [0,1000]
            % timeDelay     - how many seconds display on the screen - default 10
            % startTime     - second to start from
            % plotGroup     - only if called from PlotElectrodeGroup, default false
            
            obj.plotData.type = type;
            assert(~isempty(obj.HOrigData),'the field HOrigData with unepoched data cant be empty');
            assert(isa(psy,'struct'),'Prvy parameter musi byt struktura s datami z psychopy'); % kvoli response time
            obj.E.PsyData = CPsyData(psy); 
            
            obj.plotData.iChannels = channels; % pri priamom volani len jeden kanal!
            
            if ~exist('frequencies','var') || isempty(frequencies); obj.plotData.iFreqs = 1:numel(obj.E.Hf); % ak sme nevybrali frekvencie, zobrazia sa vsetky
                else [~,obj.plotData.iFreqs] = intersect(obj.E.Hf,frequencies); end % vypocitam indexy vybranych frekvencii
                
            if ~exist('limits','var') || isempty(limits); limits = [0,1000]; end % set default y axis limits
            obj.plotData.ylimits = limits;
            
            if ~exist('timeDelay','var'); timeDelay = 10; end % number of seconds visible on the screen
            obj.plotData.timeDelay = timeDelay*obj.E.fs;
            
            if ~exist('startTime','var'); startTime = 0; end 
                obj.plotData.iTime = startTime; % initiate time index - second to start from
            
            obj.plotData.f = figure('Name','All Frequencies','Position', [20, 100, 1000, 600]); % init the figure
            set(obj.plotData.f, 'KeyPressFcn', @obj.MovePlotFreqs); % set key press function
            
            if ~exist('plotGroup','var'); obj.plotData.plotGroup = false; else, obj.plotData.plotGroup = plotGroup; end
     
            if ~isfield(obj.plotData,'plotGroup'); obj.plotData.plotGroup = false; end 
            obj.plotData.stimuli = (obj.E.PsyData.P.data(:,obj.E.PsyData.P.sloupce.ts_podnet) - obj.E.tabs(1))*24*3600;
            obj.plotData.responses = (obj.E.PsyData.P.data(:,obj.E.PsyData.P.sloupce.ts_odpoved) - obj.E.tabs(1))*24*3600;
            obj.plotData.colors = {'black','green','red','blue'};
            
            obj.plotUnepochedData();
        end
                
        function PlotElectrodeGroup(obj, iGroup, psy, frequency, limits, timeDelay)
            contacts = obj.E.CH.chgroups{iGroup};
            obj.PlotFrequencies(psy, contacts, frequency, limits, timeDelay, true);
        end
        
        function PlotAllEpochs(obj, iCondition, channel, zlimits)
            % PlotAllEpochs(obj, iCondition, channel,ylimits) - cislo podminky, kanal, rozsah z osy
            % plots all available time x frequency epoch maps for given channel and condition (aedist: 1 = ego, 2 = allo, 0 = red)
            % Nada since 2018/01
            assert(~isempty(obj.E.HFreqEpochs),'soubor s frekvencnimi daty pro epochy neexistuje');
            if(~exist('zlimits','var')), zlimits =[-1 1]; end
            
            condition = obj.E.PsyData.CategoryName(iCondition);  %zjisti jmeno kategorie z jejiho cisla
            figure('Name',[condition, ' - channel ', num2str(channel)], 'NumberTitle', 'off');
            
            time = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.HFreqEpochs,1));
            
            epochs = find(obj.E.PsyData.P.data(:, obj.E.PsyData.P.sloupce.kategorie) == iCondition);
            errors = obj.E.PsyData.GetErrorTrials();
            correct = 0;
            
            subplot_x = ceil(sqrt(length(epochs)/15)*3); % dynamicke rozmery subplotu pre rozny pocet epoch (pre screen cca 3:5)
            subplot_y = ceil(sqrt(length(epochs)/15)*5); % vychadzam z 3x * 5x = pocet_epoch
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT TIMExFREQ OF ALL EPOCHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:length(epochs)
                
               subplot(subplot_x, subplot_y, i);
               imagesc(squeeze(obj.E.HFreqEpochs(:,channel,:,epochs(i)))', 'XData', time, 'YData', obj.E.Hf);
               hold on;
               colormap parula; %aby to bylo jasne u vsech verzi matlabu - i 2010 - vypada v takovem mnozstvi lip nez jet
               caxis(zlimits);
               set(gca,'YDir','normal');
               
               title([num2str(obj.E.PsyData.P.data(epochs(i),5)), '-Epoch (', num2str(epochs(i)), ')'],'FontSize',8);
               
               % plot rejected line
               if obj.PlotRejected(channel, time, epochs(i), sum(errors(epochs(i),:))) %jestli jde o vyrazenou epochu
                  correct = correct +1; % ukladam indexy spravnych epoch
                  if isfield(obj.E.PsyData.P.sloupce,'opakovani')
                      sloupec_opakovani = obj.E.PsyData.P.sloupce.opakovani; %pro AEDIST
                  elseif isfield(obj.E.PsyData.P.sloupce,'opakovani_obrazku')
                      sloupec_opakovani = obj.E.PsyData.P.sloupce.opakovani_obrazku; %pro PPA
                  else
                      sloupec_opakovani = 1; %soubor, skoro vzdy 0
                  end
                  title(sprintf('%d - epoch %d (%d) ', obj.E.PsyData.P.data(epochs(i),sloupec_opakovani), correct, epochs(i)),'FontSize',8) ;
               end
               
               % plot response time 
               response_time = obj.E.PsyData.P.data(epochs(i), obj.E.PsyData.P.sloupce.rt);
               plot([response_time response_time], [obj.E.Hf(1)-10 obj.E.Hf(end)+10],'black','LineWidth',2);               
               hold on;
           end
            
           hold off;
           mtit([condition, ' - channel ', num2str(channel)],'fontsize',20,'color','black');
           
           
           %%%%%%%%%%%%%%%%%%%%%%%% PLOT MEAN FREQUENICES OF ALL EPOCHS %%%%%%%%%%%%%%%%%%%%%%%%%%
           figure('Name','Mean frequencies for all epochs');
           
           correct_epochs = obj.CorrectEpochs(channel, iCondition);
           imagesc(time, 1:length(correct_epochs), squeeze(mean(squeeze(obj.E.HFreqEpochs(:,channel,:,correct_epochs)),2))'); %kreslime jen nevyrazene epochy
           
           title([condition, ' - channel ', num2str(channel)],'fontsize',20,'color','black');
           colormap parula; %aby to bylo jasne u vsech verzi matlabu - i 2016
           colorbar;
           caxis(zlimits);
        end
        
        function fig = PlotFrequencyPower(obj, channel, icondition)
            %funkcia pre vykreslenie time x frequency a frequency x power plotov 
            %pre dany channel a condition (napriklad 0=cervena, 1=vy, 2=znacka, 9=all conditions)
            %priemerne z-scored power cez cely casovy usek pre kazdu frekvenciu
            %a zvlast pre useky pred/po podnete
            assert(~isempty(obj.E.HFreqEpochs),'pole HFreqEpochs s frekvencnimi daty pro epochy neexistuje');
            correct_epochs = obj.CorrectEpochs(channel, icondition); %vyfiltruje len spravne epochy
            
            time = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.HFreqEpochs,1));
            time_before = time(time<=0); %pred podnetom
            
            names = {obj.E.CH.H.channels.ass_brainAtlas}; %cast mozgu
            labels = {obj.E.CH.H.channels.neurologyLabel}; %neurology label
            electrodes = {obj.E.CH.H.channels.name}; %nazov elektrody
            
            [mean_fs, std_fs] = obj.meanZscoredPower(1:length(time),channel,correct_epochs); %priemerna power a stredna chyba priemeru pre kazdu frekvenciu cez celu periodu
            [mean_before, std_before] = obj.meanZscoredPower(1:length(time_before),channel,correct_epochs); %priemerna power a stredna chyba priemeru pre kazdu frekvenciu pred podnetom
            [mean_after, std_after] = obj.meanZscoredPower((length(time_before)+1):length(time),channel,correct_epochs); %priemerna power a str. chyba priemeru pre kazdu frekvenciu po podnete
            
            max_std = max(max([std_fs std_before std_after]));
            mean_min = min(min([mean_fs mean_before mean_after])) - max_std; %min a max means pre zjednotenie ylimits grafov
            mean_max = max(max([mean_fs mean_before mean_after])) + max_std;
            
            fig = figure;
            
            subplot(2,2,1);
            %vykresli time x frequency x z-scored power pre dany channel
            colormap parula; 
            imagesc(squeeze(mean(obj.E.HFreqEpochs(:,channel,:,correct_epochs),4))', 'XData', time, 'YData', obj.E.Hf); 
            set(gca,'YDir','normal');
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            h = colorbar;
            ylabel(h, 'Power'); 
            title([names{channel}, ' ', labels{channel}]);
            
            subplot(2,2,2);
            %vykresli priemernu z-scored power pre kazdu frekvenciu napriec celym casom
            errorbar(mean_fs, std_fs);
            %ylim([mean_min mean_max]);
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title('-1:1s');
            
            subplot(2,2,3);
            %vykresli priemernu z-scored power pre kazdu frekvenciu pred podnetom
            errorbar(mean_before, std_before);
            ylim([mean_min mean_max]);
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title('-1:0s PRED PODNETOM');
            
            subplot(2,2,4);
            %vykresli priemernu z-scored power pre kazdu frekvenciu po podnete
            errorbar(mean_after, std_after);
            ylim([mean_min mean_max]); 
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title('0:1s PO PODNETE');
            
            if icondition == 9
                mtit(sprintf('%s PACIENT %s - CHANNEL %d, ALL CONDITIONS \n (%d epoch) ', electrodes{channel}, obj.E.CH.H.subjName, channel, length(correct_epochs)));
            else
                condition = obj.E.PsyData.CategoryName(icondition);  %zjisti jmeno kategorie z jejiho cisla
                mtit(sprintf('%s PACIENT %s - CHANNEL %d, %s \n (%d epoch) ', electrodes{channel}, obj.E.CH.H.subjName, channel, condition, length(correct_epochs)));
            end
            
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperPosition', [0 0 30 20]);
        end
        
        function fig = PlotAllFrequencyPower(obj, channel)
            %funkcia pre vykreslenie time x frequency a frequency x power plotov 
            %pre dany channel a vsetky conditions
            %priemerne z-scored power cez cely casovy usek pre kazdu frekvenciu
            %a zvlast pre useky pred/po podnete
            assert(~isempty(obj.E.HFreqEpochs),'soubor s frekvencnimi daty pro epochy neexistuje');
            
            time = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.HFreqEpochs,1));
            time_before = time(time<=0); %pred podnetom
            [katnum, katstr] = obj.E.PsyData.Categories();
            kategories = length(katnum);
            hue = 0.8;
            colorskat = {[0 0 0],[0 1 0],[1 0 0],[0 0 1]; [hue hue hue],[hue 1 hue],[1 hue hue],[hue hue 1]};
            
                    
            names = {obj.E.CH.H.channels.ass_brainAtlas}; %cast mozgu
            labels = {obj.E.CH.H.channels.neurologyLabel}; %neurology label
            electrodes = {obj.E.CH.H.channels.name}; %nazov elektrody
            
            for i = 1:kategories
                correct_epochs{i} = obj.CorrectEpochs(channel, katnum(i)); %vyfiltruje len spravne epochy pre danu condition
                [mean_fs{i}, std_fs{i}] = obj.meanZscoredPower(1:length(time),channel,correct_epochs{i}); %priemer + odchylky power pre kazdu frekvenciu cez celu periodu
                [mean_before{i}, std_before{i}] = obj.meanZscoredPower(1:length(time_before),channel,correct_epochs{i}); %priemer + odchylky power pre kazdu frekvenciu pred podnetom
                [mean_after{i}, std_after{i}] = obj.meanZscoredPower((length(time_before)+1):length(time),channel,correct_epochs{i}); %priemer + odchylky power pre kazdu frekvenciu po podnete
            
            end
            
            max_std = max(max([cell2mat(std_fs) cell2mat(std_before) cell2mat(std_after)]));
            mean_min = min(min([cell2mat(mean_fs) cell2mat(mean_before) cell2mat(mean_after)])) - max_std; %min a max means pre zjednotenie ylimits grafov
            mean_max = max(max([cell2mat(mean_fs) cell2mat(mean_before) cell2mat(mean_after)])) + max_std;
            
            zscore_min = min(min(squeeze(mean(obj.E.HFreqEpochs(:,channel,:,:),4))));
            zscore_max = max(max(squeeze(mean(obj.E.HFreqEpochs(:,channel,:,:),4))));
            
            fig = figure;
            for i = 1:kategories
                subplot(2,kategories,i); %vykresli time x frequency x z-scored power pre dany channel
                obj.subplotTimeFrequency(squeeze(mean(obj.E.HFreqEpochs(:,channel,:,correct_epochs{i}),4))', time);
                %caxis([zscore_min zscore_max])
                title(sprintf('%s (%d epoch)', katstr{i}, length(correct_epochs{i})));
                hold on;
            end
            
            subplot(2,kategories,kategories+1); %priemerna z-scored power pre kazdu frekvenciu napriec celym casom (vsetky kategorie)
            for i = 1:kategories
                plotband(obj.E.Hfmean, mean_fs{i}, std_fs{i}, colorskat{2,i});
                hold on;
                h(i) = plot(obj.E.Hfmean, mean_fs{i},'LineWidth',1,'Color',colorskat{1,i});
                hold on;
                ylim([mean_min mean_max]);
                xlabel('Frequency (Hz)');
                ylabel('z-scored power');
                title('-1:1s');
            end
            legend(h, katstr);
            hold off;
            
            subplot(2,kategories,kategories+2); %priemerna z-scored power pre kazdu frekvenciu pred podnetom (vsetky kategorie)
            for i = 1:kategories
                plotband(obj.E.Hfmean, mean_before{i}, std_before{i}, colorskat{2,i});
                hold on;
                h(i) = plot(obj.E.Hfmean, mean_before{i},'LineWidth',1,'Color',colorskat{1,i});
                hold on;
                ylim([mean_min mean_max]);
                xlabel('Frequency (Hz)');
                ylabel('z-scored power');
                title('-1:0s PRED PODNETOM');
            end
            legend(h, katstr);
            hold off;
             
            subplot(2,kategories,kategories+3); %priemerna z-scored power pre kazdu frekvenciu po podnete (vsetky kategorie)
            for i = 1:kategories
                plotband(obj.E.Hfmean, mean_after{i}, std_after{i}, colorskat{2,i});
                hold on;
                h(i) = plot(obj.E.Hfmean, mean_after{i},'LineWidth',1,'Color',colorskat{1,i});
                hold on;
                ylim([mean_min mean_max]);
                xlabel('Frequency (Hz)');
                ylabel('z-scored power');
                title('0:1s PO PODNETE');
            end
            legend(h, katstr);
            hold off;
      
            mtit(sprintf('%s PACIENT %s - CHANNEL %d \n %s - %s \n ', electrodes{channel}, obj.E.CH.H.subjName, channel, names{channel}, labels{channel}));

            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperPosition', [0 0 55 30]);
        end
        
        function PlotMovingEpochs(obj, channels, iChannel)
            %pro zadane channels kresli graf cas x frekvence pro kazdou epochu zvlast
            %sipkami se da prochazet pres epochy a kanaly - private function MovePlotEpochs
            %Nada since 2018/01
            if ~exist('iChannel', 'var'); iChannel = 1; end
            if ~exist('channels', 'var') || isempty(channels); channels = 1:obj.E.channels; end 
            assert(~isempty(obj.E.HFreqEpochs),'soubor s frekvencnimi daty pro epochy neexistuje');
            obj.plotData.f = figure('Name','All Epochs','Position', [20, 100, 1000, 600]);
            obj.plotData.channels = channels;
            set(obj.plotData.f, 'KeyPressFcn', @obj.MovePlotEpochs);
            
            obj.plotData.iChannel = iChannel; % initiate channel index
            obj.plotData.iEpoch = 1; % initiate epoch index
            obj.plotData.T = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.HFreqEpochs,1)); % time
            obj.plotData.rejectedEpochs = obj.E.PsyData.GetErrorTrials(); % get rejected epoch trials
            
            % calculate zlimits for all channels
            obj.plotData.zlimits = zeros(length(channels),2);
            for ch = 1:length(channels)
                obj.plotData.zlimits(ch, :) = obj.getZlimits(obj.plotData.channels(ch));
            end
            
            obj.plotEpochData();
        end   
        
        function [itpc, itpc_p, itpc_pmean, n_epoch] = CalculateITPC(obj, conditions, channels)
            % itpc = ch x condition x time x fq
            % itpc_p ch x condition x time x fq - p values
            % itpc_pmean ch x condition x time - p values z mean itpc pres frekvence
            if ~exist('channels', 'var') 
                channels = 1:obj.E.channels; 
            end
            if ~exist('conditions', 'var')
                conditions = [obj.E.PsyData.P.strings.podminka{:,2}];
            end
            n_fq = length(obj.E.Hf);
            n_cond = length(conditions);
            itpc = zeros(obj.E.channels,n_cond,obj.E.samples, n_fq);
            itpc_p = zeros(obj.E.channels,n_cond,obj.E.samples, n_fq); % p-value vsetkych frekvencii
            itpc_pmean = zeros(obj.E.channels,n_cond,obj.E.samples); % p-value priemeru cez vsetky frekvencie
            n_epoch = zeros(1,n_cond);
            
            for channel = channels
                %display(['Channel ' num2str(channel)]);
                for cond = 1:length(conditions)
                    if iscell(conditions)
                        icondition = conditions{cond};
                    else
                        icondition = conditions(cond);    
                    end
                    iEp = obj.CorrectEpochs(channel, icondition);
                    n_epoch(cond) = length(iEp);
                    for fq = 1:n_fq
                        itpc(channel, cond, :, fq) = abs(mean(exp(1i*squeeze(obj.E.fphaseEpochs(:,channel,fq,iEp))),2));
                        itpc_p(channel, cond, :, fq) = exp(-n_epoch(cond) * squeeze(itpc(channel, cond,:,fq)).^2);  %hladina siginifikance je zavisla na poctu epoch
                    end
                    itpc_pmean(channel, cond, :) = exp(-n_epoch(cond) * mean(squeeze(itpc(channel, cond,:,:)),2).^2);
                end
            end
        end
        
        function [ditpc, ditpc_p, ditpc_pmean] = CalculateITPCdiffs(obj, diffs,itpc, channels) %#ok<INUSL>
            %stejne jako CalculateITPC ale pro rozdil 
            if ~exist('channels', 'var') || isempty(channels)
                channels = 1:obj.E.channels; 
            end
            if ~exist('diffs', 'var')  || isempty(diffs)
                diffs = [[0 1]; [0 2]; [1 2]]; 
            end
            if ~exist('itpc','var') %muzu dat jako dalsi parametr uz spocitane itpc
                [itpc, ~, ~] = obj.CalculateITPC(channels); %#ok<ASGLU>
            end
            n_fq = length(obj.E.Hf);
            ditpc = zeros(obj.E.channels,length(diffs),obj.E.samples, n_fq);
            ditpc_p = zeros(obj.E.channels,length(diffs),obj.E.samples, n_fq); % p-value vsetkych frekvencii
            ditpc_pmean = zeros(obj.E.channels,length(diffs),obj.E.samples); % p-value priemeru cez vsetky frekvencie
                        
            [itpc, ~, ~, ~] = obj.CalculateITPC(channels);
            for channel = channels
                %display(['Channel ' num2str(channel)]);
                for diff = 1:size(diffs, 1)
                    cond1 = diffs(diff,1);
                    cond2 = diffs(diff,2);
                    ditpc(channel, diff, :, :) = squeeze(itpc(channel, cond1+1, :, :) - itpc(channel, cond2+1, :, :)); % time x freq

                    n1 = length(obj.CorrectEpochs(channel, cond1));
                    n2 = length(obj.CorrectEpochs(channel, cond2));
                    
                    ditpc_p(channel, diff, :, :) = 2*normpdf(abs(atanh(squeeze(ditpc(channel, diff, :, :)))/sqrt((1/(n1-3))+(1/(n2-3)))));
                    ditpc_pmean(channel, diff, :) = 2*normpdf(abs(atanh(mean(squeeze(ditpc(channel, diff, :, :)),2))/sqrt((1/(n1-3))+(1/(n2-3)))));
                end
            end
        end
        
        
        function [citpc, citpc_p, citpc_pmean,combinations] = CalculateITPCcontrasts(obj, contrasts, itpc, n, channels)
            %{ 
             contrasts is either an array or cell array
             [0 1 2] or {2 [1 2] [2 1 0]} etc.
             citpc = ch x contrast1 x contrast2 x time x fq
             citpc_p ch x contrast1 x contrast2 x time x fq - p values
             citpc_pmean ch x contrast1 x contrast2 x time - p values z mean itpc pres frekvence
            %}
            if ~exist('channels', 'var') 
                channels = 1:obj.E.channels; 
            end
            if ~exist('contrasts', 'var') 
                contrasts = [obj.E.PsyData.P.strings.podminka{:,2}];
            end
            if ~exist('itpc','var') || ~exist('n','var')
                [itpc, ~, ~, n] = obj.CalculateITPC(channels, contrasts); 
            end
            
            contrasts_num = length(contrasts);
            combinations = nchoosek(1:contrasts_num,2); % returns all combinations of 2 values from 1:length(contrasts)
                                                        % = the indexes of contrast combinations to be calculated
     
            citpc = nan(obj.E.channels,contrasts_num,contrasts_num,obj.E.samples, length(obj.E.Hf));
            citpc_p = nan(obj.E.channels,contrasts_num,contrasts_num,obj.E.samples, length(obj.E.Hf)); % p-value vsetkych frekvencii
            citpc_pmean = nan(obj.E.channels,contrasts_num,contrasts_num,obj.E.samples); % p-value priemeru cez vsetky frekvencie
            
            for channel = channels
                %display(['Channel ' num2str(channel)]);
                for comb = 1:length(combinations)
                                        
                    cond1 = combinations(comb,1); % index of the first element in the contrast
                    cond2 = combinations(comb,2); % index of the second element in the contrast
                    citpc(channel, cond1, cond2, :, :) = squeeze(itpc(channel, cond1, :, :) - itpc(channel, cond2, :, :)); % time x freq
                    
                    citpc_p(channel, cond1, cond2, :, :) = 2*normpdf(abs(atanh(squeeze(citpc(channel, cond1, cond2, :, :)))/sqrt((1/(n(cond1)-3))+(1/(n(cond2)-3)))));
                    citpc_pmean(channel, cond1, cond2, :) = 2*normpdf(abs(atanh(mean(squeeze(citpc(channel, cond1, cond2, :, :)),2))/sqrt((1/(n(cond1)-3))+(1/(n(cond2)-3)))));
               
                end
            end
            
        end
        
        %ISPC functions
        function [ispc] = ISPCCalculate(obj,channels,fdr)
            % inter site phase coherence = PLV - according to MikeXCohen
            % currently only: ISPC-trial, for one group of channels, for all computed frequencies, for all time points
            % ispc = time x fq
            if ~exist('channels','var') || isempty(channels), channels = 1:obj.E.channels; end %default is all channels
            if ~exist('fdr','var'), fdr = 1; end %the more strict method
            assert(isprop(obj.E,'fphaseEpochs') && ~isempty(obj.E.fphaseEpochs), 'fphaseEpochs must be computed to use ISPC');
            n_epochs = size(obj.E.fphaseEpochs,4); %number of epochs
            ispc = nan(numel(channels), numel(channels), size(obj.E.fphaseEpochs,1),size(obj.E.fphaseEpochs,3)); %chn x chn x time x fq            
            
            %this is too large field to be computed for 100 channels: Requested 100x100x64x46x650 (142.6GB) array 
            %ispc_all = nan(size(obj.E.fphaseEpochs,1),size(obj.E.fphaseEpochs,3),n_epochs);%chn x chn x time x fq x epochs                           
            
            ispc_p_min = nan(numel(channels), numel(channels));
            fprintf('ISPCCalculate:    ');
            for ichn1 = 1:numel(channels)
                if sum(obj.E.CH.RjCh==channels(ichn1))>0, continue; end %if the channel is rejected, co not compute
                fprintf('\b\b\b%3i',ichn1); %delete previous three characters and print the %i in three characters
                for ichn2 = ichn1+1:numel(channels)    
                    if sum(obj.E.CH.RjCh==channels(ichn2))>0, continue; end
                    if ~isa(obj.E, 'CHilbertMulti') || find(obj.E.els >= channels(ichn1),1) == find(obj.E.els >= channels(ichn2),1) %the find expression get the pacient number
                        %compute only for single pacient or for CHilbertMulti if the pacient is the same
                        fphaseEpochs1 = squeeze(obj.E.fphaseEpochs(:,channels(ichn1),:,:)); %time x freq x epochs
                        fphaseEpochs2 = squeeze(obj.E.fphaseEpochs(:,channels(ichn2),:,:));
                        %the cycles over frequencies and time removed - all computed at once in the following commands
                        ispc_all = exp(1i*(fphaseEpochs1 - fphaseEpochs2)); %time x freq x epochs - complex number
                        ispc_ch = abs(mean(ispc_all,3)); %abs = lenght of the average vector ,%time x freq , mean over epochs
                        ispc(ichn1,ichn2,:,:) = ispc_ch;
                        ispc_p0 = CEEGStat.ISPCBaseline(squeeze(ispc(ichn1,ichn2,:,:)),n_epochs, obj.E.epochtime,obj.E.baseline,obj.E.fs,0); %significance relative to mean baseline ispc, no fdr correction for now
                        if ~exist('ispc_p','var') %easiest way how to find the first dimension of ispc_p0
                            ispc_p = nan(numel(channels), numel(channels),size(ispc_p0,1),size(ispc_p0,2));
                        end
                        ispc_p(ichn1,ichn2,:,:) = ispc_p0;                       
                    end
                end                
            end
            fprintf('\n');
            if fdr>0
                iispc_p = ~isnan(ispc_p); %fdr correction over all non nan values (=computed actually)
                pvals = ispc_p(iispc_p); %this makes a simple array (from original 4D array)
                if fdr == 2, method = 'dep'; else method = 'pdep'; end %#ok<SEPEX>
                [~, ~, adj_p]=fdr_bh(pvals,0.05,method,'no'); %dep is more sctrict than pdep 
                ispc_p(iispc_p) = adj_p; %overwrite the corrected values by FDR corrected - only the original non-nan values
            end
            for ichn1 = 1:numel(channels)                
                for ichn2 = ichn1+1:numel(channels) 
                    ispc_p0 = squeeze(ispc_p(ichn1,ichn2,:,:));
                    iispc_p = ~isnan(ispc_p0);
                    if sum(sum(iispc_p))>0
                        ispc_p_min(ichn1,ichn2)=min(ispc_p0(iispc_p)); %find the most significan ispc value for each combination of channels                                        
                    end
                end
            end
            obj.plotISPC.ispc = ispc;            
            obj.plotISPC.ispc_all = ispc_all; %the angles computed for the last combination of channels
            obj.plotISPC.ispc_p = ispc_p;
            obj.plotISPC.fdr = fdr;
            obj.plotISPC.channels = channels;            
            obj.plotISPC.ispc_p_min = ispc_p_min;             
        end
        
        function ISPCPlot(obj,chnpair,TimeHf, updateS, ISPCPlotSelection)
            %plot all ispc values  for a specific pair of channels from computed values            
            % chnpair - pair of the channels (or rather their xy coordinates in the plot) to be higlited on the left top plot of all channels, and shown an all other figure
            % TimeHf - point in the time-frequency space highlighted on the left top plot - 
            assert(isfield(obj.plotISPC,'ispc')&& ~isempty(obj.plotISPC.ispc), 'no ispc values are computed');
            if ~exist('chnpair','var')  || isempty(chnpair) 
                if isfield(obj.plotISPC,'chnpair') && ~isempty(obj.plotISPC.chnpair)
                    chnpair = obj.plotISPC.chnpair;                     
                else
                    chnpair = 1:2; %plotISPC.chnpair is position in sortorder, plotISPC.sortorder are posisions in plotISPC.ispc ans also in  plotISPC.channels.  plotISPC.channels are actual channel numbers               
                    obj.plotISPC.chnpair = chnpair;
                end
            else
                obj.plotISPC.chnpair = chnpair; %store the state of this figure
            end
            if isfield(obj.plotISPC, 'sortorder')
                chn1 = obj.plotISPC.sortorder(chnpair(1));
                chn2 = obj.plotISPC.sortorder(chnpair(2));
            else
                chn1 = chnpair(1);
                chn2 = chnpair(2);
            end
                
            if isnan(max(max(squeeze(obj.plotISPC.ispc(chn1,chn2,:,:))))) %if the order of the channels is reversed - the ispc values are stored only once
                chn1 = obj.plotISPC.sortorder(chnpair(2));
                chn2 = obj.plotISPC.sortorder(chnpair(1));
            end
            if ~exist('TimeHf','var') || isempty(TimeHf)
                if isnan(obj.plotISPC.ispc_p_min(chn1,chn2))
                    TimeHf = [round(size(obj.plotISPC.ispc,3)/2) round(size(obj.plotISPC.ispc,4)/2)]; % point in the middle
                else
                    [iTime,iHf] = find(squeeze(obj.plotISPC.ispc_p(chn1,chn2,:,:)) == obj.plotISPC.ispc_p_min(chn1,chn2));
                    TimeHf = [iTime(1) iHf(1)];
                end
            end 
            obj.plotISPC.TimeHf = TimeHf; %store the state of this figure   
            if ~exist('updateS','var') || isempty(updateS) , updateS = [1 0 0 0]; end
            if ~exist('ISPCPlotSelection','var'), ISPCPlotSelection = 0; end
            assignhandle = false;
            if ~isfield(obj.plotISPC,'fh_ispc') || isempty(obj.plotISPC.fh_ispc) || ~ishghandle(obj.plotISPC.fh_ispc)
                obj.plotISPC.fh_ispc = figure('Name','ISPCPlot');   
                obj.plotISPC.subph = {}; %handles to subplots
                obj.plotISPC.annoth = {}; %handles to annotations
                obj.plotISPC.plotsig = [0.95 1]; % limits of significance plots
                obj.plotISPC.onlysub2 = 0; % plot only second chart with an overview of all channels
                obj.plotISPC.clf = 0; % clear the whole plot to plot all subplots after only second chart was plotted                
                obj.plotISPC.sortorder = 1:size(obj.plotISPC.ispc,1); %originally channels sorted by default
                obj.plotISPC.sortby = 0; %sorted originally, 1=sorted by brainlabel
                obj.plotISPC.els = obj.E.CH.els; %borders of the electrodes/pacients
                obj.plotISPC.aplha = 0.4; %transparency for non signif IPSC values 1=non transparent, 0= completely transparent
                updateS = [1 1 1 1]; %create all subplots
                assignhandle = true;
            else
                figure(obj.plotISPC.fh_ispc);                             
                if ~obj.plotISPC.onlysub2 && ~obj.plotISPC.clf
                    for iu = 1:numel(updateS) %clear individual subplots
                        if updateS(iu)
                            cla( obj.plotISPC.subph{iu})
                            delete(obj.plotISPC.subph{iu});                        
                            if numel(obj.plotISPC.annoth()) >= iu && ~isempty(obj.plotISPC.annoth{iu})
                                delete(obj.plotISPC.annoth{iu});
                            end
                        end                        
                    end
                else
                    clf; %vymazu obrazek    
                    if ~obj.plotISPC.onlysub2 && obj.plotISPC.clf
                        updateS = [1 1 1 1]; %restore all subplots
                        obj.plotISPC.clf = 0;
                    end
                end
            end   
            
            T = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1));
            F =  obj.E.Hfmean;                                 
            
            %subplot with ISPC values from one pair of channels
            if updateS(1) 
                if ~obj.plotISPC.onlysub2
                    obj.plotISPC.subph{1} = subplot(2,2,1);  
                end
                ispc=squeeze(obj.plotISPC.ispc(chn1,chn2,:,:))';
                im1=imagesc(T,F,ispc); %has to transpose, imagesc probably plots first dimension in rows 
                axis xy; %make the y axis increasing from bottom to top
                colormap parula;
                colorbar;
                title(sprintf(' channels %i, %s - %i, %s', ...
                    obj.plotISPC.channels(chn1),obj.E.CH.H.channels(obj.plotISPC.channels(chn1)).name,...
                    obj.plotISPC.channels(chn2),obj.E.CH.H.channels(obj.plotISPC.channels(chn2)).name));
                xlabel([num2str(T(TimeHf(1))) ' sec']);
                ylabel([num2str(F(TimeHf(2))) ' Hz']);
                hold on;
                plot(T(TimeHf(1)),F(TimeHf(2)),'o','MarkerFaceColor','red');            
                obj.plotISPC.annoth{1} = annotation(obj.plotISPC.fh_ispc,'textbox',[0 .9 .4 .1],'String',num2str(obj.plotISPC.ispc(chn1,chn2, TimeHf(1),TimeHf(2))), 'EdgeColor', 'none');
            end
            
            %subplot with ISPC values from one frequency of one pair of channels
            if updateS(3) 
                if ~obj.plotISPC.onlysub2
                    obj.plotISPC.subph{3} = subplot(2,2,3); 
                end
                plot(T,squeeze(obj.plotISPC.ispc(chn1,chn2, :,TimeHf(2))));
                hold on;
                plot(T(TimeHf(1)),obj.plotISPC.ispc(chn1,chn2, TimeHf(1),TimeHf(2)),'o');
                title([num2str(F(TimeHf(2))) ' Hz']);
                maxispc = max(max(squeeze(obj.plotISPC.ispc(chn1,chn2,:,:))));
                ylim([0 iff(isnan(maxispc),1,maxispc)]);            
            end
            
            %subplot with minimum ISPC significance across all channels
            if updateS(2)                
                if ~obj.plotISPC.onlysub2
                    obj.plotISPC.subph{2} = subplot(2,2,2); 
                end
                ispc_p_min = obj.plotISPC.ispc_p_min(obj.plotISPC.sortorder,obj.plotISPC.sortorder); %values sorted according to sortorder                
                [x,y]=find(isnan(ispc_p_min));
                ispc_p_min(sub2ind(size(ispc_p_min),x,y)) = ispc_p_min(sub2ind(size(ispc_p_min),y,x)); %make the matrix symmetrical - each value will be there twice 
                im2=imagesc(obj.plotISPC.channels,obj.plotISPC.channels,1-ispc_p_min',obj.plotISPC.plotsig); 
                imAlpha=ones(size(ispc_p_min'));  %https://www.mathworks.com/matlabcentral/answers/81938-set-nan-as-another-color-than-default-using-imagesc
                imAlpha(isnan(ispc_p_min'))=0;
                im2.AlphaData = imAlpha; %set transparency for nan values
                set(gca,'color',0*[1 1 1]); %black background
                axis xy; colorbar; 
                hold on;
                plot(chnpair(1),chnpair(2),'o','MarkerFaceColor','red');
                for e=1:numel(obj.plotISPC.els)-1
                    line([obj.plotISPC.els(e)+.5 obj.plotISPC.els(e)+.5], [1 numel(obj.plotISPC.channels)]);
                    line([1 numel(obj.plotISPC.channels)],[obj.plotISPC.els(e)+.5 obj.plotISPC.els(e)+.5]);
                end
                if obj.plotISPC.sortby ==1 %if ordered by brain labels
                    labels=unique({obj.E.CH.brainlabels(:).label});
                    els = [1 obj.plotISPC.els(1:end-1)+1];
                    for e=1:numel(els)
                        text(els(e)+5,els(e)+2,labels{e},'Color','red');                        
                    end
                end
                t1 = sprintf('channel %i (%s) x %i (%s)',obj.plotISPC.channels(chn1),obj.E.CH.H.channels(obj.plotISPC.channels(chn1)).name,...
                    obj.plotISPC.channels(chn2), obj.E.CH.H.channels(obj.plotISPC.channels(chn2)).name);
                t2 = [obj.E.CH.H.channels(obj.plotISPC.channels(chn1)).neurologyLabel '-' obj.E.CH.H.channels(obj.plotISPC.channels(chn2)).neurologyLabel '; ' ...
                    obj.E.CH.brainlabels(obj.plotISPC.channels(chn1)).label '-' obj.E.CH.brainlabels(obj.plotISPC.channels(chn2)).label]; 
                title( iff(obj.plotISPC.onlysub2,{t1,t2},t1)); %two line title for this large plot
            end
            
            %subplot with significance of ISPC of this pair of channels
            if updateS(4) 
                if ~obj.plotISPC.onlysub2
                    obj.plotISPC.subph{4} = subplot(2,2,4); 
                end
                T4 = linspace(0, obj.E.epochtime(2), size(obj.plotISPC.ispc_p,3));
                ispc_p = 1-squeeze(obj.plotISPC.ispc_p(chn1,chn2,:,:))';
                im4=imagesc(T4,F,ispc_p,obj.plotISPC.plotsig);  %0.5
                imAlpha=ones(size(ispc_p));  %https://www.mathworks.com/matlabcentral/answers/81938-set-nan-as-another-color-than-default-using-imagesc
                imAlpha(ispc_p<obj.plotISPC.plotsig(1))=obj.plotISPC.aplha;
                im4.AlphaData = imAlpha; %set transparency for nan values
                imAplhaISPC = cat(2,repmat(obj.plotISPC.aplha,size(imAlpha,1),size(ispc,2)-size(ispc_p,2)),imAlpha);
                im1.AlphaData = imAplhaISPC;
                axis xy;                  
                colorbar;
                hold on;
                plot(T(TimeHf(1)),F(TimeHf(2)),'o','MarkerFaceColor','red');  
                title(sprintf('channel %i (%s) x %i (%s)',obj.plotISPC.channels(chn1),obj.E.CH.H.channels(obj.plotISPC.channels(chn1)).neurologyLabel,...
                    obj.plotISPC.channels(chn2), obj.E.CH.H.channels(obj.plotISPC.channels(chn2)).neurologyLabel));
            end
                    
            %another optional plot all phase angles for selected time and freq
            if ~isempty(TimeHf(2)) && ISPCPlotSelection && chnpair(1) == obj.plotISPC.channels(end-1) && chnpair(2) == obj.plotISPC.channels(end)
                
                ispc_sel = squeeze(obj.plotISPC.ispc_all(TimeHf(1),TimeHf(2),:));                 
                if ~isfield(obj.plotISPC,'fh_ispc_sel') || isempty(obj.plotISPC.fh_ispc_sel) || ~ishghandle(obj.plotISPC.fh_ispc_sel)
                    obj.plotISPC.fh_ispc_sel = figure('Name','ISPCPlot Selection');                
                else
                    figure(obj.plotISPC.fh_ispc_sel);
                    clf; %vymazu obrazek                 
                end        
                circle(1,0,0,'k');
                axis equal;
                hold on;
                line([zeros(size(ispc_sel))' ; real(ispc_sel)'],[zeros(size(ispc_sel))' ;  imag(ispc_sel)' ]);
                title(['ispc = ' num2str(obj.plotISPC.ispc(chn1,chn2, TimeHf(1),TimeHf(2)))]);                
            end 
                        
            %figure(obj.plotISPC.fh_ispc); %put the figure to top
            if assignhandle
                methodhandle = @obj.hybejPlotISPC;
                set(obj.plotISPC.fh_ispc,'KeyPressFcn',methodhandle);
                set(obj.plotISPC.fh_ispc, 'WindowButtonDownFcn', @obj.hybejPlotISPCClick);
            end
        end
                    
        function ISPCSave(obj)
            assert(~isempty(obj.E),'no original CiEEGData data exist');
            assert(~isempty(obj.plotISPC),'no plotISPC data exist');
            fname = obj.filenameISPC(obj.E.filename);
            plotISPC = obj.plotISPC; %#ok<PROP>
            plotISPC.fh_ispc = []; %#ok<PROP>
            plotISPC.fh_ispc_sel = []; %#ok<PROP>
            plotISPC.annoth = {};%#ok<PROP>
            plotISPC.subph = {};%#ok<PROP>
            save(fname,'plotISPC','-v7.3'); 
            disp(['saved to ' fname ]);
        end
        function ISPCLoad(obj)
            assert(~isempty(obj.E),'no original CiEEGData data exist');
            fname = obj.filenameISPC(obj.E.filename);
            if exist(fname,'file')
                V = load(fname);
                obj.plotISPC = V.plotISPC;
                disp(['loaded ' fname ]);
            else
                disp(['not found: ' fname ]);
            end
        end
        
        
        function fig = PlotITPCsig(obj, ch, sig, diffs)
            % itpc (ch, cond, time, fq)
            % itpc_p (ch, cond, time)
            
            if ~exist('diffs', 'var') 
                diffs = [[0 1]; [0 2]; [1 2]]; 
            end
            p_val = 0.05;
            names = {obj.E.CH.H.channels.ass_brainAtlas}; %cast mozgu
            labels = {obj.E.CH.H.channels.neurologyLabel}; %neurology label
            electrodes = {obj.E.CH.H.channels.name}; %nazov elektrody
            time = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1));
            
            fig = figure; colormap(jet)
            [itpc, ~, p] = obj.CalculateITPC(ch);
            [ditpc, ~, dp] = obj.CalculateITPCdiffs(ch, diffs);
            
            conditions = [obj.E.PsyData.P.strings.podminka{:,2}];
            numCond = length(conditions);
            for cond = 1:numCond
                % subplot individual conditions
                subplot(2,numCond, cond);
                if ~sig % subplot time x frequency for every condition
                    imagesc(time, obj.E.Hf, squeeze(itpc(ch, cond, :, :))');
                    xlabel('Time (s)'); ylabel('Frequency (Hz)'); colorbar;
                    set(gca,'YDir','normal'); caxis([0 0.7]);
                    title([obj.E.PsyData.P.strings.podminka{cond,1}]);
                    hold on;
                else % subplot time x (p-value of mean itpc across frequencies) for every condition
                    plot(time, mean(squeeze(itpc(ch, cond, :, :)),2), 'g');
                    hold on;
                    plot(time,squeeze(p(ch, cond, :)));
                    ylim([0 1]);
                    xlabel('Time (s)'); ylabel('ITPC p-value');
                    title([obj.E.PsyData.P.strings.podminka{cond,1} '(' num2str(length(obj.CorrectEpochs(ch, obj.E.PsyData.P.strings.podminka{cond,2}))) ')']);
                    hold on;
                    plot(time, p_val*ones(size(time)),'r');
                    hold on;
                end
            end
            
            for diff = 1:size(diffs, 1)
                subplot(2, numCond, diff + numCond);
                if ~sig % subplot time x frequency for condition differences
                    imagesc(time, obj.E.Hf, abs(squeeze(ditpc(ch, diff, :, :)))');
                    xlabel('Time (s)'); ylabel('Frequency (Hz)'); colorbar;
                    set(gca,'YDir','normal'); caxis([0 0.3]);
                    title([obj.E.PsyData.P.strings.podminka{diffs(diff,1)+1,1} '-' obj.E.PsyData.P.strings.podminka{diffs(diff,2)+1,1}]);
                    hold on;
                else % subplot time x (p-value of mean itpc across frequencies) for every condition difference
                    plot(time,abs(mean(squeeze(ditpc(ch, diff, :, :)),2)), 'g'); hold on;
                    plot(time,squeeze(dp(ch, diff, :))); ylim([0 1]);
                    xlabel('Time (s)'); ylabel('ITPC p-value');
                    title([obj.E.PsyData.P.strings.podminka{diffs(diff,1)+1,1} '-' obj.E.PsyData.P.strings.podminka{diffs(diff,2)+1,1}]);
                    hold on;
                    plot(time, p_val*ones(size(time)),'r');
                    hold on;
                end
            end
            mtit(sprintf('%s - Channel %d; %s - %s \n', electrodes{ch}, ch, labels{ch}, names{ch}));
            
        end
                
        function fH = PlotITPCallEpochs(obj, channel, iHf, sortbycondition)
            % modro zlte plots synchronizacie fazy image            
            if ~exist('iHf','var'), iHf = 1:numel(obj.E.Hfmean); end %all frequencies by default
            if ~exist('sortbycondition','var'), sortbycondition = 0; end %if to sorf epochs by condition
            assert(numel(channel)<=3 || numel(iHf)==1, 'needs to plot either one channel or one frequency');
            fH= figure('Name','PlotITPCallEpochs');
            time = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1));
            for ich = 1:numel(channel)
                for iFq = 1:numel(iHf)                    
                    if numel(channel) == 1 %show multiple frequencies for only one channel
                        subplot(2,ceil(length(iHf)/2),iFq);
                        xoff = iff( rem(ceil(length(iHf)/2), 2)==1  ,.3,0); %for only even number of columns of subplots - move the mtit a bit
                    elseif numel(channel) == 2 || numel(channel) == 3 %compare multiple frequencies in 2 or 3 channels
                        subplot(numel(channel),length(iHf),(ich-1)*numel(iHf)+iFq); 
                        xoff = iff( rem(length(iHf), 2)==1  ,.3,0);  %for only even number of columns of subplots - move the mtit a bit
                    else %show one frequency for multiple channels
                        subplot(2,ceil(numel(channel)/2),ich);
                        xoff = iff( rem(ceil(numel(channel)/2),2)==1,.3,0); %for only even number of columns of subplots - move the mtit a bit
                    end                    
                    % sort epochs by condition
                    epochsByCond = [[obj.E.epochData{:,2}]' (1:length(obj.E.epochData))']; % colmuns: stimulus condition, epoch no               
                    if sortbycondition
                        epochsByCond = sortrows(epochsByCond,1); %sort epochs by stimulus condition
                    end 
                    phases = squeeze(obj.E.frealEpochs(:,channel(ich),iHf(iFq),:))'; %epochs x time - the filtered signal
                    image(phases(epochsByCond(:,2),:), 'XData', time); % yellow = signal above zero, blue = signal below zero. 
                    hold on;
                    set(gca,'YDir','normal')
                    xlabel('Time (s)'); ylabel('Epoch');
                    if numel(channel) == 1
                        title([num2str(round(obj.E.Hfmean(iHf((iFq))),2)), ' Hz']); %title for this subplot
                    elseif numel(channel) == 2 || numel(channel) == 3
                        title([num2str(round(obj.E.Hfmean(iHf((iFq))),2)), ' Hz']); 
                    else
                        title([num2str(channel(ich)), obj.E.CH.H.channels(channel(ich)).name]); %title for this subplot
                    end
                    if iFq == 1 && (numel(channel) == 2 || numel(channel) == 3)
                        ylabel([num2str(channel(ich)), obj.E.CH.H.channels(channel(ich)).name]);
                    end
                end                
            end           
            if numel(channel) == 1
                mtit(sprintf('Channel %d, %s', channel(ich),obj.E.CH.H.channels(channel(ich)).name),'xoff',xoff);
            elseif numel(channel) >3 
                mtit(sprintf('Fq %2.2d Hz', round(obj.E.Hfmean(iHf((iFq))),2)),'xoff',xoff);
            end
            hold off;
        end
        
        function fig = PlotITPCallEpochsFiltered(obj, channel, iFq)
            fig = figure;
            time = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1));
            conditions = [obj.E.PsyData.P.strings.podminka{:,2}];
            ncond = length(conditions);
            
            for cond = 1:ncond
                subplot(2,ncond,cond);
                icondition = conditions(cond);
                iEp = obj.CorrectEpochs(channel, icondition);
                plot(time,squeeze(obj.E.frealEpochs(:,channel,iFq,iEp)), 'LineWidth',0.0001);
                title(obj.E.PsyData.P.strings.podminka{cond,1});
                hold on;
                
                subplot(2,ncond,cond+ncond);
                phases = squeeze(obj.E.frealEpochs(:,channel,iFq,iEp))';
                image(phases, 'XData', time);
                hold on;
                set(gca,'YDir','normal')
                xlabel('Time (s)'); ylabel('Epoch');
            end
            hold off;
            xlabel('Time (s)');
            mtit(sprintf('Channel %d , %0.1f Hz \n', channel, obj.E.Hf(iFq)));
        end
        
    end
    
    methods (Static,Access = public)        
        function filename2 = filenameISPC(filename)
             %vraci jmeno souboru s daty teto tridy
           filename=strrep(filename,'_CHilb',''); %odstranim pripony vytvorene pri save
           filename=strrep(filename,'_CiEEG','');
           filename=strrep(filename,'_CHMult','');
           [pathstr,fname,ext] = CiEEGData.matextension(filename);            
           filename2 = fullfile(pathstr,[fname '_ISPC' ext]);
        end
        
    end
    methods (Access = private)
        
        function itpc = ITPCpermute(n_permutes, itpcOrig, eegphase)
            % itpcOrig (time x fq)
            % eegphase (time x fq x epoch)
            eegtemp = zeros(size(eegphase));
            eegperm = zeros(n_permutes, size(eegphase,1), size(eegphase,2));
            
            for permi = 1:n_permutes
                for epoch = 1:size(eegphase, 3)
                    cutpoint = randsample(2:size(eegphase,1)-2, 1);
                    eegtemp(:,:,epoch) = eegphase([cutpoint:end 1:cutpoint-1], :, epoch);
                end
                eegperm(permi,:,:) = squeeze(abs(mean(eegtemp, 3)));
            end
            % pre kazdy prvok mapy sa spocita z-score cez vsetky permutacie
            zmap = (itpcOrig-squeeze(mean(eegperm,1))) ./ squeeze(std(eegperm,1));
            itpc = itpcOrig;
            itpc(abs(zmap)<norminv(1-0.05))=0;
            
        end
        
        function obj = MovePlotFreqs(obj,~,eventDat)
            %zpracovava stlaceni klavesy pro graf PlotUnepoched
            switch eventDat.Key
                case 'rightarrow' % +1 time window
                    obj.plotData.iTime = obj.plotData.iTime + 1; %min([obj.plotData.iTime + 1, size(obj.HFreqEpochs,4)]);
                case 'leftarrow'  % -1 time window
                    obj.plotData.iTime = max([obj.plotData.iTime - 1, 0]);%max([obj.plotData.iEpoch - 1, 1]);
                otherwise  
                   display(['key pressed: ' eventDat.Key]); %vypise stlacenou klavesu
            end
            obj.plotUnepochedData();
        end
       
        function plotUnepochedData(obj)
            % Plot all unepoched transformed data - power/original transformed signal and phase
            % one subplot = one frequency
            % called from PlotUnepoched()
            assert(~isempty(obj.HOrigData),'the field HOrigData with unepoched data cant be empty');
            time = (obj.plotData.iTime*obj.E.fs+1):(obj.plotData.iTime*obj.E.fs+obj.plotData.timeDelay); % moving time window
            x = time./obj.E.fs; % current x axis in seconds
            % determine number of subplots based on selected channels or frequencies
            if obj.plotData.plotGroup; numSubplot = length(obj.plotData.iChannels); % plot electrode channels 
            else numSubplot = length(obj.plotData.iFreqs)+2; end                     % plot frequencies + original eeg data

            for i = 1:numSubplot-1
                
                %%%%%%%%%%%%%%%%%%% ORIGINAL EEG DATA %%%%%%%%%%%%%%%%%%%
                if ~obj.plotData.plotGroup && i > numSubplot-2 % if plotting frequencies, display also original signal on the bottom
                    y = obj.HOrigData(time,obj.plotData.iChannels); % original eeg values
                    % LEFT
                    subplot(numSubplot,2,[i*2-1, i*2+1]);
                    plot(x, y, 'Color', 'black'); hold on;
                    xlim([x(1) x(end)]); % set x axis limits
                    ylim([-80,80]); % set y axis limits
                    % RIGHT
                    subplot(numSubplot,2,[i*2, i*2+2]);
                    plot(x, y, 'Color', 'black'); hold on;
                    xlim([x(1) x(end)]); % set x axis limits
                    ylim([-80,80]); % set y axis limits
                else
                    %%%%%%%%%%%%%%%%%%% POWER or ORIGINAL FILTERED data %%%%%%%%%%%%%%%%%%%
                    subplot(numSubplot,2,i*2-1)
                    if obj.plotData.plotGroup % data to be subplotted (given channel or frequency)
                        if obj.plotData.type %%%%%%%%%%%% ORIGINAL FILTERED DATA
                            y = squeeze(obj.E.freal(time,obj.plotData.iChannels(i),obj.plotData.iFreqs)); % channel original filtered values
                        else %%%%%%%%%%%%%%%%%%%%%%%%%%%% POWER DATA
                            y = squeeze(obj.E.HFreq(time,obj.plotData.iChannels(i),obj.plotData.iFreqs)); % channel power values
                        end
                        names = {obj.E.CH.H.channels.ass_brainAtlas}; %cast mozgu
                        labels = {obj.E.CH.H.channels.neurologyLabel}; %neurology label
                    else
                        if obj.plotData.type %%%%%%%%%%%% ORIGINAL FILTERED DATA
                            y = squeeze(obj.E.freal(time,obj.plotData.iChannels,obj.plotData.iFreqs(i))); % frequency original filtered values
                        else %%%%%%%%%%%%%%%%%%%%%%%%%%%% POWER DATA
                            y = squeeze(obj.E.HFreq(time,obj.plotData.iChannels,obj.plotData.iFreqs(i))); % channel power values
                        end
                        
                    end
                    if ~obj.plotData.type
                        neg = y<0; % negative frequency power values
                        plot(x(~neg),y(~neg),'b.',x(neg),y(neg),'r.','markers',4); hold on; % plots positive values with blue, negative with red
                    else
                        plot(x,y,'b');
                    end
                    xlim([x(1) x(end)]); % set x axis limits
                    ylim(obj.plotData.ylimits); % set y axis limits

                    iStimuli = find(obj.plotData.stimuli >= x(1) & obj.plotData.stimuli <= x(end));
                    t = '';
                    for iPodnet = 1:numel(iStimuli) % plot colored stimuli and responses
                        podnet = obj.E.PsyData.P.data(iStimuli(iPodnet),obj.E.PsyData.P.sloupce.kategorie);
                        plot([obj.plotData.stimuli(iStimuli(iPodnet)) obj.plotData.stimuli(iStimuli(iPodnet))]', ylim', 'LineWidth',3,'Color', obj.plotData.colors{podnet+1}); hold on;    
                        t = strcat(t,'{\color{', obj.plotData.colors{podnet+1},'}',obj.E.PsyData.P.strings.podminka(podnet+1),'}, '); % title
                    end

                    responses = obj.plotData.responses(obj.plotData.responses >= x(1) & obj.plotData.responses <= x(end));
                    plot([responses responses]', repmat(ylim,length(responses),1)','k', 'LineWidth',2); % plot responses
                    if obj.plotData.plotGroup
                        title(strcat(num2str(obj.plotData.iChannels(i)), '. channel (', t, ') ',names{obj.plotData.iChannels(i)},'-',labels{obj.plotData.iChannels(i)}));
                    else
                        title(strcat(num2str(obj.E.Hf(obj.plotData.iFreqs(i))), ' Hz (', t, ')'));
                    end

                    %%%%%%%%%%%%%%%%%%% PHASE %%%%%%%%%%%%%%%%%%% 
                    subplot(numSubplot,2,i*2)
                    % phase values
                    if obj.plotData.plotGroup
                        y = squeeze(obj.E.fphase(time,obj.plotData.iChannels(i),obj.plotData.iFreqs)); % channel phase values
                        disp('correct, more')
                    else
                        y = squeeze(obj.E.fphase(time,obj.plotData.iChannels,obj.plotData.iFreqs(i))); % frequency phase values
                    end 
                    %phasemap('rad')

                    %imagesc([x(1) x(end)],[0,5],repmat(y',[5,1])); hold on;

                    plot(x,y, 'Color', 'b'); hold on;
                    xlim([x(1) x(end)]); % set x axis limits
                    ylim([-5,5]); % set y axis limits
                    iStimuli = find(obj.plotData.stimuli >= x(1) & obj.plotData.stimuli <= x(end));
                    t = '';
                    for iPodnet = 1:numel(iStimuli) % plot colored stimuli and responses
                        podnet = obj.E.PsyData.P.data(iStimuli(iPodnet),obj.E.PsyData.P.sloupce.kategorie);
                        plot([obj.plotData.stimuli(iStimuli(iPodnet)) obj.plotData.stimuli(iStimuli(iPodnet))]', ylim', 'LineWidth',3,'Color', obj.plotData.colors{podnet+1}); hold on;    
                        t = strcat(t,'{\color{', obj.plotData.colors{podnet+1},'}',obj.E.PsyData.P.strings.podminka(podnet+1),'}, '); % title
                    end

                    responses = obj.plotData.responses(obj.plotData.responses >= x(1) & obj.plotData.responses <= x(end));
                    plot([responses responses]', repmat(ylim,length(responses),1)','k', 'LineWidth',2); % plot responses
                end
                
            end
           
        end
        
        function correct = PlotRejected(obj, channel, time, epoch, sum_chyby)
           %vraci true, pokud epocha neni oznacena jako vyrazena
           %pouziva se v PlotAllEpochs
           correct = false;           
           if obj.E.PsyData.P.data(epoch, obj.E.PsyData.P.sloupce.zpetnavazba) % mark trening answers with red
               hold on; plot([time(1) time(end)], [obj.E.Hf(end) obj.E.Hf(1)],'red','LineWidth',6)
           
           elseif ~obj.E.PsyData.P.data(epoch, 3) || sum_chyby > 0 % mark wrong answers with black
               hold on; plot([time(1) time(end)], [obj.E.Hf(end) obj.E.Hf(1)],'black','LineWidth',6)
           
           elseif ismember(epoch, obj.E.RjEpoch) % mark rejected epochs with green
               hold on; plot([time(1) time(end)], [obj.E.Hf(end) obj.E.Hf(1)],'green','LineWidth',6)
           
           elseif obj.E.RjEpochCh(channel, epoch) % mark rejected channel's epochs with blue
               hold on; plot([time(1) time(end)], [obj.E.Hf(end) obj.E.Hf(1)],'blue','LineWidth',6)
           else
              correct = true;
           end
        end
        
        function iEp = CorrectEpochs(obj, channel, iconditions)
            %pouzite v PlotFrequencyPower
            %vrati indexy vsetkych spravnych epoch pre dany channel a condition
            %napr. aedist condition 0=cervena, 1=vy, 2=znacka, 9=all

            iEp = obj.E.GetEpochsExclude();
            iEp = find(squeeze(iEp(channel,:))); %indexy nenulovych prvkov na vyradenie

            if iconditions(1,1) ~= 9 %ak chceme len jednu konkretnu condition
                iEp = intersect(iEp,find(sum(obj.E.PsyData.P.data(:,obj.E.PsyData.P.sloupce.kategorie)' == iconditions',1)));
            end

            iEp = setdiff(iEp, find(obj.E.RjEpochCh(channel,:)));
        end
        
        function [fmean, fstd] = meanZscoredPower(obj, time, channel, epochs)
            %vrati priemernu power a strednu chybu priemeru pre kazdu
            %frekvenciu cez cas a epochy
            mean_time = squeeze(mean(obj.E.HFreqEpochs(time,channel,:,epochs),1)); %priemerna power pre kazdu frekvenciu cez cas -> ostane freq x epoch
            fmean = mean(mean_time,2); %priemerna power pre kazdu frekvenciu cez epochy
            fstd = std(mean_time,0,2)./sqrt(length(epochs)); %stredna chyba priemeru pre kazdu frekvenciu cez epochy -> std()/sqrt(n)
        end
        
        function subplotTimeFrequency(obj, data, time)
            colormap parula; 
            imagesc(data, 'XData', time, 'YData', obj.E.Hf); 
            set(gca,'YDir','normal');
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            h = colorbar;
            ylabel(h, 'z-scored power'); 
           % title([names{channel}, ' ', labels{channel}]);
        end
        
        function plotEpochData(obj)
            %vykresli vlavo time x frequency graf pre danu epochu 
            %vpravo average power cez vsetky frekvencie danej epochy 
            %pouziva sa v PlotMovingEpochs
            assert(~isempty(obj.E.HFreqEpochs),'soubor s frekvencnimi daty pro epochy neexistuje');
            subplot(1,2,1) % subplot time x frequency power for given epoch
            imagesc(squeeze(obj.E.HFreqEpochs(:,obj.plotData.channels(obj.plotData.iChannel),:,obj.plotData.iEpoch))', 'XData', obj.plotData.T, 'YData', obj.E.Hf);
            caxis(obj.plotData.zlimits(obj.plotData.iChannel,:));xlabel('Time (s)'); ylabel('Frequency (Hz)'); 
            colormap parula; %aby to bylo jasne u vsech verzi matlabu - i 2016
            set(gca,'YDir','normal');
            colorbar;
            
            hold on; % plot rejected line 
            obj.PlotRejected(obj.plotData.channels(obj.plotData.iChannel), obj.plotData.T, obj.plotData.iEpoch, sum(obj.plotData.rejectedEpochs(obj.plotData.iEpoch,:)));
            response_time = obj.E.PsyData.P.data(obj.plotData.iEpoch, 4);
            hold on; % plot response time 
            plot([response_time response_time], [obj.E.Hf(1)-10 obj.E.Hf(end)+10],'black','LineWidth',4);
            hold on;
            
            subplot(1,2,2) % subplot mean power across all frequencies
            plot(obj.plotData.T, obj.E.d(:,obj.plotData.channels(obj.plotData.iChannel), obj.plotData.iEpoch)');
            
            hold on; % plot response time 
            plot([response_time response_time], obj.plotData.zlimits(obj.plotData.iChannel,:), 'black', 'LineWidth', 4);
            ylim(obj.plotData.zlimits(obj.plotData.iChannel,:)); % y axis = zlimits (default/specified by user)
            xlim([obj.plotData.T(1) obj.plotData.T(end)]); % x axis = time
            xlabel('Time (s)'); ylabel('Power');
            title(sprintf('%s - channel %d epoch %d', obj.E.epochData{obj.plotData.iEpoch,1}, obj.plotData.channels(obj.plotData.iChannel), obj.plotData.iEpoch), 'FontSize', 12);
            hold off;
        end
        
        function zlimits = getZlimits(obj, ch)
            %vypocita minimalnu a maximalnu power pre dany channel cez
            %vsetky epochy a frekvencie
            %pouziva sa v PlotMovingEpochs
            ymin = min(min(obj.E.d(:, ch, :)));
            ymax = max(max(obj.E.d(:, ch, :)));
            zlimits = [ymin ymax];
        end
        
        function obj = MovePlotEpochs(obj,~,eventDat)
            %zpracovava stlaceni klavesy pro graf PlotMovingEpochs
            switch eventDat.Key
                case 'rightarrow' % +1 epoch
                    obj.plotData.iEpoch = min([obj.plotData.iEpoch + 1, size(obj.E.HFreqEpochs,4)]);
                case 'leftarrow'  % -1 epoch
                    obj.plotData.iEpoch = max([obj.plotData.iEpoch - 1, 1]);
                case 'uparrow'    % -1 channel
                    obj.plotData.iChannel = max([obj.plotData.iChannel - 1, 1]);
                case 'downarrow'  % +1 channel
                    obj.plotData.iChannel = min([obj.plotData.iChannel + 1, length(obj.plotData.channels)]);
                case 'numpad6' % nasledujuca s rovnakou condition
                    obj.plotData.iEpoch = min([obj.getNextCondition(1), size(obj.E.HFreqEpochs,4)]);
                case 'numpad4' % predchadzajuca s rovnakou condition
                    obj.plotData.iEpoch = max([obj.getLastCondition(1), 1]);
                case 'pageup' % nasledujuca condition
                    obj.plotData.iEpoch = min([obj.getNextCondition(0), size(obj.E.HFreqEpochs,4)]);
                case 'pagedown' % predchadzajuca condition
                    obj.plotData.iEpoch = max([obj.getLastCondition(0), 1]);
                case {'multiply','8'} %hvezdicka na numericke klavesnici
                   %dialog na vlozeni minima a maxima osy y
                   answ = inputdlg('Enter ymax and min:','Yaxis limits', [1 50], {num2str(obj.plotData.zlimits(obj.plotData.iChannel,:))});
                   if numel(answ)>0  %odpoved je vzdy cell 1x1 - pri cancel je to cell 0x0
                       if isempty(answ{1}) || any(answ{1}=='*') %pokud vlozim hvezdicku nebo nic, chci znovy spocitat max a min
                           obj.plotData.zlimits(obj.plotData.iChannel,:) = obj.getZlimits(obj.plotData.channels(obj.plotData.iChannel));
                       else %jinak predpokladam dve hodnoty
                           data = str2num(answ{:});  %#ok<ST2NM>
                           if numel(data)>= 2 %pokud nejsou dve hodnoty, nedelam nic
                             obj.plotData.zlimits(obj.plotData.iChannel,:) = [data(1) data(2)];
                           end
                       end
                   end
                   obj.plotEpochData(); %prekreslim grafy
                otherwise  
                   display(['key pressed: ' eventDat.Key]); %vypise stlacenou klavesu
            end
            obj.plotEpochData();
        end
        
        function last = getLastCondition(obj, same)
            %najde index poslednej najblizsej epochy 
            %rovnakej(same=1)/rozdielnej(same=0) kategorie
            %pouziva sa v PlotMovingEpochs pri numpad4/pagedown
            if same
                last = find(obj.E.PsyData.P.data(1:(obj.plotData.iEpoch-1),7) == obj.E.epochData{obj.plotData.iEpoch,2}, 1, 'last');
            else 
                last = find(obj.E.PsyData.P.data(1:(obj.plotData.iEpoch-1),7) ~= obj.E.epochData{obj.plotData.iEpoch,2}, 1, 'last');
            end  
            if isempty(last)
                last = obj.plotData.iEpoch;
            end
        end
        
        function next = getNextCondition(obj, same)
            %najde index najblizsej epochy 
            %rovnakej(same=1)/rozdielnej(same=0) kategorie
            %pouziva sa v PlotMovingEpochs pri numpad6/pageup
            if same
                next = obj.plotData.iEpoch + find(obj.E.PsyData.P.data((obj.plotData.iEpoch+1):end,obj.E.PsyData.P.sloupce.kategorie) == obj.E.epochData{obj.plotData.iEpoch,2}, 1);
            else 
                next = obj.plotData.iEpoch + find(obj.E.PsyData.P.data((obj.plotData.iEpoch+1):end,obj.E.PsyData.P.sloupce.kategorie) ~= obj.E.epochData{obj.plotData.iEpoch,2}, 1);
            end
            if isempty(next)
                next = obj.plotData.iEpoch;
            end
        end
        
        function obj = hybejPlotISPC(obj,~,eventDat)  
            %reacts to events in figure ISPCPlot
            if ~isempty(eventDat.Modifier) && strcmp(eventDat.Modifier{:},'shift')
                shift = 10;
            else
                shift = 1;
            end
            switch eventDat.Key
                case {'numpad6','d'} %next time sample
                    if obj.plotISPC.TimeHf(1) <= size(obj.plotISPC.ispc,3)-shift %compare to num of time samples
                        TimeHf = [(obj.plotISPC.TimeHf(1) + shift), obj.plotISPC.TimeHf(2)];
                        obj.ISPCPlot(obj.plotISPC.chnpair, TimeHf);
                    end
                case {'numpad4','a'} %previous time sample   
                    if obj.plotISPC.TimeHf(1) > shift
                        TimeHf = [(obj.plotISPC.TimeHf(1) - shift), obj.plotISPC.TimeHf(2)];
                        obj.ISPCPlot(obj.plotISPC.chnpair, TimeHf)
                    end         
                case {'numpad8','w'} %next frequency
                    if obj.plotISPC.TimeHf(2) <= size(obj.plotISPC.ispc,4) - shift   %compare to num of freq 
                        TimeHf = [obj.plotISPC.TimeHf(1) , (obj.plotISPC.TimeHf(2)+shift)];
                        obj.ISPCPlot(obj.plotISPC.chnpair, TimeHf)
                    end
                case {'numpad2','s'} %previous frequency 
                    if obj.plotISPC.TimeHf(2) > shift
                        TimeHf = [obj.plotISPC.TimeHf(1) , (obj.plotISPC.TimeHf(2)-shift)];
                        obj.ISPCPlot(obj.plotISPC.chnpair, TimeHf)
                    end                
                case {'home'}                     
                    obj.ISPCPlot(obj.plotISPC.chnpair, []);
                case {'return'}
                    if obj.plotISPC.onlysub2 == 1
                       obj.plotISPC.onlysub2 = 0;
                       obj.plotISPC.clf = 1;
                    end
                    obj.ISPCPlot(obj.plotISPC.chnpair, [], [1 1 1 1]);
                case {'rightarrow'} %move the channel pair to right
                    if obj.plotISPC.chnpair(1) <= numel(obj.plotISPC.channels) -shift
                        chnpair = [obj.plotISPC.chnpair(1)+shift obj.plotISPC.chnpair(2)];
                        obj.ISPCPlot(chnpair, obj.plotISPC.TimeHf, [0 1 0 0]);
                    end
                case {'leftarrow'} % move the channel pair to left
                    if obj.plotISPC.chnpair(1) > shift
                        chnpair = [obj.plotISPC.chnpair(1)-shift obj.plotISPC.chnpair(2)];
                        obj.ISPCPlot(chnpair, obj.plotISPC.TimeHf, [0 1 0 0]);
                    end
                case {'uparrow'} % zvyseni cisla aktivni statistiky
                    if obj.plotISPC.chnpair(2) <= numel(obj.plotISPC.channels) -shift
                        chnpair = [obj.plotISPC.chnpair(1) obj.plotISPC.chnpair(2)+shift];
                        obj.ISPCPlot(chnpair, obj.plotISPC.TimeHf, [0 1 0 0]);
                    end
                case {'downarrow'}
                    if obj.plotISPC.chnpair(2) > shift
                        chnpair = [obj.plotISPC.chnpair(1) obj.plotISPC.chnpair(2)-shift];
                        obj.ISPCPlot(chnpair, obj.plotISPC.TimeHf, [0 1 0 0]);
                    end
                case 'q'
                    if obj.plotISPC.plotsig(1) == 0.5
                        obj.plotISPC.plotsig = [0.95 1];
                    elseif obj.plotISPC.plotsig(1) == 0.95 
                        obj.plotISPC.plotsig = [0.99 1];
                    else
                        obj.plotISPC.plotsig = [0.5 0.95];
                    end
                    obj.ISPCPlot(obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [0 1 0 1]);
                case {'o'}    %plots all epochs for this frequency - their phase
                    chnpair = [obj.plotISPC.channels(obj.plotISPC.sortorder(obj.plotISPC.chnpair(1))) ...
                            obj.plotISPC.channels(obj.plotISPC.sortorder(obj.plotISPC.chnpair(2)))];
                    obj.PlotITPCallEpochs(chnpair,obj.plotISPC.TimeHf(2));
                case {'p'}  %plot optional plot all phase angles for selected time and freq
                    obj.ISPCoPlot(obj.plotISPC.chnpair, obj.plotISPC.TimeHf,[],1)
                case 'space' %zoom/unzoom the top right plot
                    obj.plotISPC.onlysub2 = 1-obj.plotISPC.onlysub2;
                    if obj.plotISPC.onlysub2 == 0
                       obj.plotISPC.clf = 1;
                    end
                    obj.ISPCPlot(obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [0 1 0 0]); 
                case 'v' %change sorting of channels in the top right figure - original vs. brainlabels
                    if obj.plotISPC.sortby == 0
                        labels={obj.E.CH.brainlabels(:).label};
                        [slabels,sortindex]=sort(labels);
                        [~,ia] =unique(slabels,'stable');
                        obj.plotISPC.els = [ia(2:end)'-1 numel(obj.plotISPC.channels)];
                        obj.plotISPC.sortorder = sortindex;
                        obj.plotISPC.sortby = 1;
                    else
                        obj.plotISPC.els=obj.E.CH.els;
                        obj.plotISPC.sortorder = 1:size(obj.plotISPC.ispc,1);
                        obj.plotISPC.sortby = 0;
                    end
                    obj.ISPCPlot(obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [0 1 0 0]); %update the all channel plot
                case 'subtract'
                    if obj.plotISPC.aplha > 0
                        obj.plotISPC.aplha = obj.plotISPC.aplha - 0.1;
                        obj.ISPCPlot(obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [1 0 0 1]);
                    end
                case 'add'
                    if obj.plotISPC.aplha < 1
                        obj.plotISPC.aplha = obj.plotISPC.aplha + 0.1;
                        obj.ISPCPlot(obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [1 0 0 1]);
                    end
            end
        end
        function hybejPlotISPCClick(obj,~,~)     
            mousept = get(gca,'currentPoint');
            if mousept(1,1)<1
                %probably click in the time freq chart
                
                T = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1));
                F =  obj.E.Hfmean; 
                x = find(T>=mousept(1,1),1); 
                y = find(F>=mousept(1,2),1); %coordinated of the click
                if x <= size(obj.plotISPC.ispc,3) && x>=0 && y <= size(obj.plotISPC.ispc,4) && y >= 0%compare to num of time samples                        
                        obj.ISPCPlot(obj.plotISPC.chnpair, [x y],[1 0 1 1]);
                end                
            else %click in the channel plot
                x = round(mousept(1,1)); y = round(mousept(1,2)); %coordinated of the click
                if x <= numel(obj.plotISPC.channels) && x>=1 && y<=numel(obj.plotISPC.channels) && y >=1                        
                            updateS = iff(obj.plotISPC.onlysub2 == 0,[1 1 1 1],[0 1 0 0]);
                            obj.ISPCPlot([x y], obj.plotISPC.TimeHf,updateS);
                end
            end
            
        end
         
    end
    
end

