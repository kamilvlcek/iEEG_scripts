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
        
        %ISPC functions (for computing PLV between two channels)
        function [ispc] = ISPCCalculate(obj,channels,fdr)
            % inter site phase clustering = PLV - according to MikeXCohen
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
                fprintf('\b\b\b\b\b%5i',ichn1); %delete previous three characters and print the %i in three characters
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
        
        function ISPCFilterChanPair(obj, byROI) % Sofia from 6.6.2022 
            % filter pairs of chan which have ISPC values significant for each condition separately or by ROI - to leave only between ROI pairs of chan (not within)
            % byROI = 0 to filter pairs of chan with significant ISPC values; 1 to filter pairs of chan between ROIs   
             
            if isfield(obj.plotISPC,'channels') && ~isempty(obj.plotISPC.channels)
                channels = obj.plotISPC.channels;
            else
                channels = 1:obj.E.channels; % defaul - all channels
            end
            
            if byROI == 0 % filter pairs of chan with ISPC signif for each condition separately
                n_kats = numel(obj.plotISPC.ispc_sign_chanPair); % total number of conditions with computed ispc
                chanPair_sign_filtered = cell(n_kats,1);
                for ik = 1:n_kats                    
                    [ch1, ch2] = find(obj.plotISPC.ispc_sign_chanPair{ik} == 1); % find pair of chan with significant ispc after cluster correction
                    chanPair = [ch1, ch2];
                    chanPair_sign_filtered{ik} = chanPair;                    
                end
                obj.plotISPC.chanPair_sign_filtered = chanPair_sign_filtered;
                disp('channel pairs with significant ISPC were filtered' );
                
            else  % filter only between ROI pairs of chan
                chanPair = ones(length(channels)); % matrix of all possible pairs
                for ichn1 = 1:length(channels)
                    for ichn2 = 1:length(channels)
                        if strcmp(obj.E.CH.brainlabels(channels(ichn1)).label,obj.E.CH.brainlabels(channels(ichn2)).label)
                             chanPair(ichn1, ichn2) = 0; % if brain labels of that chan pair are the same, filter that pair out
                        end
                    end
                end
                [ch1, ch2] = find(triu(chanPair) == 1); % find pairs of chan with 1 (different brain labels), without repetition 
                chanPair = [ch1, ch2];   
                obj.plotISPC.roi_chanpair = chanPair;
                disp('between ROI pairs of channels were filtered' );
            end
        end
        
        function [ispc_cats, ispc_cats_signif, ispc_cats_clust_corr] = ISPCCalculateCats(obj, kats, channels, ChanPair_betweenROI, exclude_ChanPair_lowEpochs, baseline)  % Sofia from 6.6.2022 
            % inter site phase clustering = phase-locking value (PLV) is computed for each condition with permutation statistics; incorrect and epi epochs are excluded
            %%% input: kats - indexes of conditions, e.g. [2,3]
            % channels - channels, for which we want to compute PLV, e.g. [1:30]
            % ChanPair_betweenROI = 1 to compute ispc only for between ROI pairs of chan, 0 - for all chan pairs (default) 
            % exclude_ChanPair_lowEpochs = 1; not to compute ispc for chan pairs with number of epochs < 30; 0 - to compute
            % baseline - baseline period to substract to normalize ispc values
            %%% output: ispc_cats = {kats}(chn x chn x time x fq) plv without statistics
            % ispc_cats_signif = {kats}(chn x chn x time x fq) plv after permutation statistics (uncorrected), with values below permuted threshold set to 0
            % ispc_cats_clust_corr = {kats}(chn x chn x time x fq) plv with significant values after cluster correction
            
            tic
            if ~exist('channels', 'var') || isempty(channels)
                channels = 1:obj.E.channels; % default - all channels
            end
            if ~exist('ChanPair_betweenROI', 'var') || isempty(ChanPair_betweenROI)
                ChanPair_betweenROI = 0; % default all chan pairs
            end                    % if 1 computes only for between ROI pairs of chan, but first call function ISPCFilterChanPair
            
            if ~exist('exclude_ChanPair_lowEpochs', 'var') || isempty(exclude_ChanPair_lowEpochs)
                exclude_ChanPair_lowEpochs = 1; % default - exclude chan pairs with number of epochs < 30 (not compute ispc = NaN), according to Mike X Cohen’s book p.340-341, ISPC-trials is sensitive to the number of trials, with lower ustable results
            end 
            
            if ~exist('baseline', 'var') || isempty(baseline)
                baseline = obj.E.baseline; % default - the same baseline as used in time-frequency analysis           
            end
                      
            % define categories = conditions
            if ~exist('kats', 'var') || isempty(kats)
                if isfield(obj.E.plotRCh,'kategories') && ~isempty(obj.E.plotRCh.kategories)
                    kats = obj.E.plotRCh.kategories;
                else
                    kats = obj.E.Wp(obj.E.WpActive).kats;
                end
            end
            obj.plotISPC.ispc_cats_names{1} = obj.E.PsyData.CategoryName(kats, []); % save names of conditions for which ISPC was computed
            obj.plotISPC.ispc_cats_names{2} = kats; % but also their original numbers - needed for ploting then
            ispc_cats = cell(numel(kats),1); %{kats}(chn x chn x time x fq)
            ispc_cats_signif = cell(numel(kats),1); %{kats}(chn x chn x time x fq)
            ispc_cats_clust_corr = cell(numel(kats),1); %{kats}(chn x chn x time x fq)
            ispc_cats_n_epochs = cell(numel(kats),1); %{kats}(chn x chn) number of epochs for each pair of chan, over which ispc was computed
            ispc_sign_chan = cell(numel(kats),1); %{kats}(chn x chn)  % according to the old code of Kamil - ispc_p_min         
               
            for i = 1:numel(kats) % initialize all matrices inside the cells
                ispc_cats{i} = NaN(numel(channels), numel(channels), size(obj.E.fphaseEpochs, 1), size(obj.E.fphaseEpochs, 3)); % NaN will be for non computed chan pairs
                ispc_cats_signif{i} = NaN(numel(channels), numel(channels), size(obj.E.fphaseEpochs, 1), size(obj.E.fphaseEpochs, 3));
                ispc_cats_clust_corr{i} = NaN(numel(channels), numel(channels), size(obj.E.fphaseEpochs, 1), size(obj.E.fphaseEpochs, 3));                
                ispc_cats_n_epochs{i} = NaN(numel(channels), numel(channels));
                ispc_sign_chan{i} = zeros(numel(channels), numel(channels)); % to distinguish non-significant and not computed values
            end
            
            fprintf('ISPCCalculateCats:    ');
            for ichn1 = 1:numel(channels)
                if sum(obj.E.CH.RjCh==channels(ichn1))>0, continue; end % if the channel is rejected, not compute
                fprintf('\b\b\b\b\b%5i',ichn1); % delete previous three characters and print the %i in three characters
                
                for ichn2 = ichn1:numel(channels) % iterate only through the upper triangular part of the matrix, not to compute the same value twice for each chan pair
                    if sum(obj.E.CH.RjCh==channels(ichn2))>0, continue; end
                    if ~isa(obj.E, 'CHilbertMulti') || find(obj.E.els >= channels(ichn1),1) == find(obj.E.els >= channels(ichn2),1) %find expression to get the pacient number 
                        for k = 1:numel(kats) % over all conditions
                            if ChanPair_betweenROI==1 && ~(ismember([channels(ichn1), channels(ichn2)], obj.plotISPC.roi_chanpair,'rows')|| ismember([channels(ichn2), channels(ichn1)], obj.plotISPC.roi_chanpair,'rows')) % if the current pair of chan is not between ROI pair
                                ispc_sign_chan{k}(ichn1,ichn2) = NaN;
                                ispc_sign_chan{k}(ichn2,ichn1) = NaN; % symmetrical position
                            else
                                % correct epochs for chan 1
                                iEp1 = obj.CorrectEpochs(channels(ichn1), kats(k)); % only correct epochs for each condition (without errors, training and epi epochs)
                                % correct epochs for chan 2
                                iEp2 = obj.CorrectEpochs(channels(ichn2), kats(k));
                                % for computing ISPC, the number of epochs should be the same in both channels
                                iEpBoth = intersect(iEp1, iEp2);
                                if exclude_ChanPair_lowEpochs && numel(iEpBoth) < 30 % not to compute ispc for chan pairs with a low number of epochs
                                    ispc_sign_chan{k}(ichn1,ichn2) = NaN;
                                    ispc_sign_chan{k}(ichn2,ichn1) = NaN; % symmetrical position
                                    continue
                                end
                                fphaseEpochsCh1 = squeeze(obj.E.fphaseEpochs(:,channels(ichn1),:,iEpBoth)); % time x freq x epochs in ch1
                                fphaseEpochsCh2 = squeeze(obj.E.fphaseEpochs(:,channels(ichn2),:,iEpBoth)); % time x freq x epochs in ch2
                                
                                % compute ispc - over all frequencies and time at once
                                ispc_all = fphaseEpochsCh1 - fphaseEpochsCh2; % time x freq x epochs - just phase angle difference between 2 chan (we need them for permut stat)                               
                                ispc_ch = abs(mean(exp(1i*ispc_all),3)); % time x freq, computation of ispc over epochs                                                                              
                                ispc_cats_n_epochs{k}(ichn1,ichn2) = numel(iEpBoth); % number of epochs for each pair of chan
                                
                                % compute statistics
                                [ispc_norm, ispc_signif, ispc_clust_corr] = obj.ISPCPermut(ispc_ch, ispc_all, baseline); % permutation statistics with cluster correction
                                ispc_cats{k}(ichn1,ichn2,:,:) = ispc_norm; % put normalized ispc after baseline substraction in cell array with all conditions
                                ispc_cats_signif{k}(ichn1,ichn2,:,:) = ispc_signif; % put thresholded ispc values in cell array with all conditions
                                ispc_cats_clust_corr{k}(ichn1,ichn2,:,:) = ispc_clust_corr; % cluster corrected significant ispc values
                                if sum(sum(ispc_clust_corr))>0                              
                                    ispc_sign_chan{k}(ichn1,ichn2) = 1; % assign 1 to pair of chan if their significant ispc values survived after cluster correction
                                end
                                
                                % assign the computed values to the symmetrical positions
                                ispc_cats{k}(ichn2,ichn1,:,:) = ispc_norm;
                                ispc_cats_signif{k}(ichn2,ichn1,:,:) = ispc_signif;
                                ispc_cats_clust_corr{k}(ichn2,ichn1,:,:) = ispc_clust_corr;
                                ispc_sign_chan{k}(ichn2,ichn1) = ispc_sign_chan{k}(ichn1,ichn2);
                                ispc_cats_n_epochs{k}(ichn2,ichn1) = ispc_cats_n_epochs{k}(ichn1,ichn2);
                            end
                        end
                    end
                end
            end
            obj.plotISPC.ispc_cats = ispc_cats;
            obj.plotISPC.ispc_cats_n_epochs = ispc_cats_n_epochs;
            obj.plotISPC.ispc_cats_signif = ispc_cats_signif;
            obj.plotISPC.ispc_cats_clust_corr = ispc_cats_clust_corr;
            obj.plotISPC.channels = channels;            
            obj.plotISPC.ispc_sign_chanPair = ispc_sign_chan;
            toc
        end
        
        function [ispc_norm, ispc_signif, ispc_clust_corr] = ISPCPermut(obj, ispc_computed, ispc_all_trials, baseline, n_permutes, threshold) % Sofia from 18.10.2023
            % computes permutation-based statistics for PLV values (adapted code from Mike's Cohen book, figure 34.9) for one pair of channels
            %%% Mike's Cohen book, p. 489: the null hypothesis is that shifting the time course of the phase angle time series by a random amount would not affect the final ISPC strength. 
            %%% The null hypothesis can be generated by shifting the time series of each trial by a random amount         
            % ispc_computed  - already computed PLV (ISPC) data for one pair of chan: time x freq 
            % ispc_all_trials - phase differences from all trials: time x freq x trials
            % baseline - baseline period to substract to normalize ispc values
            % n_permutes - number of permutations (default = 200)
            % threshold - to consider computed PLV significant (default = 0.05)
            % output: ispc_norm - normalized ispc values after baseline substraction
            % ispc_signif - thresholded PLV (ISPC) data for one pair of chan: time x freq with values below the threshold set to zero
            % ispc_clust_corr - time x freq, significant PLV (ISPC) after the cluster correction (method to deal with multiple comparisons in permutat stat) 
            
            if ~exist('n_permutes', 'var') || isempty(n_permutes)
                n_permutes = 200; % default, if more e.g. 500, it is becoming very time-consuming
            end
            
            if ~exist('threshold', 'var') || isempty(threshold)
                threshold = 0.05; % default
            end
            
            iepochtime = round(obj.E.epochtime(1:2).*obj.E.fs); % epochtime in samples
            ibaseline =  round(baseline.*obj.E.fs); % baseline in samples 
            time_interval = abs(iepochtime(1)-ibaseline(2))+1 : size(ispc_computed,1); % time interval in samples over which compute stats (after stimulus period)                                
            ispc_all_trials_afterStim = ispc_all_trials(time_interval, :, :); % take the data only after stimulus
            ispc_computed_afterStim = ispc_computed(time_interval, :);
            
            % initialize null hypothesis matrices
            eegtemp = zeros(size(ispc_all_trials_afterStim)); % time x freq x trials
            permuted_vals = zeros(n_permutes, size(ispc_all_trials_afterStim,1), size(ispc_all_trials_afterStim,2)); % n_permutes x time x freq
            time_points = size(ispc_all_trials_afterStim,1);
            
            for permi = 1:n_permutes
                for triali = 1:size(ispc_all_trials_afterStim, 3)                   
                    % shifting the time series of each trial by a random amount
                    cutpoint = randsample(2:time_points-2, 1);
                    eegtemp(:,:,triali) = ispc_all_trials_afterStim([cutpoint:end 1:cutpoint-1], :, triali);
                end
                permuted_vals(permi,:,:) = squeeze(abs(mean(exp(1i*eegtemp),3))); % mean over epochs
%                 cutpoint = randsample(2:time_points-2,1); % to replace the loop above
%                 permuted_vals(permi,:,:) = squeeze(abs(mean(exp(1i*ispc_all_trials([cutpoint:end 1:cutpoint-1],:,:)),3))); % permuted plv
            end
            
            % transform to z values and set to zero values below thershold
            zmap = (ispc_computed_afterStim-squeeze(mean(permuted_vals,1))) ./ squeeze(std(permuted_vals,1));
            ispc_signif = ispc_computed_afterStim;
            ispc_signif(abs(zmap)<norminv(1-threshold))=0;                       
           
            % the cluster correction on the permuted data (to deal with multiple comparisons, adapted code from Mike's Cohen book, figure 34.9)
            max_clust_info = zeros(n_permutes,1);
            for permi = 1:n_permutes               
                % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
                fakecorrsz = squeeze((permuted_vals(permi,:,:)-mean(permuted_vals,1)) ./ std(permuted_vals,[],1) );
                fakecorrsz(abs(fakecorrsz)<norminv(1-threshold))=0;
                
                % get number of elements in largest supra-threshold cluster
                clustinfo = bwconncomp(fakecorrsz);
                max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
                % using cellfun here eliminates the need for a slower loop over cells
            end
            
            % apply cluster-level corrected threshold
            ispc_clust_corr = ispc_signif;            
            % find islands and remove those smaller than cluster size threshold
            clustinfo = bwconncomp(ispc_clust_corr);
            clust_info = cellfun(@numel,clustinfo.PixelIdxList);
            clust_threshold = prctile(max_clust_info,100-threshold*100);
            
            % identify clusters to remove
            whichclusters2remove = find(clust_info<clust_threshold);
            
            % remove clusters
            for i=1:length(whichclusters2remove)
                ispc_clust_corr(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
            end
            
            % add zeros in baseline period for both matrices with significant values to make them the
            % same size as computed ispc (for easier ploting in the future)
            ispc_signif = padarray(ispc_signif, [time_interval(1)-1 0], 0, 'pre');
            ispc_clust_corr = padarray(ispc_clust_corr, [time_interval(1)-1 0], 0, 'pre');
            
            % baseline substraction - normalization of ispc values
            ispc_baseline = mean(ispc_computed(ibaseline(1)-iepochtime(1)+1 : ibaseline(2)-iepochtime(1),:),1); % mean over time (baseline period)
            ispc_norm = bsxfun(@minus,ispc_computed,ispc_baseline); % time x freq, substract mean baseline
            % replace the values which are significant by normalized ispc values
            ispc_signif(ispc_signif ~= 0) = ispc_norm(ispc_signif ~= 0);
            ispc_clust_corr(ispc_clust_corr~=0) = ispc_norm(ispc_clust_corr~=0);
            
        end
        
        function [bins_n_chanPair, max_bins] = ISPCFindSignBins(obj, freq_bin, time_bin, how_aggregate, n_max_bins) % Sofia since 30.10.2023
            % aggregates significant ispc data to enlarged bins with the spicified size of bins: freq_bin - for freq in Hz, time_bin - for time in sec,
            % either by using mean of smaller bins (how_aggregate = 0) or max value (how_aggregate = 1)
            % also counts the number of channel pairs that have significant values in each enlarged bin (output - bins_n_chanPair)
            % and returns the time window and frequency range in which most chan pairs show the significance (max_bins = {kats}[itime ifreq time freq n_chan_pairs])
            if ~exist('freq_bin', 'var') || isempty(freq_bin)
                freq_bin = 1; % default, freq bin = 1 Hz       
            end
            
            if ~exist('time_bin', 'var') || isempty(time_bin)
                time_bin = 0.5; % default, time bin = 0.5 sec         
            end
            
            if ~exist('how_aggregate', 'var') || isempty(how_aggregate)
                how_aggregate = 0; % how we want aggregate data in each bin: default=0 - by taking mean of smaller bins; 1 - by taking the maximum value
            end
            
            if ~exist('n_max_bins', 'var') || isempty(n_max_bins)
                n_max_bins = 5; % how many bins we want to find with the largest number of significant chan pairs
            end
            
            if ~isfield(obj.plotISPC,'PlotISPCbins')
                obj.plotISPC.PlotISPCbins = struct; % create a separate struct to save all binned data here
            end

            ispc_cats_clust_corr = obj.plotISPC.ispc_cats_clust_corr; % {kats}(chn x chn x time x fq),  significant ispc values after cluster correction for each condition            
            nkats = size(ispc_cats_clust_corr, 1); % number of conditions
            nchannels = size(ispc_cats_clust_corr{1},1); % number of channels
            
            freq = obj.E.Hf(1) : freq_bin : obj.E.Hf(end); 
            obj.plotISPC.PlotISPCbins.freq = freq; % new enlarged frequencies bins
            obj.plotISPC.PlotISPCbins.time = obj.E.baseline(2):time_bin:obj.E.epochtime(2); % new enlarged time bins           
            [~, ifreq] = ismember(freq, obj.E.Hf); % indexes of enlarged frequencies bins        
            iepochtime = round(obj.E.epochtime(1:2).*obj.E.fs); % epochtime in samples
            ibaseline =  round(obj.E.baseline.*obj.E.fs); % baseline in samples
            itime_bin = round(time_bin*obj.E.fs); % size of time bin in samples
            itime = abs(iepochtime(1)-ibaseline(2))+1 : itime_bin : size(ispc_cats_clust_corr{1},3); % indexes of enlarged time bins 
            itime(end+1) = size(ispc_cats_clust_corr{1},3); % because the lenght of epoch is not even number (memact, delay = 3.9 sec)
            
            ispc_cats_bins = cell(nkats,1); %{kats}(chn1 x chn2 x time bins x freq bins) new binned ispc data
                        
            % first get average or max ispc for each new bin, each channel pair and each condition
            for ichn1 = 1:nchannels
                for ichn2 = 1:nchannels
                    for k = 1:nkats
                        for ifr = 1:numel(ifreq)-1
                            for it = 1:numel(itime)-1
                                if ~isnan(obj.plotISPC.ispc_sign_chanPair{k}(ichn1, ichn2)) % if this chan pair has computed ispc
                                    if how_aggregate == 0
                                        temp = mean(mean(squeeze(ispc_cats_clust_corr{k}(ichn1,ichn2,itime(it):itime(it+1),ifreq(ifr):ifreq(ifr+1))),1)); % mean over selected frequency and time bins
                                    else
                                        temp = max(max(squeeze(ispc_cats_clust_corr{k}(ichn1,ichn2,itime(it):itime(it+1),ifreq(ifr):ifreq(ifr+1))),[],1)); % max value from selected frequency and time bins                                       
                                    end
                                    ispc_cats_bins{k}(ichn1, ichn2,it,ifr) = temp; % one average or max ispc value for this new time-freq bin for this channel pair
                                else
                                    ispc_cats_bins{k}(ichn1, ichn2,it,ifr) = NaN; % for not computed channel pairs, leave NaN
                                end
                            end
                        end
                    end
                end
            end
            
            % then count the number of chan pairs with non-zero values in each bin
            bins_n_chanPair = cell(nkats, 1); % {kats} time x freq, store a number of channels for each enlarged bin
            bins_idx_chanPairs = cell(nkats, numel(itime)-1, numel(ifreq)-1); % {kats x time x freq} [ch1 ch2] store all indexes of channels in pairs for each enlarged bin
            for k = 1:nkats
                for ifr = 1:numel(ifreq)-1
                    for it = 1:numel(itime)-1
                        bins_n_chanPair{k}(it,ifr) = sum(sum(squeeze(ispc_cats_bins{k}(:, :, it, ifr))~=0 & ~isnan(squeeze(ispc_cats_bins{k}(:, :, it, ifr)))))/2; % as matrix is symmetrical, we shoud divide by 2 
                        [ch1, ch2] = find(triu(squeeze(ispc_cats_bins{k}(:, :, it, ifr))~=0 & ~isnan(squeeze(ispc_cats_bins{k}(:, :, it, ifr))))); % find pairs of chan with non-zero and not NaN values, without repetition (upper triangle of the matrix) 
                        bins_idx_chanPairs{k,it,ifr} = [ch1, ch2];
                    end
                end
            end
            
            % find first N bins with the maximum number of chan pairs with significant values
            max_bins = cell(nkats,1); %{kats}[itime ifreq time freq n_chan_pairs]
            for k = 1:nkats
                bins_n_chanPair_kat = bins_n_chanPair{k};
                [n_chan_pairs, imax] = maxk(bins_n_chanPair_kat(:), n_max_bins); % linear indices
                [itime_max, ifreq_max] = ind2sub(size(bins_n_chanPair_kat), imax); % converting linear indices to row and column indices
                max_bins{k} = table();
                max_bins{k}.itime = itime_max;
                max_bins{k}.ifreq = ifreq_max; % these indexes then help to find exact chan indexes in bins_idx_chanPairs  
                max_bins{k}.time_sec = transpose(round((itime(itime_max)-abs(iepochtime(1)))/obj.E.fs,1)); % real time value in sec relative to stimulus time
                max_bins{k}.freq_hz = transpose(obj.E.Hf(ifreq(ifreq_max))); % real freq in Hz
                max_bins{k}.N_chan_pairs = n_chan_pairs; % number of chan pairs in this bin
            end
            
            obj.plotISPC.PlotISPCbins.ispc_cats_bins = ispc_cats_bins; % save all to the structure
            obj.plotISPC.PlotISPCbins.bins_n_chanPair = bins_n_chanPair;
            obj.plotISPC.PlotISPCbins.bins_idx_chanPairs = bins_idx_chanPairs;
            obj.plotISPC.PlotISPCbins.max_bins = max_bins;           
        end
      
        function [ispc_cats_intime] = ISPCAverageFreq(obj, freq, signChanPair, only_significant_values)  % Sofia from 13.6.2022 
            % averages time-frequency matrices of ISPC across frequencies for each condition - to obtain one curve in time for each pair of chan            
            % freq - frequency range over which to average, e.g. [4 8]
            % signChanPair - if 1, averages only for channel pairs with significant ispc values (each condition separately), if 0 - averages across all between ROI chan pairs
            % only_significant_values = if 1, averages only significant ispc values in obj.plotISPC.ispc_cats_clust_corr; 
            % if 0, averages all ispc values in obj.plotISPC.ispc_cats          
            % output argument ispc_cats_intime - {kats}(chn1, chn2, time points) averaged ispc data across frequencies for each channel pair
            
            if ~isfield(obj.plotISPC,'PlotChPair')
                obj.plotISPC.PlotChPair = struct; % create a separate struct to save all data related to the ISPCPlotChPair plot here
            end
                        
            if ~exist('freq', 'var') || isempty(freq)
                freq = [obj.E.Hf(1) obj.E.Hf(end)]; % default, average over all freq
            end
            obj.plotISPC.PlotChPair.freq = freq; % save to use it in the future plot
            
            if ~exist('signChanPair', 'var') || isempty(signChanPair)
                signChanPair = 1; % default, averages across only channel pairs with significant ispc values
            end
            
            if ~exist('only_significant_values', 'var') || isempty(only_significant_values)
                only_significant_values = 1; % default, averages only significant ispc values                
            end
            if only_significant_values
                ispc_cats = obj.plotISPC.ispc_cats_clust_corr; % {kats}(chn x chn x time x fq) significant ispc values for each condition (non-signif = 0)
            else
                ispc_cats = obj.plotISPC.ispc_cats; % {kats}(chn x chn x time x fq) raw ispc values for each condition
            end
            
            [~, ifreq] = ismember(freq, obj.E.Hf); % indexes of desired freq in the data
            nchannels = size(ispc_cats{1},1); % number of channels
            nkats = size(ispc_cats, 1); % number of conditions
            ispc_cats_intime = cell(nkats,1); %{kats}(chn1, chn2, time points)
            i = 1; % index for each chan pair
            for ichn1 = 1:nchannels
                for ichn2 = 1:nchannels
                    if ismember([ichn1, ichn2], obj.plotISPC.roi_chanpair,'rows') % if chan pair is between ROI
                        for k = 1:nkats
                            if signChanPair && ~(ismember([ichn1, ichn2],obj.plotISPC.chanPair_sign_filtered{k},'rows') || ismember([ichn2, ichn1],obj.plotISPC.chanPair_sign_filtered{k},'rows'))
                                continue % if the current pair of chan is not signif, skip it
                            end
                            temp = mean(squeeze(ispc_cats{k}(ichn1,ichn2,:,ifreq(1):ifreq(2))),2); % mean over selected frequency bins
                            ispc_cats_intime{k}(i,:) = [ichn1, ichn2, temp']; % chan pair: ISPC in time
                        end
                        i = i+1;
                    end
                end
            end
            
            % find rows with all zeros (chan pairs which were skipped)
            for k = 1:nkats
                skipped_chanpairs = all(ispc_cats_intime{k} == 0, 2);               
                ispc_cats_intime{k}(skipped_chanpairs, :) = []; % remove them
            end            
            obj.plotISPC.PlotChPair.ispc_cats_intime = ispc_cats_intime;
        end
        
        function ISPCPlotChPair(obj, chpair) 
            % Sofiia since 14.6.2022
            % plots mean ISPC over frequency bins for each channel pair and condition
            % chpair - index of channel pair which we want to plot
           
            assert(isfield(obj.plotISPC.PlotChPair,'ispc_cats_intime') && ~isempty(obj.plotISPC.PlotChPair.ispc_cats_intime), 'no ispc values were averaged across frequencies, first call PN.ISPCAverageFreq()');

            if ~isfield(obj.plotISPC.PlotChPair,'chpairs')
                chpairs = obj.plotISPC.PlotChPair.ispc_cats_intime{1, 1}(:,[1 2]); % all channel pairs with ISPC averaged
                obj.plotISPC.PlotChPair.chpairs = chpairs;
            end
            
            if ~exist('chpair','var') || isempty(chpair) % selected channel pair
               if isfield(obj.plotISPC.PlotChPair,'chpair') && ~isempty(obj.plotISPC.PlotChPair.chpair) 
                   chpair = obj.plotISPC.PlotChPair.chpair; 
                else                     
                   chpair = 1; % the fisrt pair of channels
                   obj.plotISPC.PlotChPair.chpair = chpair; 
                end
            else
                obj.plotISPC.PlotChPair.chpair = chpair;
            end
            
            if ~isfield(obj.plotISPC.PlotChPair,'kats')
                kats = obj.plotISPC.ispc_cats_names{1,2}; % original numbers of categories
                obj.plotISPC.PlotChPair.kats = kats;
            else
                kats = obj.plotISPC.PlotChPair.kats;
            end
            
            if ~isfield(obj.plotISPC.PlotChPair,'kats_legend')
                kats_legend = obj.plotISPC.ispc_cats_names{1,1}; % names of categories
                obj.plotISPC.PlotChPair.kats_legend = kats_legend;
            else
                kats_legend = obj.plotISPC.PlotChPair.kats_legend;
            end
            
            T = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1)); % time on x-axis            
            
            % ploting
            if isfield(obj.plotISPC.PlotChPair,'h') && ~isempty(obj.plotISPC.PlotChPair.h) && ishandle(obj.plotISPC.PlotChPair.h)
                figure(obj.plotISPC.PlotChPair.h) % use the existing plot
                clf(obj.plotISPC.PlotChPair.h); % clear it
            else
                obj.plotISPC.PlotChPair.h = figure('Name',['ISPCPlotChPair mean over frequency ' num2str(obj.plotISPC.PlotChPair.freq(1)) '-' num2str(obj.plotISPC.PlotChPair.freq(2)) 'Hz']); % or new plot
            end
            
            if numel(kats)>3 % to get colours for all conditions
                colors = distinguishable_colors(numel(kats));
            else
                colors = cell2mat(obj.E.colorskat(1:end)');
            end
            
            h_kat = zeros(numel(kats),1); % handles of plots               
            for k = 1 : numel(kats)
                kk = obj.E.KatIndex(kats,k); % index of condition with correct color code
                ISPCchpair_k = obj.plotISPC.PlotChPair.ispc_cats_intime{k}(chpair,3:end); % plot curve
                h_kat(k,1) = plot(T,ISPCchpair_k,'LineWidth',2,'Color',colors(kk,:));
                hold on                
            end
            
            if ~isfield(obj.plotISPC.PlotChPair,'ylim') || isempty(obj.plotISPC.PlotChPair.ylim)
                yrange = [0 max(obj.plotISPC.PlotChPair.ispc_cats_intime{1}(chpair,3:end))+0.15];
                obj.plotISPC.PlotChPair.ylim = yrange;
            else
                yrange = obj.plotISPC.PlotChPair.ylim; % if we changed the y axis manually on the plot by pressing *
            end            
            set(gca, 'Ylim', yrange);
            ylabel('ISPC value')
            xlabel('time, sec')
            xrange = get(gca,'XLim');           
            plot([0 0], yrange, ':', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.5); % plot stimulus time
            hold on
            legend(h_kat, kats_legend,'Interpreter', 'none','Location','best')
            text(xrange(1)+0.1,yrange(2)-0.025, [obj.E.CH.H.channels(obj.plotISPC.PlotChPair.chpairs(chpair,1)).neurologyLabel '-' obj.E.CH.H.channels(obj.plotISPC.PlotChPair.chpairs(chpair,2)).neurologyLabel]) % neurology labels of channels
            text(xrange(1)+0.1,yrange(2)-0.04, [obj.E.CH.brainlabels(obj.plotISPC.PlotChPair.chpairs(chpair,1)).label '-' obj.E.CH.brainlabels(obj.plotISPC.PlotChPair.chpairs(chpair,2)).label]) % brain labels of channels
            text(xrange(1)+0.1,yrange(2)-0.01,sprintf(' channels %i, %s - %i, %s', ...
                    obj.plotISPC.PlotChPair.chpairs(chpair,1),obj.E.CH.H.channels(obj.plotISPC.PlotChPair.chpairs(chpair,1)).name,...
                    obj.plotISPC.PlotChPair.chpairs(chpair,2),obj.E.CH.H.channels(obj.plotISPC.PlotChPair.chpairs(chpair,2)).name)); % number of channel and its name
            title(['channel pair ' num2str(chpair) '/' num2str(length(obj.plotISPC.PlotChPair.chpairs))], 'Interpreter', 'none'); % number of chan pair
            
            methodhandle = @obj.hybejPlotChPairISPC; % switch chan pairs 
            set(obj.plotISPC.PlotChPair.h,'KeyPressFcn',methodhandle);
            figure(obj.plotISPC.PlotChPair.h); % activates this figure again
        end
        
        function ISPCPlotCatROIMean(obj, ROI1, ROI2) % Sofiia since 13.6.2022
            % plots mean ISPC (which was averaged also across freq) and std err of mean over all pairs of chan between 2 ROIs for each condition
            % ROI1 - the name of ROI1, should correspond to the name in brainlabels, e.g. 'LTC'
            % ROI2 - the name of ROI2, should correspond to the name in brainlabels, e.g. 'IPL'

            assert(exist('ROI1', 'var') && ~isempty(ROI1) && exist('ROI2', 'var') && ~isempty(ROI2), 'specify ROI1 and ROI2')
            [labels_uniq,~] = unique({obj.E.CH.brainlabels(:).label}, 'stable'); % all ROI=brainlabels in the file
            labels_uniq = lower(labels_uniq); % make all names lower
            ROI1 = lower(ROI1);
            ROI2 = lower(ROI2);
            assert(ismember(ROI1, labels_uniq) & ismember(ROI2, labels_uniq), 'ROI should correspond to the names in brainlabels');
            
            if ~isfield(obj.plotISPC.PlotChPair,'chpairs')
                chpairs = obj.plotISPC.PlotChPair.ispc_cats_intime{1, 1}(:,[1 2]); % all channel pairs with ISPC averaged
                obj.plotISPC.PlotChPair.chpairs = chpairs;
            else
                chpairs = obj.plotISPC.PlotChPair.chpairs; % npairs x 2: [ch1 ch2] 
            end
            
            labels = {obj.E.CH.brainlabels(:).label}'; % all brainlabels for all channels                                      
            ROI1_chan_idx = find(contains(lower(labels),ROI1)); % channel indexes for ROI1         
            ROI2_chan_idx = find(contains(lower(labels),ROI2)); % channel indexes for ROI2
            
            % all possible pairs between 2 ROIs
            ROI_pairs = zeros(length(ROI1_chan_idx) * length(ROI2_chan_idx), 2);
            count = 1;            
            for i = 1:length(ROI1_chan_idx)
                for j = 1:length(ROI2_chan_idx)
                    ROI_pairs(count, :) = [ROI2_chan_idx(j) ROI1_chan_idx(i)];
                    count = count + 1;
                end
            end
            [~, ipairs2average] = ismember(ROI_pairs, chpairs, 'rows');
            ipairs2average = sort(ipairs2average(ipairs2average ~= 0)); % filter out all zeros and sort
            
            nChanPairs = numel(ipairs2average); % number of pairs of channels
            nkats = size(obj.plotISPC.PlotChPair.ispc_cats_intime, 1);  % number of conditions         
            ispc_cats_roi_mean = zeros(nkats,size(obj.plotISPC.PlotChPair.ispc_cats_intime{1},2)-2); % mean, kats x time (chan pairs will be averaged)
            ispc_cats_roiE = zeros(nkats,size(obj.plotISPC.PlotChPair.ispc_cats_intime{1},2)-2); % std err of mean, kats x time
            
            for k = 1:nkats
                isps_time_k = obj.plotISPC.PlotChPair.ispc_cats_intime{k}(ipairs2average, 3:end); % only time points without number of chan
                ispc_cats_roi_mean(k,:) = mean(isps_time_k, 1); % average across selected pairs of chan, kats x time points
                ispc_cats_roiE(k,:) = std(isps_time_k,[],1)/sqrt(nChanPairs); % std err of mean
            end
            
            % ploting
            T = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1)); % time on x-axis
            figure('Name','ISPCPlotCatMeanROI')
            kats = obj.plotISPC.ispc_cats_names{1,2}; % original numbers of categories
            hue = 0.8;
            colorsErrorBars = cellfun(@(a) min(a+hue, 1), obj.E.colorskat, 'UniformOutput', false);
            
            ploth = zeros(1,nkats); % handles of plots
            for k=1:nkats % for each category
                kk = obj.E.KatIndex(kats,k);
                colorkatk = [obj.E.colorskat{kk} ; colorsErrorBars{kk}]; % two colours, for the main curve and error bar
                ciplot(ispc_cats_roi_mean(k,:)+ispc_cats_roiE(k,:), ispc_cats_roi_mean(k,:)-ispc_cats_roiE(k,:), T, colorkatk(2,:)); % plot std err of mean
                hold on
                ploth(k) = plot(T,ispc_cats_roi_mean(k,:),'LineWidth',1,'Color',colorkatk(1,:)); % plot mean
                hold on
            end           
            set(gca, 'Ylim', [min(min(ispc_cats_roi_mean))-0.01 max(max(ispc_cats_roi_mean))+0.01]);       
            ylabel('ISPC value')
            xlabel('time, sec')
            yrange = get(gca, 'Ylim');
            plot([0 0], yrange, ':', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.5); % plot stimulus time
            legend(ploth,obj.plotISPC.ispc_cats_names{1,1},'Interpreter', 'none','Location','best');           
            title(['ISPC mean over ' num2str(nChanPairs) ' pairs of channels in ' upper(ROI1) ' - ' upper(ROI2) ' in frequencies: ' num2str(obj.plotISPC.PlotChPair.freq(1)) '-' num2str(obj.plotISPC.PlotChPair.freq(2)) ' Hz'])
        end
        
        function ISPCPlot(obj, cat, chnpair, TimeHf, updateS, ISPCPlotSelection)
            % plot all ispc values  for a specific pair of channels from computed values for one category
            % cat - one category, which we want to plot (one number)
            % chnpair - pair of the channels (or rather their xy coordinates in the plot) to be higlited on the left top plot of all channels, and shown an all other figure
            % TimeHf - point in the time-frequency space highlighted on the left top plot - 
           
            assert(isfield(obj.plotISPC,'ispc_cats')&& ~isempty(obj.plotISPC.ispc_cats), 'no ispc values are computed');
            
            % define one category = condition, which we want to plot
            if ~exist('cat', 'var') || isempty(cat)
                if isfield(obj.plotISPC,'ispc_cats_names') && ~isempty(obj.plotISPC.ispc_cats_names)
                    cat = obj.plotISPC.ispc_cats_names{1, 2}(1); % take the first category of all existing
                else
                    cat = obj.E.Wp(obj.E.WpActive).kats(1); % take the first category of all existing
                end
            end 
            obj.plotISPC.current_cat = cat; % save category to data to use later when updating the figure
            ik = find(cat==obj.plotISPC.ispc_cats_names{1, 2}); % find index of category in the data           
                    
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
                
%             if isnan(max(max(squeeze(obj.plotISPC.ispc_cats{ik}(chn1,chn2,:,:))))) %if the order of the channels is reversed - the ispc values are stored only once
%                 chn1 = obj.plotISPC.sortorder(chnpair(2));
%                 chn2 = obj.plotISPC.sortorder(chnpair(1));
%             end
            if ~exist('TimeHf','var') || isempty(TimeHf)
                if obj.plotISPC.ispc_sign_chanPair{ik}(chn1,chn2)==0 || isnan(obj.plotISPC.ispc_sign_chanPair{ik}(chn1,chn2))
                    TimeHf = [round(size(obj.plotISPC.ispc_cats{ik},3)/2) round(size(obj.plotISPC.ispc_cats{ik},4)/2)]; % point in the middle
                else
                    [iTime,iHf] = find(squeeze(obj.plotISPC.ispc_cats_clust_corr{ik}(chn1,chn2,:,:)) == max(max(squeeze(obj.plotISPC.ispc_cats_clust_corr{ik}(chn1,chn2,:,:)))),1);
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
%                 obj.plotISPC.plotsig = [0.95 1]; % limits of significance plots
                obj.plotISPC.onlysub2 = 0; % plot only second chart with an overview of all channels
                obj.plotISPC.clf = 0; % clear the whole plot to plot all subplots after only second chart was plotted                
                obj.plotISPC.sortorder = 1:size(obj.plotISPC.ispc_cats{ik},1); %originally channels sorted by default
                obj.plotISPC.sortby = 0; %sorted originally, 1=sorted by brainlabel
                obj.plotISPC.els = obj.E.CH.els; %borders of the electrodes/pacients
%                 obj.plotISPC.aplha = 0.4; %transparency for non signif IPSC values 1=non transparent, 0= completely transparent
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
 
            %subplot with ISPC values from one pair of channels (by countour, significant uncorrected ispc values are highlighed)
            if updateS(1)
                if ~obj.plotISPC.onlysub2
                    obj.plotISPC.subph{1} = subplot(2,2,1);
                end
                ispc=squeeze(obj.plotISPC.ispc_cats{ik}(chn1,chn2,:,:))'; % all ispc values
                ispc_signif = squeeze(obj.plotISPC.ispc_cats_signif{ik}(chn1,chn2,:,:))'; % only significant uncorrected ispc values
                % im1=imagesc(T,F,ispc); %has to transpose, imagesc probably plots first dimension in rows
                contourf(T,F,ispc,40,'linecolor','none'); % plot all ispc values 
                hold on
                contour(T,F,ispc_signif,1,'linecolor','k', 'LineWidth', 0.01) % on the same chart, highlight significant uncorrected ispc values
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
                obj.plotISPC.annoth{1} = annotation(obj.plotISPC.fh_ispc,'textbox',[0 .9 .4 .1],'String',num2str(obj.plotISPC.ispc_cats{ik}(chn1,chn2, TimeHf(1),TimeHf(2))), 'EdgeColor', 'none');
            end
            
            %subplot with ISPC values from one frequency of one pair of channels
            if updateS(3) 
                if ~obj.plotISPC.onlysub2
                    obj.plotISPC.subph{3} = subplot(2,2,3); 
                end
                plot(T,squeeze(obj.plotISPC.ispc_cats{ik}(chn1,chn2, :,TimeHf(2))));
                hold on;
                plot(T(TimeHf(1)),obj.plotISPC.ispc_cats{ik}(chn1,chn2, TimeHf(1),TimeHf(2)),'o');
                title([num2str(F(TimeHf(2))) ' Hz']);
                if isnan(obj.plotISPC.ispc_cats{ik}(chn1,chn2,:,:));
                    maxispc = 0;
                else
                    maxispc = max(max(squeeze(obj.plotISPC.ispc_cats{ik}(chn1,chn2,:,:))));
                end
                ylim([-0.2 iff(maxispc==0,1,maxispc)]);            
            end
            
            %subplot with minimum ISPC significance across all channels
            if updateS(2)                
                if ~obj.plotISPC.onlysub2
                    obj.plotISPC.subph{2} = subplot(2,2,2); 
                end
                ispc_sign_chanPair = obj.plotISPC.ispc_sign_chanPair{ik}(obj.plotISPC.sortorder,obj.plotISPC.sortorder); %values sorted according to sortorder                
%                 [x,y]=find(isnan(obj.plotISPC.ispc_p_min));
%                 ispc_p_min(sub2ind(size(ispc_p_min),x,y)) = ispc_p_min(sub2ind(size(ispc_p_min),y,x)); %make the matrix symmetrical - each value will be there twice 
%                 im2=imagesc(obj.plotISPC.channels,obj.plotISPC.channels, ispc_p_min',obj.plotISPC.plotsig); 
                im2=imagesc(obj.plotISPC.channels,obj.plotISPC.channels, ispc_sign_chanPair');
                % adding a color for NaN values
                imAlpha=ones(size(ispc_sign_chanPair'));  %https://www.mathworks.com/matlabcentral/answers/81938-set-nan-as-another-color-than-default-using-imagesc
                imAlpha(isnan(ispc_sign_chanPair'))=0;
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
                t3 = obj.plotISPC.ispc_cats_names{1, 1}{ik}; % name of category
                title( iff(obj.plotISPC.onlysub2,{t1,t2, t3},t1)); %two line title for this large plot
            end
            
            %subplot with significance of ISPC of this pair of channels after cluster correction
            if updateS(4)
                if ~obj.plotISPC.onlysub2
                    obj.plotISPC.subph{4} = subplot(2,2,4);
                end
%                 T4 = linspace(0, obj.E.epochtime(2), size(obj.plotISPC.ispc_cats{ik},3));
                %                 ispc_p = 1-squeeze(obj.plotISPC.ispc_p(chn1,chn2,:,:))';
                ispc_clust_corr = squeeze(obj.plotISPC.ispc_cats_clust_corr{ik}(chn1,chn2,:,:)); % plot significant ispc values after cluster correction
%                 im4=imagesc(T4,F,ispc_p,obj.plotISPC.plotsig);  %0.5
%                 imAlpha=ones(size(ispc_p));  %https://www.mathworks.com/matlabcentral/answers/81938-set-nan-as-another-color-than-default-using-imagesc
%                 imAlpha(ispc_p<obj.plotISPC.plotsig(1))=obj.plotISPC.aplha;
%                 im4.AlphaData = imAlpha; %set transparency for nan values
%                 imAplhaISPC = cat(2,repmat(obj.plotISPC.aplha,size(imAlpha,1),size(ispc,2)-size(ispc_p,2)),imAlpha);
%                 im1.AlphaData = imAplhaISPC;
                contourf(T,F,ispc_clust_corr',40,'linecolor','none')
                xlabel('Time (s)'), ylabel('Frequency (Hz)')
                
                axis xy;
                colorbar;
                hold on;
                plot(T(TimeHf(1)),F(TimeHf(2)),'o','MarkerFaceColor','red');
                title(sprintf('channel %i (%s) x %i (%s)',obj.plotISPC.channels(chn1),obj.E.CH.H.channels(obj.plotISPC.channels(chn1)).neurologyLabel,...
                    obj.plotISPC.channels(chn2), obj.E.CH.H.channels(obj.plotISPC.channels(chn2)).neurologyLabel));
            end
            
            % Adding a title above all subplots, showing the condition 
            axes( 'Position', [0, 0.95, 1, 0.05] ) ;
            set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
            text( 0.5, 0, obj.plotISPC.ispc_cats_names{1, 1}{ik}, 'Interpreter', 'none', 'FontSize', 14', 'FontWeight', 'Bold', ...
            'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;          
            
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
        function ISPCPlotEnlargedBins(obj, chpair, ikat, imax_bin) % Sofiia since 3.11.2023
            % plots new enlarged time-frequency bins of significant ispc values for each channel pair (only one condition)
            % chpair - index of channel pair which we want to plot
            % ikat - index of condition
            % imax_bin - index in obj.plotISPC.PlotISPCbins.max_bins{ikat} to choose which chan pairs to plot
            assert(isfield(obj.plotISPC,'PlotISPCbins')&& ~isempty(obj.plotISPC.PlotISPCbins), 'first call ISPCFindSignBins to create enlarged bins');
            
            if ~exist('chpair','var') || isempty(chpair) % selected channel pair
                if isfield(obj.plotISPC.PlotISPCbins,'chpair') && ~isempty(obj.plotISPC.PlotISPCbins.chpair)
                    chpair = obj.plotISPC.PlotISPCbins.chpair;
                else
                    chpair = 1; % the fisrt pair of channels
                    obj.plotISPC.PlotISPCbins.chpair = chpair;
                end
            else
                obj.plotISPC.PlotISPCbins.chpair = chpair;
            end
            
            if ~exist('ikat','var') || isempty(ikat) % selected condition - its index               
                if isfield(obj.plotISPC.PlotISPCbins,'ikat')
                    ikat = obj.plotISPC.PlotISPCbins.ikat;
                else
                    ikat = 1; % first condition
                    obj.plotISPC.PlotISPCbins.ikat = ikat;
                end
            else
                obj.plotISPC.PlotISPCbins.ikat = ikat;
            end
            
            if ~exist('imax_bin','var') || isempty(imax_bin) % index of window with max number of chan pairs               
                if isfield(obj.plotISPC.PlotISPCbins,'imax_bin')
                    imax_bin = obj.plotISPC.PlotISPCbins.imax_bin;
                else
                    imax_bin = 1; % the first max bin
                    obj.plotISPC.PlotISPCbins.imax_bin = imax_bin;
                end
            else
                obj.plotISPC.PlotISPCbins.imax_bin = imax_bin;
            end
            
            % take all chan pairs in this bin
            itime = obj.plotISPC.PlotISPCbins.max_bins{ikat}.itime(imax_bin);
            ifreq = obj.plotISPC.PlotISPCbins.max_bins{ikat}.ifreq(imax_bin);
            chpairs = obj.plotISPC.PlotISPCbins.bins_idx_chanPairs{ikat,itime,ifreq}; % {kats x time x freq} [ch1 ch2]; all channel pairs with ISPC averaged for this bin
            obj.plotISPC.PlotISPCbins.chpairs = chpairs;
                        
            % ploting
            freq_bin = obj.plotISPC.PlotISPCbins.freq(2) - obj.plotISPC.PlotISPCbins.freq(1); % size of freq-time bin
            time_bin = obj.plotISPC.PlotISPCbins.time(2) - obj.plotISPC.PlotISPCbins.time(1);
            if isfield(obj.plotISPC.PlotISPCbins,'h') && ~isempty(obj.plotISPC.PlotISPCbins.h) && ishandle(obj.plotISPC.PlotISPCbins.h)
                figure(obj.plotISPC.PlotISPCbins.h) % use the existing plot
                clf(obj.plotISPC.PlotISPCbins.h); % clear it
            else
                obj.plotISPC.PlotISPCbins.h = figure('Name',['ISPCPlotEnlargedBins with bin size:  ' num2str(freq_bin) ' Hz and ' num2str(time_bin) ' sec']); % or new plot
            end            
            
            T = obj.plotISPC.PlotISPCbins.time;
            F = obj.plotISPC.PlotISPCbins.freq(1:end-1);
            chn1 = obj.plotISPC.PlotISPCbins.chpairs(chpair,1);
            chn2 = obj.plotISPC.PlotISPCbins.chpairs(chpair,2);
            ispc_bins = squeeze(obj.plotISPC.PlotISPCbins.ispc_cats_bins{ikat}(chn1,chn2,:,:)); % binned data for this chan pair
            imagesc(T, F, ispc_bins');               
            axis xy;
            set(gca, 'XTick',T-time_bin/2, 'XTickLabel',T) % shift ticks to the left
            set(gca, 'YTick',F-freq_bin/2, 'YTickLabel',F)
            xlabel('Time (s)'), ylabel('Frequency (Hz)')  
            colorbar;  
            
            yrange = get(gca, 'Ylim');
            xrange = get(gca,'XLim');
            
            text(xrange(1)+0.1,yrange(2)-0.2, sprintf('channel %i %s x %i %s',obj.plotISPC.channels(chn1),obj.E.CH.H.channels(obj.plotISPC.channels(chn1)).neurologyLabel,...
                obj.plotISPC.channels(chn2), obj.E.CH.H.channels(obj.plotISPC.channels(chn2)).neurologyLabel), 'Color' ,'r'); % neurology label
            text(xrange(1)+0.1,yrange(2)-1, [obj.E.CH.brainlabels(obj.plotISPC.PlotISPCbins.chpairs(chpair,1)).label '-' obj.E.CH.brainlabels(obj.plotISPC.PlotISPCbins.chpairs(chpair,2)).label], 'Color' ,'r') % brain labels of channels
            text(xrange(1)+0.1,yrange(2)-0.6,sprintf(' %s x %s', ...
                    obj.E.CH.H.channels(obj.plotISPC.PlotISPCbins.chpairs(chpair,1)).name, obj.E.CH.H.channels(obj.plotISPC.PlotISPCbins.chpairs(chpair,2)).name), 'Color' ,'r'); %  channels name   
            title(['channel pair ' num2str(chpair) '/' num2str(length(obj.plotISPC.PlotISPCbins.chpairs)) ', condition: ' obj.plotISPC.ispc_cats_names{1, 1}{ikat}], 'Interpreter', 'none'); % number of chan pair and the name of condition
           
            methodhandle = @obj.hybejISPCPlotEnlargedBins; % switch chan pairs 
            set(obj.plotISPC.PlotISPCbins.h,'KeyPressFcn',methodhandle);
            figure(obj.plotISPC.PlotISPCbins.h); % activates this figure again             
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
        
        function obj = hybejPlotChPairISPC(obj,~,eventDat) %Sofiia 14.6.2022
            % reacts to events in figure ISPCPlotChPair
            switch eventDat.Key
                case {'rightarrow'} % next channel pair by ->
                    chpair = min([obj.plotISPC.PlotChPair.chpair + 1 , length(obj.plotISPC.PlotChPair.chpairs)]);
                    obj.ISPCPlotChPair(chpair);
                case {'leftarrow'} % previous channel pair by <-
                    chpair = max([obj.plotISPC.PlotChPair.chpair - 1 , 1]);
                    obj.ISPCPlotChPair(chpair);
                case 'pagedown' % move by 10 ch pairs forward
                    chpair = min([obj.plotISPC.PlotChPair.chpair + 10 , length(obj.plotISPC.PlotChPair.chpairs)]);
                    obj.ISPCPlotChPair(chpair);
                case 'pageup' % move by 10 ch pairs back
                    chpair = max([obj.plotISPC.PlotChPair.chpair - 10 , 1]);
                    obj.ISPCPlotChPair(chpair);
               case 'home' % choose the first chan pair
                    obj.ISPCPlotChPair(1);                    
                case 'end' % choose the last chan pair
                    obj.ISPCPlotChPair(length(obj.plotISPC.PlotChPair.chpairs)); 
                case 'o' % specify the chan pair by input number
                     answ = inputdlg('Enter channel pair number:','Go to channel pair', 1,{num2str(obj.plotISPC.PlotChPair.chpair)});
                     if numel(answ)>0
                         chpair = str2double(answ{1}); 
                         if ~isempty(chpair), obj.ISPCPlotChPair(chpair); end
                     end
                case {'multiply','8'} % asterisk on the numerical keyboard, or asterisk above 8
                    % window to choose min and max value of y axis
                    answ = inputdlg('Enter ymax and min:','Yaxis limits', [1 50],{num2str(obj.plotISPC.PlotChPair.ylim)});
                    if numel(answ)>0  % answer is cell 1x1 - if it was cancelled cell 0x0
                        if isempty(answ{1}) || any(answ{1}=='*') % if it is empty or contains *, recompute ylim again
                            obj.plotISPC.PlotChPair.ylim = [];
                        else 
                            data = str2num(answ{:});  % otherwise change ylim according to input numbers
                            if numel(data)>= 2 && data(1)< data(2)
                                obj.plotISPC.PlotChPair.ylim = [data(1) data(2)];
                            end
                        end
                    end
                    obj.ISPCPlotChPair();
            end
        end
        function obj = hybejISPCPlotEnlargedBins(obj,~,eventDat) % Sofiia since 3.11.2023
            % reacts to events in figure ISPCPlotEnlargedBins
            switch eventDat.Key
                case {'rightarrow'} % next channel pair by ->
                    chpair = min([obj.plotISPC.PlotISPCbins.chpair + 1 , length(obj.plotISPC.PlotISPCbins.chpairs)]);
                    obj.ISPCPlotEnlargedBins(chpair);
                case {'leftarrow'} % previous channel pair by <-
                    chpair = max([obj.plotISPC.PlotISPCbins.chpair - 1 , 1]);
                    obj.ISPCPlotEnlargedBins(chpair);
                case 'pagedown' % move by 10 ch pairs forward
                    chpair = min([obj.plotISPC.PlotISPCbins.chpair + 10 , length(obj.plotISPC.PlotISPCbins.chpairs)]);
                    obj.ISPCPlotEnlargedBins(chpair);
                case 'pageup' % move by 10 ch pairs back
                    chpair = max([obj.plotISPC.PlotISPCbins.chpair - 10 , 1]);
                    obj.ISPCPlotEnlargedBins(chpair);
               case 'home' % choose the first chan pair
                    obj.ISPCPlotEnlargedBins(1);                    
                case 'end' % choose the last chan pair
                    obj.ISPCPlotEnlargedBins(length(obj.plotISPC.PlotISPCbins.chpairs)); 
                case 'o' % specify the chan pair by input number
                     answ = inputdlg('Enter channel pair number:','Go to channel pair', 1,{num2str(obj.plotISPC.PlotISPCbins.chpair)});
                     if numel(answ)>0
                         chpair = str2double(answ{1}); 
                         if ~isempty(chpair), obj.ISPCPlotEnlargedBins(chpair); end
                     end
            end            
        end
        
        function obj = hybejPlotISPC(obj,~,eventDat)  
            %reacts to events in figure ISPCPlot
            if ~isempty(eventDat.Modifier) && strcmp(eventDat.Modifier{:},'shift')
                shift = 10;
            else
                shift = 1;
            end
            ik = find(obj.plotISPC.current_cat==obj.plotISPC.ispc_cats_names{1, 2}); % index of current category
            switch eventDat.Key
                case {'numpad6','d'} %next time sample
                    if obj.plotISPC.TimeHf(1) <= size(obj.plotISPC.ispc_cats{ik},3)-shift %compare to num of time samples
                        TimeHf = [(obj.plotISPC.TimeHf(1) + shift), obj.plotISPC.TimeHf(2)];
                        obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, TimeHf);
                    end
                case {'numpad4','a'} %previous time sample   
                    if obj.plotISPC.TimeHf(1) > shift
                        TimeHf = [(obj.plotISPC.TimeHf(1) - shift), obj.plotISPC.TimeHf(2)];
                        obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, TimeHf)
                    end         
                case {'numpad8','w'} %next frequency
                    if obj.plotISPC.TimeHf(2) <= size(obj.plotISPC.ispc_cats{ik},4) - shift   %compare to num of freq 
                        TimeHf = [obj.plotISPC.TimeHf(1) , (obj.plotISPC.TimeHf(2)+shift)];
                        obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, TimeHf)
                    end
                case {'numpad2','s'} %previous frequency 
                    if obj.plotISPC.TimeHf(2) > shift
                        TimeHf = [obj.plotISPC.TimeHf(1) , (obj.plotISPC.TimeHf(2)-shift)];
                        obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, TimeHf)
                    end                
                case {'home'}                     
                    obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, []);
                case {'return'}
                    if obj.plotISPC.onlysub2 == 1
                       obj.plotISPC.onlysub2 = 0;
                       obj.plotISPC.clf = 1;
                    end
                    obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, [], [1 1 1 1]);
                case {'rightarrow'} %move the channel pair to right
                    if obj.plotISPC.chnpair(1) <= numel(obj.plotISPC.channels) -shift
                        chnpair = [obj.plotISPC.chnpair(1)+shift obj.plotISPC.chnpair(2)];
                        obj.ISPCPlot(obj.plotISPC.current_cat, chnpair, obj.plotISPC.TimeHf, [0 1 0 0]);
                    end
                case {'leftarrow'} % move the channel pair to left
                    if obj.plotISPC.chnpair(1) > shift
                        chnpair = [obj.plotISPC.chnpair(1)-shift obj.plotISPC.chnpair(2)];
                        obj.ISPCPlot(obj.plotISPC.current_cat, chnpair, obj.plotISPC.TimeHf, [0 1 0 0]);
                    end
                case {'uparrow'} % zvyseni cisla aktivni statistiky
                    if obj.plotISPC.chnpair(2) <= numel(obj.plotISPC.channels) -shift
                        chnpair = [obj.plotISPC.chnpair(1) obj.plotISPC.chnpair(2)+shift];
                        obj.ISPCPlot(obj.plotISPC.current_cat, chnpair, obj.plotISPC.TimeHf, [0 1 0 0]);
                    end
                case {'downarrow'}
                    if obj.plotISPC.chnpair(2) > shift
                        chnpair = [obj.plotISPC.chnpair(1) obj.plotISPC.chnpair(2)-shift];
                        obj.ISPCPlot(obj.plotISPC.current_cat, chnpair, obj.plotISPC.TimeHf, [0 1 0 0]);
                    end
                case 'q'
                    if obj.plotISPC.plotsig(1) == 0.5
                        obj.plotISPC.plotsig = [0.95 1];
                    elseif obj.plotISPC.plotsig(1) == 0.95 
                        obj.plotISPC.plotsig = [0.99 1];
                    else
                        obj.plotISPC.plotsig = [0.5 0.95];
                    end
                    obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [0 1 0 1]);
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
                    obj.ISPCPlot(obj.plotISPC.current_cat,obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [0 1 0 0]); 
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
                        obj.plotISPC.sortorder = 1:size(obj.plotISPC.ispc_cats{ik},1);
                        obj.plotISPC.sortby = 0;
                    end
                    obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [0 1 0 0]); %update the all channel plot
                case 'subtract'
                    if obj.plotISPC.aplha > 0
                        obj.plotISPC.aplha = obj.plotISPC.aplha - 0.1;
                        obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [1 0 0 1]);
                    end
                case 'add'
                    if obj.plotISPC.aplha < 1
                        obj.plotISPC.aplha = obj.plotISPC.aplha + 0.1;
                        obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, obj.plotISPC.TimeHf, [1 0 0 1]);
                    end
            end
        end
        function hybejPlotISPCClick(obj,~,~)     
            mousept = get(gca,'currentPoint');
            ik = find(obj.plotISPC.current_cat==obj.plotISPC.ispc_cats_names{1, 2}); % index of current category
            if mousept(1,1)<1
                %probably click in the time freq chart
                
                T = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1));
                F =  obj.E.Hfmean; 
                x = find(T>=mousept(1,1),1); 
                y = find(F>=mousept(1,2),1); %coordinated of the click
                if x <= size(obj.plotISPC.ispc_cats{ik},3) && x>=0 && y <= size(obj.plotISPC.ispc_cats{ik},4) && y >= 0%compare to num of time samples                        
                        obj.ISPCPlot(obj.plotISPC.current_cat, obj.plotISPC.chnpair, [x y],[1 0 1 1]);
                end                
            else %click in the channel plot
                x = round(mousept(1,1)); y = round(mousept(1,2)); %coordinated of the click
                if x <= numel(obj.plotISPC.channels) && x>=1 && y<=numel(obj.plotISPC.channels) && y >=1                        
                            updateS = iff(obj.plotISPC.onlysub2 == 0,[1 1 1 1],[0 1 0 0]);
                            obj.ISPCPlot(obj.plotISPC.current_cat, [x y], obj.plotISPC.TimeHf,updateS);
                end
            end
            
        end
         
    end
    
end

