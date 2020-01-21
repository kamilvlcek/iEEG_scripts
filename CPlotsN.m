classdef CPlotsN < handle
    %CPlotsN pokrocile funkce na zobrazeni frekvencnich dat, pouziti napr. PN = CPlotsN(E)
    %   Nada Bednarova 2018/04
    
    properties (Access = public)
        E;                      % CiEEGData object
        HOrigData;              % original untransformed EEG data
        plotData = struct;      % contains figure and plot info
    end
    
    methods (Access = public)
        
        function obj = CPlotsN(E, HOrigData)
            obj.E = E;
            if ~exist('HOrigData', 'var') || isempty(HOrigData)
                obj.HOrigData = downsample(E.d,E.decimatefactor); %do funkce plotUnepochedData
            else
                obj.HOrigData = HOrigData;
            end
        end
        
        function PlotUnepoched(obj, type, psy, channels, frequencies, limits, timeDelay, startTime, plotGroup)
            % Plots powers and phases of transformed eeg signal for :
            %   * all selected frequencies for one channel
            %   * or one frequency for all channels if called from PlotElectrodeGroup
            % INPUT:
            % type          - 0 = power, 1 = original filtered data
            % psy           - psy structure, aka aedist for stimuli and response time
            % channels      - if called directly, only one channel
            % frequencies   - frequencies to plot (not indexes)
            % limits        - y axis limits, default [0,1000]
            % timeDelay     - how many seconds display on the screen - default 10
            % startTime     - second to start from
            % plotGroup     - only if called from PlotElectrodeGroup, default false
            
            obj.plotData.type = type;
            
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
            
            if ~exist('plotGroup','var'); obj.plotData.plotGroup = false; else obj.plotData.plotGroup = plotGroup; end
     
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
                  title(sprintf('%d - epoch %d (%d) ', obj.E.PsyData.P.data(epochs(i),obj.E.PsyData.P.sloupce.opakovani), correct, epochs(i)),'FontSize',8) ;
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
            assert(~isempty(obj.E.HFreqEpochs),'soubor s frekvencnimi daty pro epochy neexistuje');
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
                
        % modro zlte plots synchronizacie fazy image
        function fig = PlotITPCallEpochs(obj, channel)
            fig = figure;
            time = linspace(obj.E.epochtime(1), obj.E.epochtime(2), size(obj.E.fphaseEpochs,1));
            for iFq = 1:length(obj.E.Hf)
                subplot(2,ceil(length(obj.E.Hf)/2),iFq);
                % sort epochs by condition
                epochsByCond = sortrows([[obj.E.epochData{:,2}]' (1:length(obj.E.epochData))'],1);
                phases = squeeze(obj.E.frealEpochs(:,channel,iFq,:))';
                image(phases(epochsByCond(:,2),:), 'XData', time);
                hold on;
                set(gca,'YDir','normal')
                xlabel('Time (s)'); ylabel('Epoch');
                title([num2str(obj.E.Hf(iFq)), ' Hz']);
            end
            mtit(sprintf('Channel %d \n', channel));
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
                    if obj.plotData.plotGroup; 
                        title(strcat(num2str(obj.plotData.iChannels(i)), '. channel (', t, ') ',names{obj.plotData.iChannels(i)},'-',labels{obj.plotData.iChannels(i)}));
                    else
                        title(strcat(num2str(obj.E.Hf(obj.plotData.iFreqs(i))), ' Hz (', t, ')'));
                    end

                    %%%%%%%%%%%%%%%%%%% PHASE %%%%%%%%%%%%%%%%%%% 
                    subplot(numSubplot,2,i*2)
                    % phase values
                    if obj.plotData.plotGroup; 
                        y = squeeze(obj.E.fphase(time,obj.plotData.iChannels(i),obj.plotData.iFreqs)); % channel phase values
                        display('correct, more')
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
         
    end
    
end

