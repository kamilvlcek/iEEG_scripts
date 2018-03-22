classdef CPlotsN < handle
    %CPlotsN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        E;                      % CiEEGData object
        plotFreqs = struct;     % contains figure and plot info
    end
    
    methods (Access = public)
        function obj = CPlotsN(E)
           obj.E = E;
        end
        
        function PlotFrequencies(obj, psy, channels, frequencies, limits, timeDelay, plotGroup)
            % plots powers and phases of transformed eeg signal for :
            %   * all selected frequencies for one channel
            %   * or one frequency for all channels
            assert(isa(psy,'struct'),'Prvy parameter musi byt struktura s datami z psychopy');
            obj.E.PsyData = CPsyData(psy); 
            
            obj.plotFreqs.iChannels = channels;
            
            if ~exist('frequencies','var') || isempty(frequencies); obj.plotFreqs.iFreqs = 1:numel(obj.E.Hf); % ak sme nevybrali frekvencie, zobrazia sa vsetky
                else [~,obj.plotFreqs.iFreqs] = intersect(obj.E.Hf,frequencies); end % indexy vybranych frekvencii
                
            if ~exist('limits','var') || isempty(limits); limits = [0,1000]; end % set y axis limits
            obj.plotFreqs.ylimits = limits;
            
            if ~exist('timeDelay','var'); timeDelay = 10; end % seconds visible on the screen
            obj.plotFreqs.timeDelay = timeDelay*obj.E.fs;
            
            obj.plotFreqs.f = figure('Name','All Frequencies','Position', [20, 100, 1000, 600]); % init the figure
            set(obj.plotFreqs.f, 'KeyPressFcn', @obj.MovePlotFreqs); % set key press function
            
            if ~exist('plotGroup','var'); obj.plotFreqs.plotGroup = false; else obj.plotFreqs.plotGroup = plotGroup; end
            
            obj.plotFreqs.iTime = 0; % initiate time index
     
            if ~isfield(obj.plotFreqs,'plotGroup'); obj.plotFreqs.plotGroup = false; end 
            obj.plotFreqs.stimuli = (obj.E.PsyData.P.data(:,obj.E.PsyData.P.sloupce.ts_podnet) - obj.E.tabs(1))*24*3600;
            obj.plotFreqs.responses = (obj.E.PsyData.P.data(:,obj.E.PsyData.P.sloupce.ts_odpoved) - obj.E.tabs(1))*24*3600;
            obj.plotFreqs.colors = {'black','green','red','blue'};
            
            obj.plotFreqData();
        end
                
        function PlotElectrodeGroup(obj, iGroup, psy, frequency, limits, timeDelay)
            contacts = obj.E.CH.chgroups{iGroup};
            obj.PlotFrequencies(psy, contacts, frequency, limits, timeDelay, true);
        end
        
    end
    
    
    methods (Access = private)
        
        function obj = MovePlotFreqs(obj,~,eventDat)
            %zpracovava stlaceni klavesy pro graf PlotFrequencies
            switch eventDat.Key
                case 'rightarrow' % +1 time window
                    obj.plotFreqs.iTime = obj.plotFreqs.iTime + 1; %min([obj.plotEpochs.iTime + 1, size(obj.HFreqEpochs,4)]);
                case 'leftarrow'  % -1 time window
                    obj.plotFreqs.iTime = max([obj.plotFreqs.iTime - 1, 0]);%max([obj.plotEpochs.iEpoch - 1, 1]);
                otherwise  
                   display(['key pressed: ' eventDat.Key]); %vypise stlacenou klavesu
            end
            obj.plotFreqData();
        end
       
        function plotFreqData(obj)
            % plot all unepoched transformed data
            % one subplot = one frequency
            % called from PlotFrequencies()
            
            time = (obj.plotFreqs.iTime*obj.E.fs+1):(obj.plotFreqs.iTime*obj.E.fs+obj.plotFreqs.timeDelay); % moving time window
            x = time./obj.E.fs; % current x axis in seconds
            % determine number of subplots based on selected channels or frequencies
            if obj.plotFreqs.plotGroup; numSubplot = length(obj.plotFreqs.iChannels); % plot electrode channels 
            else numSubplot = length(obj.plotFreqs.iFreqs)+1; end                     % plot frequencies + original eeg data

            for i = 1:numSubplot 
                %%%%%%%%%%%%%%%%%%% POWER %%%%%%%%%%%%%%%%%%% 
                
                if ~obj.plotFreqs.plotGroup && i > numSubplot-1 % if plotting frequencies, display also original signal on the bottom
                    subplot(numSubplot,2,[i*2-1,i*2])
                    y = obj.E.HOrigData(time,obj.plotFreqs.iChannels); % original eeg values
                    plot(x, y, 'Color', 'black'); hold on;
                    xlim([x(1) x(end)]); % set x axis limits
                    ylim([-80,80]); % set y axis limits
                else
                    subplot(numSubplot,2,i*2-1)
                    if obj.plotFreqs.plotGroup; % data to be subplotted (given channel or frequency)
                        y = squeeze(obj.E.HFreq(time,obj.plotFreqs.iChannels(i),obj.plotFreqs.iFreqs)); % channel power values
                        names = {obj.E.CH.H.channels.ass_brainAtlas}; %cast mozgu
                        labels = {obj.E.CH.H.channels.neurologyLabel}; %neurology label
                    else
                        y = squeeze(obj.E.HFreq(time,obj.plotFreqs.iChannels,obj.plotFreqs.iFreqs(i))); % frequency power values
                    end
                    neg = y<0; % negative frequency power values
                    plot(x(~neg),y(~neg),'b.',x(neg),y(neg),'r.','markers',4); hold on; % plots positive values with blue, negative with red

                    xlim([x(1) x(end)]); % set x axis limits
                    ylim(obj.plotFreqs.ylimits); % set y axis limits

                    iStimuli = find(obj.plotFreqs.stimuli >= x(1) & obj.plotFreqs.stimuli <= x(end));
                    t = '';
                    for iPodnet = 1:numel(iStimuli) % plot colored stimuli and responses
                        podnet = obj.E.PsyData.P.data(iStimuli(iPodnet),obj.E.PsyData.P.sloupce.kategorie);
                        plot([obj.plotFreqs.stimuli(iStimuli(iPodnet)) obj.plotFreqs.stimuli(iStimuli(iPodnet))]', ylim', 'LineWidth',3,'Color', obj.plotFreqs.colors{podnet+1}); hold on;    
                        t = strcat(t,'{\color{', obj.plotFreqs.colors{podnet+1},'}',obj.E.PsyData.P.strings.podminka(podnet+1),'}, '); % title
                    end

                    responses = obj.plotFreqs.responses(obj.plotFreqs.responses >= x(1) & obj.plotFreqs.responses <= x(end));
                    plot([responses responses]', repmat(ylim,length(responses),1)','k', 'LineWidth',2); % plot responses
                    if obj.plotFreqs.plotGroup; 
                        title(strcat(num2str(obj.plotFreqs.iChannels(i)), '. channel (', t, ') ',names{obj.plotFreqs.iChannels(i)},'-',labels{obj.plotFreqs.iChannels(i)}));
                    else
                        title(strcat(num2str(obj.E.Hf(obj.plotFreqs.iFreqs(i))), ' Hz (', t, ')'));
                    end

                    %%%%%%%%%%%%%%%%%%% PHASE %%%%%%%%%%%%%%%%%%% 
                    subplot(numSubplot,2,i*2)
                    % phase values
                    if obj.plotFreqs.plotGroup; 
                        y = squeeze(obj.E.fphase(time,obj.plotFreqs.iChannels(i),obj.plotFreqs.iFreqs)); % channel phase values
                        display('correct, more')
                    else
                        y = squeeze(obj.E.fphase(time,obj.plotFreqs.iChannels,obj.plotFreqs.iFreqs(i))); % frequency phase values
                    end 
                    %phasemap('rad')

                    %imagesc([x(1) x(end)],[0,5],repmat(y',[5,1])); hold on;

                    plot(x,y, 'Color', 'b'); hold on;
                    xlim([x(1) x(end)]); % set x axis limits
                    ylim([-5,5]); % set y axis limits
                    iStimuli = find(obj.plotFreqs.stimuli >= x(1) & obj.plotFreqs.stimuli <= x(end));
                    t = '';
                    for iPodnet = 1:numel(iStimuli) % plot colored stimuli and responses
                        podnet = obj.E.PsyData.P.data(iStimuli(iPodnet),obj.E.PsyData.P.sloupce.kategorie);
                        plot([obj.plotFreqs.stimuli(iStimuli(iPodnet)) obj.plotFreqs.stimuli(iStimuli(iPodnet))]', ylim', 'LineWidth',3,'Color', obj.plotFreqs.colors{podnet+1}); hold on;    
                        t = strcat(t,'{\color{', obj.plotFreqs.colors{podnet+1},'}',obj.E.PsyData.P.strings.podminka(podnet+1),'}, '); % title
                    end

                    responses = obj.plotFreqs.responses(obj.plotFreqs.responses >= x(1) & obj.plotFreqs.responses <= x(end));
                    plot([responses responses]', repmat(ylim,length(responses),1)','k', 'LineWidth',2); % plot responses
                end
                
            end
           
        end
         
    end
    
    
    
end

