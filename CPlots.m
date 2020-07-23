classdef CPlots < matlab.mixin.Copyable
    %CPLOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Eh; %handle to CiEEGData object       
        PlotRChMean; %data for  PlotResponseChMean        
        plotTimeInt = struct; %stavove udaje o grafu TimeIntervals
    end
    methods (Access = public)
        function obj = CPlots(E) %constructor
            obj.Eh = E; %save handle to main object
        end 
        function PlotResponseChMean(obj,kategories,ylimits, labels, channels,fdrlevel)
            if ~exist('kategories','var') || isempty(kategories)     %stimulus categories to plot            
                if isfield(obj.Eh.plotRCh,'kategories') && ~isempty(obj.Eh.plotRCh.kategories) %by default use the categories selected in PlotResponseCh
                    kategories = obj.Eh.plotRCh.kategories;                           
                elseif ~isempty(obj.Eh.Wp) && isfield(obj.Eh.Wp(obj.Eh.WpActive), 'kats') %second is to use the categories from stat
                    kategories = obj.Eh.Wp(obj.Eh.WpActive).kats; 
                else
                    kategories = obj.Eh.PsyData.Categories(); %or all categoires
                end
            end
            if ~exist('ylimits','var')             %limits of the y axis     
                if isfield(obj.PlotRChMean,'ylimits')
                    ylimits = obj.PlotRChMean.ylimits;
                else
                    ylimits = [];    
                end
            end
            
            if ~exist('channels','var') || isempty(channels)  %selection of channels to plot
                channels = obj.Eh.CH.plotCh2D.chshow;                
                figuretitle = [obj.Eh.CH.plotCh2D.chshowstr ' chns: ' num2str(numel(channels))];                
            else
                figuretitle = ['channels: ' num2str(numel(channels))];                
            end
            
            if ~exist('labels','var') || isempty(labels) % to distinguish the brain labels in the plot - only for one category to be plotted
                if ~isfield(obj.PlotRChMean,'labels') || isempty(obj.PlotRChMean.labels)
                   labels = {'no'}; %default value, meaning to plot all labels together
                else
                   labels = obj.PlotRChMean.labels;               
                end
            elseif isnumeric(labels) && labels == -1
                labels = {'no'}; %default value, meaning to plot all labels together
            end
            if ~exist('fdrlevel','var') || isempty(fdrlevel) 
                fdrlevel = 1;  %less strict fdr correction, 2 is more strict                      
            end
            
            obj.PlotRChMean.channels = channels;
            obj.PlotRChMean.kategories = kategories;
            obj.PlotRChMean.ylimits = ylimits;
            obj.PlotRChMean.labels = labels;
            if ~isfield(obj.PlotRChMean,'plotband'), obj.PlotRChMean.plotband = 0; end %default ciplot possible to copy to CorelDraw
            assert(numel(kategories)==1 || numel(labels)==1,'cannot plot multiple labels for multiple categories');
            
%             chnshow = mat2str(channels(1:(min(20,numel(channels)))));
%             if numel(channels) > 20, chnshow = [chnshow ' ...']; end
                        
            
            if strcmp(labels{1},'nrj') %plot the labels saved in CM.CH.channelPlot.plotCh3D.brainlabelColors
                [labels,barvy]=obj.GetBrainlabelsSaved();
                figuretitle = [figuretitle ' - ' obj.Eh.PsyData.CategoryName(kategories(1))];
            elseif ~strcmp(labels{1},'no') %plotting individual labels
                barvy = distinguishable_colors(numel(labels));               
                figuretitle = [figuretitle ' - ' obj.Eh.PsyData.CategoryName(kategories(1))];
            else %plotting categories
                barvy = cell2mat(obj.Eh.colorskat(1:end)'); %from cell array of rgb to matrix
            end
            HSV = rgb2hsv(barvy);
            HSV(:,2) = 0.15; %decrease saturation for error bands
            colorsErrorBars = hsv2rgb(HSV); % barvycellfun(@(a) min(a+hue, 1), barvy, 'UniformOutput', false);

            katlinewidth = 2;
            if isfield(obj.PlotRChMean,'fh') && (verLessThan('matlab','9.0') || isvalid(obj.PlotRChMean.fh)) %isvalid je od verze 2016
                figure(obj.PlotRChMean.fh); %pouziju uz vytvoreny graf
                clf(obj.PlotRChMean.fh); %graf vycistim
            else
                obj.PlotRChMean.fh = figure('Name','PlotResponseChMean','CloseRequestFcn', @obj.tearDownFigCallbackPlotResponseChMean);
            end
                        
            T = linspace(obj.Eh.epochtime(1),obj.Eh.epochtime(2),size(obj.Eh.d,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            ymax = 0;
            CHM = nan(size(obj.Eh.d,1),numel(kategories),numel(channels)); %time x kategories x channels. For multiple labels, each has different number of channels. So only part is filled
            if ~isempty(ylimits)
                ymin = ylimits(1);
            else
                ymin = -0.1;
            end
            ploth = zeros(1,max(numel(kategories),numel(labels))); % handles of plots for legend, we expect 1 category or 1 label.
            legendStr = cell(size(ploth)); 
            for k = 1 : numel(kategories) %index 1-3 (nebo 4)
                katnum = kategories(k);               
                for l = 1:numel(labels)
                    if strcmp(labels{l},'no') % all brainlabels, plotting categories
                        chsel = channels; %all channels
                        colorkatk = [barvy(katnum+1,:) ; colorsErrorBars(katnum+1,:)]; %dve barvy, na caru a stderr plochu kolem                            
                        legendStr{k} = obj.Eh.PsyData.CategoryName(katnum); % array of names of categories for a legend
                    else % plotting individual brainlabels for one stimulus category
                        brailabels = lower({obj.Eh.CH.brainlabels.label});
                        chsel = intersect(find(contains(brailabels,labels{l})),channels);
                        colorkatk = [barvy(l,:) ; colorsErrorBars(l,:)]; %dve barvy, na caru a stderr plochu kolem  
                        legendStr{l} = [labels{l} ' x' num2str(numel(chsel))];
                    end
                    [katdata,~,RjEpCh] = obj.Eh.CategoryData(katnum,[],[],chsel);%katdata now time x channels x epochs                                
                    for ich = 1:numel(chsel)
                        CHM(:,k,ich) = mean(katdata(:,chsel(ich),~RjEpCh(ich,:)),3); %mean over valid epochs from katdata
                    end
                    M = mean(CHM(:,k,1:numel(chsel)),3);  %mean over channels, 
                    E = std(CHM(:,k,1:numel(chsel)),[],3)/sqrt(numel(chsel)); %std err of mean over channels          
                    ymax = max(ymax,max(M+E));
                    if obj.PlotRChMean.plotband
                        plotband(T, M, E, colorkatk(2,:)); % transparent but copying to CorelDraw does not work (its all black)
                    else
                        ciplot(M+E, M-E, T, colorkatk(2,:)); % not transparent 
                    end
                    hold on;
                    ploth(max(k,l)) = plot(T,M,'LineWidth',katlinewidth,'Color',colorkatk(1,:));  %store handle to use in legend                                          
                end
            end 
            if numel(kategories) == 2         
                paired = 1; %paired
                Wp = CStat.Wilcox2D(CHM(:,1,:), CHM(:,2,:),0,fdrlevel,'kat 1 vs 2',[],[],paired); %more strict dep method
                iWp = Wp <= 0.05;      % Wilcoxon signed rank         
                y = ymin*0.95;
                WpLims = [find(iWp,1) find(iWp,1,'last')];
                if ~isempty(WpLims)
                    plot([T(WpLims(1)) T(WpLims(2))],[y y],'LineWidth',5,'Color',[255 99 71]./255); %full line between start and end of significance
                    if mean(iWp(WpLims(1) : WpLims(2)))<1 %if there are non signifant parts 
                        y = ymin*0.85;
                        plot(T(iWp),ones(1,sum(iWp))*y, '*','Color','red'); %stars
                    end
                end
            end
            xlim(obj.Eh.epochtime(1:2)); 
            if ~isempty(ylimits), ylim(ylimits); end
            legend(ploth,legendStr,'Location','best');
            title(figuretitle);          
            xlabel('time [s]');
            ylabel('BGA power change');
            if ~isfield(obj.PlotRChMean,'filterListener') || isempty(obj.PlotRChMean.filterListener) %function to be calle when event CHHeader.FilterChanged happens
                obj.PlotRChMean.filterListener = addlistener(obj.Eh.CH, 'FilterChanged', @obj.filterChangedCallbackPlotResponseChMean);
            end
        end
        function TimeIntervals(obj,vch,intervals,ylimits, nofile,store) % Sofiia 2020
            %  average time intervals for one channel or for a vector of channels (e.g. across a particular structure)
            %  numbers of channels can be obtained, for example, after applying CM.CH.FilterChannels({'lobe','precun'})
            %  or if it's empty it will use numbers of channels after current filtering (from CM.CH.sortoder)
            %  default intervals = [0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8];
            %  nofile - if 1, do not save any xls file, only plot figure
            
            %process arguments
            if ~exist('vch','var') || isempty(vch) %vector of channels
                vch = obj.Eh.CH.sortorder; % will use numbers of channels after current filtering
                chshowstr = obj.Eh.CH.plotCh2D.chshowstr; %description of the selected channels - name from CM.CH.FilterChannels();
            else
                chshowstr = [ num2str(length(vch)) 'channels'];
            end 
            if ~exist('intervals','var') || isempty(intervals) 
                if ~isfield(obj.plotTimeInt,'intervals')                               
                    obj.plotTimeInt.intervals = [0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8]; %store the default intervals
                end
                intervals = obj.plotTimeInt.intervals;
            elseif any(intervals == -1) %-1 to reset default values
                obj.plotTimeInt.intervals = [0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8]; %store the default intervals
                intervals = obj.plotTimeInt.intervals;
            else
                obj.plotTimeInt.intervals = intervals; %store the limits given in function call
            end
            
            if ~exist('ylimits','var') || isempty(ylimits)  %limits for the y axis
                if ~isfield(obj.plotTimeInt,'ylimits')
                    obj.plotTimeInt.ylimits = []; %store the empty limits - the default value
                end            
            elseif any(ylimits == -1) %-1 to reset default values
                obj.plotTimeInt.ylimits = []; %store the empty limits - the default value
            else
                obj.plotTimeInt.ylimits = ylimits; %store the limits given in function call
            end 
            ylimits = obj.plotTimeInt.ylimits;
            if ~exist('nofile','var') || isempty(nofile)
                if ~isfield(obj.plotTimeInt,'nofile') %save the file by default
                     obj.plotTimeInt.nofile = 0; %the default value - save the file
                end
            else
                 obj.plotTimeInt.nofile = nofile; %store the limits given in function call
            end 
        
            if isfield(obj.Eh.plotRCh,'kategories') && ~isempty(obj.Eh.plotRCh.kategories)
                kats = obj.Eh.plotRCh.kategories;
            else
                kats = obj.Eh.Wp(obj.WpActive).kats; % define categories
            end
            
            % initialize matrix for all channels
            chanMeans = zeros (size(intervals,1), numel(kats), length(vch)); % intervals x kats x channels  -  for statistica - all time intervals and category in columns
            
            %get data from obj.CategoryData for all channels at onces - to make the whole function faster
            categorydata = cell(numel(kats),2); %store precomputed category data here 
            legendStr = cell(1,numel(kats)); % strings for legend and xls table
            for k = 1:numel(kats)
                [d,~,RjEpCh,~]= obj.Eh.CategoryData(kats(k),[],[],vch); % d:time x channels x epochs - get data for all channels
                categorydata(k,:) = {d(:,vch,:),RjEpCh}; % save d (time x channels x epochs) and RjEpCh (channels x epochs)
                legendStr{k} = obj.Eh.PsyData.CategoryName(kats(k)); % array of names of categories for a legend
            end
            
            for ch = 1:size(vch,2) % for each channel                
                % initialize matrix for an individual channel
                allmeans = cell(size(intervals,1),numel(kats)); %  time intervals x kategories (x epochs)
                categories = NaN(size(obj.Eh.d,3), 1); % vector to define categories for all epochs for an individual channel - for xls table
                T = (0 : 1/obj.Eh.fs : (size(obj.Eh.d,1)-1)/obj.Eh.fs)' + obj.Eh.epochtime(1); % time in sec
                nInd = 1; % initialize index for variable categories
                
                for k = 1:numel(kats) % for each category%                   
                    RjEpCh = categorydata{k,2}(ch,:); %previously stored RjEpCh
                    d = squeeze(categorydata{k,1}(:,ch,~RjEpCh));  %time x epochs, previously stored d value                  
                    categories(nInd:(nInd+size(d,2)-1)) = kats(k); % to establish the number of category for all epochs for an individual channel - for xls table
                    nInd = size(d,2)+nInd;
                    
                    for int = 1:size(intervals,1)  % for each interval
                        index = T>intervals(int,1) & T <=intervals(int,2);
                        meanOverT = mean(d(index, :),1); % mean across time  according to time intervals - for each epoch
                        allmeans{int,k} = meanOverT; % put in a cell array; time intervals x epochs
                    end
                    
                end
                
                % means and std err over epochs
                MOverEpoch = zeros(size(intervals,1),numel(kats)); % all means over epochs in one channel - for plot
                EOverEp = zeros(size(intervals,1),numel(kats)); % all std err over epochs in one channel - for plot
                for k = 1:numel(kats) % for each category
                    catmeans = cell2mat(allmeans(:,k)); % time invervals x epochs
                    MOverEpoch(:,k) = mean(catmeans,2); % mean over epochs
                    EOverEp(:,k) = std(catmeans,[],2)/sqrt(size(catmeans,2)); % std err of mean over epochs
                    chanMeans(:,k,ch) = MOverEpoch(:,k)'; % to put means over epochs for each category and channel - for xls table
                end
            end
                  
            %store or retrieve the plotted data
            if exist('store','var') && ~isempty(store) && isstruct(store)
                if isfield(store,'store') && store.store==1
                        if ~isfield(obj.plotTimeInt,'data'), obj.plotTimeInt.data = struct(); end
                        if ~isfield(store,'label') || isempty(store.label) %any label, e.g. brainlabel
                            store.label = 'NaN';
                        end
                        if isfield(obj.plotTimeInt.data,'label')                        
                            labels = {obj.plotTimeInt.data.label}; %labels stored until now
                            nlabel = find(contains(labels,store.label));
                            if isempty(nlabel)
                                nlabel = numel(obj.plotTimeInt.data)+1;
                            end
                        else                                
                            nlabel = 1;
                        end
                            
                        if ~isfield(store,'mark') || isempty(store.mark)  %any mark, e.g. name of the current marking (fghjkl)
                            store.mark = 'NaN';                            
                        end
                        if isfield(obj.plotTimeInt.data,'mark')                        
                            marks = {obj.plotTimeInt.data.mark}; %labels stored until now
                            nmark = find(contains(marks,store.mark));
                            if isempty(nmark)
                                nmark = numel(obj.plotTimeInt.data)+1;
                            end
                        else                                
                            nmark = 1;
                        end
                            
                        n = intersect(nmark,nlabel);
                        if isempty(n), n = numel(obj.plotTimeInt.data)+1; end
                        obj.plotTimeInt.data(n).label = store.label;
                        obj.plotTimeInt.data(n).mark = store.mark;
                        obj.plotTimeInt.data(n).chanMeans = chanMeans; % intervals x kats x channels
                        obj.plotTimeInt.data(n).chshowstr = chshowstr; %current filter label
                        obj.plotTimeInt.data(n).legendStr = legendStr; %current filter label
                        obj.plotTimeInt.data(n).chnum = numel(vch); %number of channels
                        obj.plotTimeInt.data(n).channels = vch; %channel numbers
                        obj.plotTimeInt.data(n).intervals = intervals; %number of channels
                        disp(['data stored as "' store.label ' & ' store.mark '"']);
                   
                elseif isfield(store,'kat') && store.kat >0 && store.kat <= size(obj.plotTimeInt.data(1).chanMeans,2)
                    %number of category to plot
                    if isfield(store,'mark')
                        marks = {obj.plotTimeInt.data.mark};
                        idata = contains(marks,store.mark); %logical index
                        idata1 = find(idata,1);
                    else
                        store.mark = 'NaN';
                        idata = true(numel(obj.plotTimeInt.data),1);
                        idata1 = 1;
                    end
                    assert(sum(idata)>0, 'no data to plot');
                    if ~isfield(store,'normalize') %normalize each label interval values to it maximal interval mean - to be able to compare the speed of response
                        store.normalize = 0;
                    elseif store.normalize == 1
                        ylimits = [-.1 1.1]; %the values will be [0 1], temporary change 
                    end
                    vch = 1:max([obj.plotTimeInt.data(idata).chnum]); %use the stored channels as channels
                    chanMeans = nan(size(obj.plotTimeInt.data(idata1).chanMeans,1), numel(obj.plotTimeInt.data(idata)) , max(vch)); %intervals x labels x channels
                    legendStr = cell(1,numel(obj.plotTimeInt.data(idata))); 
                    kats = 1:numel(obj.plotTimeInt.data(idata)); %mimic different stored labels as categories
                    if strcmp(store.colors, 'nrj') %if to use colors from CM.CH.channelPlot.plotCh3D.brainlabelColors 
                        [brainlabels,labels_barvy] = obj.GetBrainlabelsSaved();
                        barvy = zeros(numel(kats),3);
                    end                    
                    datan = find(idata); %rows if plotTimeInt.data
                    for k = kats %rows of plotTimeInt.data
                        chanMeans(:,k,1:obj.plotTimeInt.data(datan(k)).chnum) = obj.plotTimeInt.data(datan(k)).chanMeans(:,store.kat,:);
                        if store.normalize
                            katmax = max(mean(chanMeans(:,k,1:obj.plotTimeInt.data(datan(k)).chnum),3));
                            chanMeans(:,k,1:obj.plotTimeInt.data(datan(k)).chnum) = chanMeans(:,k,1:obj.plotTimeInt.data(datan(k)).chnum)./katmax;
                        end
                        legendStr{k} = [ obj.plotTimeInt.data(datan(k)).label ' x' num2str(obj.plotTimeInt.data(datan(k)).chnum)];
                        if strcmp(store.colors, 'nrj')
                            ilabel = contains(brainlabels,obj.plotTimeInt.data(datan(k)).label);
                            barvy(k,:) = labels_barvy(ilabel,:);
                        end
                    end
                    chshowstr = [obj.plotTimeInt.data(idata1).legendStr{store.kat} ': mark=' store.mark  ]; %assumes same category accross data fields                   
                    
                end
            else
                store = [];
            end
            
            
            % plotting
            if isfield(obj.plotTimeInt,'fh') && ishandle(obj.plotTimeInt.fh)
                figure(obj.plotTimeInt.fh); %pouziju uz vytvoreny graf
                clf(obj.plotTimeInt.fh); %graf vycistim
            else
                obj.plotTimeInt.fh = figure('Name','TimeIntervals','CloseRequestFcn', @obj.tearDownFigCallbackTimeIntervals);                 
            end                          
      
         
            % plot mean and std across epochs for one channel
            if length(vch)==1 && isempty(store) && (~isfield(obj.plotTimeInt,'nosinglechannel') || obj.plotTimeInt.nosinglechannel == 0)
                M=MOverEpoch; %set variable to be plotted - means
                E=EOverEp;    %set variable to be plotted - errorbars
                % the number of the channel in the title of plot
                if ~isempty(obj.Eh.CH.sortedby) %pokud jsou kanaly serazene jinak nez podle cisla kanalu
                    chstr = [ num2str(vch(ch)) '(' obj.Eh.CH.sortedby  num2str(obj.Eh.plotRCh.ch) ')' ];
                elseif numel(obj.Eh.CH.sortorder) < obj.Eh.channels %pokud jsou kanaly nejak vyfiltrovane pomoci obj.Eh.CH.FilterChannels();
                    chstr = [ num2str(vch(ch)) '(' num2str(obj.Eh.plotRCh.ch) '/'  num2str(numel(obj.Eh.CH.sortorder)) ')' ];
                else
                    chstr = num2str(vch(ch));
                end
                title(['channel ' chstr '/' num2str(obj.Eh.channels) ' - ' obj.Eh.PacientID()], 'Interpreter', 'none'); % v titulu obrazku bude i pacientID napriklad p132-VT18                
            else
                
                % plot mean and std across channels
                MOverChan = zeros(size(intervals,1),numel(kats)); % all means over channels - for plot
                EOverChan = zeros(size(intervals,1),numel(kats)); % all std err over channels - for plot
                for k = 1:numel(kats) % for each category
                    MOverChan(:,k) = nanmean(squeeze(chanMeans(:,k,:)),2); % mean over channels
                    EOverChan(:,k) = nanstd(squeeze(chanMeans(:,k,:)),[],2)/sqrt(size(chanMeans,3)); % std err of mean over channels
                end
                M=MOverChan; %set variable to be plotted - means
                E=EOverChan; %set variable to be plotted - errorbars
                title([num2str(length(vch)) ' channels: ' chshowstr ]); %title for multiple channels
            end            
            
            ploth = zeros(1,numel(kats)); % handles of plots for legend                        
            if numel(kats)>3                 
                if ~exist('barvy','var')                
                    barvy = distinguishable_colors(numel(kats));                
                end
            else
                barvy = cell2mat(obj.Eh.colorskat(2:end)');
            end
            for k=1:numel(kats)
                hold on
                errorbar(intervals(:, 2),M(:,k),E(:,k),'.','Color', barvy(kats(k),:)); % plot errors
                hold on;
                h_mean = plot(intervals(:, 2), M(:,k),'LineWidth',2,'Color',barvy(kats(k),:)); % plot means
                ploth(k) = h_mean; % plot handle                
            end
            if numel(kats)==2
                paired = 1; %paired
                fdrlevel= 1;
                if ~isempty(ylimits) 
                    y = ylimits(1)*0.8;
                else
                    y = 0;
                end
                Wp = CStat.Wilcox2D(chanMeans(:,1,:), chanMeans(:,2,:),0,fdrlevel,'kat 1 vs 2',[],[],paired); %more strict dep method
                iWp = Wp <= 0.05;
                plot(intervals(iWp, 2),ones(1,sum(iWp))*y, '*','Color','red'); %stars
            end
            legend(ploth,legendStr,'Location','best');
            xlim([intervals(1,2)-.05, max(max(intervals))+.05]);
            xticks([intervals(1, 1) (intervals( : , 2))'])
            txtinterv = num2str(intervals,'%.1f-%.1f');
            xticklabels({'0', txtinterv(1:end, :)})            
            xlabel('time intervals, s');
            if ~isempty(ylimits), ylim(ylimits); end
            if ~isfield(obj.plotTimeInt,'filterListener') || isempty(obj.plotTimeInt.filterListener) %function to be calle when event CHHeader.FilterChanged happens
                obj.plotTimeInt.filterListener = addlistener(obj.Eh.CH, 'FilterChanged', @obj.filterChangedCallbackTimeIntervals);
            end
                        
            % export data in xls table
            if obj.plotTimeInt.nofile==0 && isempty(store) 
                obj.TimeIntervals2XLS(vch,kats , legendStr, intervals, chanMeans,chshowstr, categories, allmeans);
            end
        end
        function TimeIntervalsColect(obj,marks,labels)
            %collects data using TimeIntervals from multiple filter settings
            %labels is cell array of brain labels, filtered by label = ''
            %marks is cell array of 1. marking sets 2. markings in two columns 
            %the intervals and ylimis needs to be set before
            if ~exist('labels','var') || isempty(labels) 
                labels = {obj.Eh.CH.channelPlot.plotCh3D.brainlabelColors.label};
            elseif ischar(labels) && strcmp(labels,'clear')
                obj.plotTimeInt.data = struct; %clear all stored data
                disp('all TimeIntervals data cleared');
                return;
            end
            if exist('marks','var')
                obj.plotTimeInt.nosinglechannel = 1; %to not use the special code for single channel plot
                for l = 1:numel(labels)
                    for m = 1:size(marks,1)
                        obj.Eh.SetSelChActive(marks{m,1});
                        klavesy = 'fghjkl';                 
                        iklavesy = ismember(klavesy,marks{m,2});
                        if sum(iklavesy)==1
                            markname = obj.Eh.plotRCh.selChNames{iklavesy};
                        else
                            markname = [];
                        end
                        obj.Eh.CH.FilterChannels({'label',labels{l}},{},[marks{m,2} 'n']);
                        obj.TimeIntervals([],[],[],1,struct('store',1,'label',labels{l},'mark',markname));
                    end
                end            
            end
        end
        function TimeIntervals2XLS(obj,varargin)
            if nargin > 2 %when called directly from TimeIntervals()
               %vch,kats , legendStr, intervals, chanMeans,chshowstr, categories, allmeans
               vch = varargin{1}; 
               kats = varargin{2}; 
               legendStr = varargin{3};
               intervals = varargin{4};
               chanMeans = varargin{5};
               chshowstr = varargin{6};
               categories = varargin{7};
               allmeans = varargin{8};
               if length(vch)==1  && (~isfield(obj.plotTimeInt,'nosinglechannel') || obj.plotTimeInt.nosinglechannel == 0)
                    % export data for one channel in xls table
                    categories(isnan(categories))=[];  % remove all NaN
                    tableX = num2cell([categories, (cell2mat(allmeans))']); % table for export, epochs over all categories x (kategory num, intervals)
                    intervStr = cellstr(num2str(intervals,'%.1f-%.1f'));
                    titles4table = ['category', intervStr']; % column names
                    tableX = [titles4table; tableX]; % final table
                    strch = ['channel-' num2str(vch)];  % to distinguish in filename more just one channel
               else
                   strch = ['channels-' regexprep(chshowstr,{'{','}','[',']',' ','&'},{'','','','','-',''})]; % to distinguish in filename more than one channel
                   [tableX,titles4table] = obj.TimeIntervals2XLSTable(vch,chanMeans,kats,legendStr,intervals);
               end
            else
                if nargin == 2
                    normalize = varargin{1}; %if to normalize the values
                else
                    normalize = 0; %default value
                end
                row = 1;
                for n = 1:numel(obj.plotTimeInt.data)
                    vch = obj.plotTimeInt.data(n).channels;
                    chanMeans = obj.plotTimeInt.data(n).chanMeans; %intervals x kats x channels
                    if normalize
                       katmax = max(mean(chanMeans(:,:,1:obj.plotTimeInt.data(n).chnum),3),[],1);  %mean over channel, max over intervals
                       chanMeans = chanMeans ./ katmax;
                    end
                    if n == 1
                        kats = 1:numel(obj.plotTimeInt.data(n).legendStr);
                        legendStr = (obj.plotTimeInt.data(n).legendStr);
                        intervals = (obj.plotTimeInt.data(n).intervals);
                        [tbl,titles4table]= obj.TimeIntervals2XLSTable(vch, chanMeans,kats,legendStr,intervals);
                        titles4table = [{'n','label','mark'},titles4table]; %#ok<AGROW>
                        tableX = cell(sum([obj.plotTimeInt.data.chnum]),numel(titles4table));  %empty table for all channels                      
                    else
                        [tbl]= obj.TimeIntervals2XLSTable(vch, chanMeans);         
                    end
                    tableX(row:row+numel(vch)-1,:) = [num2cell(row:row+numel(vch)-1)', ...
                            repmat({obj.plotTimeInt.data(n).label},size(tbl,1),1), ...
                            repmat({obj.plotTimeInt.data(n).mark},size(tbl,1),1), ...
                            tbl ...
                            ];
                    row = row + numel(vch);
                end
                strch = ['Data' iff(normalize,'Norm','Orig')];
            end
    
            xlsfilename = ['logs\TimeIntervals_' strch '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.xls'];
            xlswrite(xlsfilename,[titles4table; tableX]);
            disp([ 'xls tables saved: ' xlsfilename]);
           
        end
    end
    methods  (Access = private)
        function [tableX,titles4table]= TimeIntervals2XLSTable(obj,vch, chanMeans,kats , legendStr, intervals)                   
            % export data for all channels (means over epochs and time intervals) in xls table   
            if ~isempty(obj.Eh.CH.brainlabels) % if CM.CH.brainlabels contains names for brain labels, they will be put in a final table according to the number of channel
                labels = horzcat({obj.Eh.CH.brainlabels(vch).class}',{obj.Eh.CH.brainlabels(vch).label}', {obj.Eh.CH.brainlabels(vch).lobe}');
            else
                labels = cell(length(vch),3);% otherwise labels will be just empty cell array
            end
            names = {obj.Eh.CH.H.channels(vch).name}'; %channel names
            pacients = obj.Eh.CH.PacientTag(vch); %pacient tags
            chanMeans = permute(chanMeans,[3 1 2]); %make channels the first dim
            chanMeans = reshape(chanMeans,size(chanMeans,1), size(chanMeans,2)*size(chanMeans,3)); %concat the 2. and 3. dim
            tableX = [num2cell(vch'), names, pacients, labels, num2cell(chanMeans)]; % final table
            
            if exist('kats','var')  %if to return also table header
                % long names for columns in STATISTICA (category + time interval)
                intervStr = cellstr(num2str(intervals,'%.1f-%.1f'));
                names4col = cell(1,numel(kats)*size(intervals,1));
                col=1;
                for k=1:numel(kats)
                    for int = 1:size(intervals,1)
                        Str = {legendStr{k},intervStr{int}};
                        names4col(:,col) = join(Str, ': ');
                        col=col+1;
                    end
                end
                titles4table = ['channel','channelname', 'pacient', 'class','label','lobe',names4col]; % column names 
            end
        end
        function [labels,barvy]=GetBrainlabelsSaved(obj)
            order = [obj.Eh.CH.channelPlot.plotCh3D.brainlabelColors.n];
            [~,idx] = sort(order);
            labels = {obj.Eh.CH.channelPlot.plotCh3D.brainlabelColors(idx).label}; %in the correct order
            barvy = vertcat(obj.Eh.CH.channelPlot.plotCh3D.brainlabelColors(idx).color);
        end
        
        function filterChangedCallbackPlotResponseChMean(obj,~,~)   %update chart if the filter is changed         
            if isfield(obj.PlotRChMean,'fh') && ishandle(obj.PlotRChMean.fh)
                %only if the PlotRChMean plot is shown
                if ~isequal(obj.Eh.CH.plotCh2D.chshow, obj.PlotRChMean.channels) %jinak se to z nejakeho duvodu po zmene filtru vola porad dokolecka
                    obj.PlotResponseChMean(obj.PlotRChMean.kategories);
                end
            end
        end
        function filterChangedCallbackTimeIntervals(obj,src,~)   %update chart if the filter is changed                                
            %src is the source object of CHHeader
            if isfield(obj.plotTimeInt,'fh') && ishandle(obj.plotTimeInt.fh)
                %only if the TimeIntervals plot is shown
                if ~isfield(obj.plotTimeInt,'chshowstr') || ~isequal(obj.plotTimeInt.chshowstr, src.plotCh2D.chshowstr)
                    %if the filter stored is different from the new filter
                    obj.TimeIntervals();
                    obj.plotTimeInt.chshowstr = src.plotCh2D.chshowstr;
                end
            end
        end
        function tearDownFigCallbackPlotResponseChMean(obj,src,~)
            delete(obj.PlotRChMean.filterListener);
            delete(src);
        end  
        function tearDownFigCallbackTimeIntervals(obj,src,~)
            delete(obj.plotTimeInt.filterListener);
            delete(src);
        end 
    end
    methods (Static,Access = public)       
        function PlotElectrodeEpiEvents(els,RjCh,DE,objtabs,tabs_orig,epochs,samples,epochtime,timeax,timeaxy,sec,time_n,elmaxmax,shift,iD)
            %27.4. - nahrazuju cast funkce CiEEGData.PlotElectrode toutu funkci, kvuli zkraceni, ale predavam strasne moc argumentu
            %treba bych mohl uz driv predat hodne z nich? A mit je jako properties tehle tridy?
            epieventsum = 0;  
            weights = [];
            DE.Clear_iDEtabs();
            for ch = els
                if ~ismember(ch,RjCh)
                    if epochs <= 1 %neepochovana data
                        [epitime, weight] = DE.GetEvents( [objtabs(iD(1)) objtabs(iD(2))],ch,tabs_orig(1)); 
                    else
                        epochy = sec : min(sec+ceil(time_n/samples)-1 , epochs ); %cisla zobrazenych epoch
                        tabs = [ objtabs(1,epochy)' objtabs(end,epochy)' ]; %zacatky a konce zobrazenych epoch
                        [epitime, weight] = DE.GetEvents( tabs,ch,tabs_orig(1)) ;                            
                        if numel(epitime) > 0
                            epitime(:,1) = epitime(:,1) + epochtime(1);
                        end
                    end
                    if size(epitime,1) > 0
                        plot(epitime(:,1),shift(elmaxmax-ch+1,1),'o', 'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',6);
                        epieventsum = epieventsum + size(epitime,1);
                        weights = [weights; round(weight*100)/100]; %#ok<AGROW>
                    end
                end
            end             
            text(timeax(1)+((timeax(end)-timeax(1))/8),timeaxy,['epileptic events: ' num2str(epieventsum) ' (' num2str(min(weights)) '-' num2str(max(weights)) ')']);         
        end
        function spatne = EpochsEpi(RjEpochCh,objels,CHH)
            figure('Name','Rejected epochs in individual channels');
            pocty = sum(RjEpochCh,2)/size(RjEpochCh,2); %pocty epoch v jednotlivych kanalech / celkovym poctem epoch
            plot(pocty,'.-');
            for el = 1:numel(objels)-1
                line([objels(el) objels(el)]+1,[0 1],'Color',[0.5 0.5 0.5]);
            end
            line([1 CHH.H.selCh_H(end)],[0.30 0.30],'Color','red');
            spatne = find(pocty >= 0.30); 
            disp(['kanaly s podilem vyrazenych epoch >= 0.3:' mat2str(spatne)]);                
        end
        function PlotEpochsKatPPA(E,ch,s)
            %udaje z PPAconditionsEEG_info.xlsx
            %unikatni cisla 325 obrazku, tak jak sly za sebou
            %ch - kanal k vykresleni
            %s - sample k vykresleni (okamzik v case)
            jpgnum = [168 46 121 231 303 4 102 178 203 283 70 181 149 199 317 99 2 76 44 226 210 175 81 64 80 165 143 47 247 126 240 190 185 180 318 306 251 26 272 96 69 273 59 278 195 188 289 218 184 41 74 48 161 3 63 205 182 235 8 100 141 248 298 255 157 284 109 106 40 145 65 164 67 130 51 228 103 319 111 155 139 287 73 146 34 62 5 253 68 75 214 71 325 232 250 20 204 276 224 305 91 201 233 55 163 200 60 260 77 294 206 215 245 236 266 194 187 11 85 110 220 54 167 186 316 312 118 137 95 135 323 144 254 264 244 92 241 58 173 208 45 324 134 290 90 9 209 160 97 222 18 14 219 300 93 280 243 38 36 176 52 17 147 56 309 170 191 35 119 33 249 216 279 322 292 72 258 89 275 133 50 117 108 107 13 189 152 39 193 177 263 285 22 88 87 24 66 313 131 304 230 301 153 297 259 125 286 113 239 140 198 267 227 308 257 84 223 183 269 171 94 83 132 86 158 142 211 150 154 129 21 270 237 42 242 25 234 28 217 15 197 23 192 256 202 29 315 225 123 1 101 311 265 321 16 114 7 179 57 151 53 27 196 271 320 282 159 31 274 37 212 172 299 262 19 207 112 120 238 295 82 116 12 156 246 302 79 61 261 104 281 30 138 314 32 148 310 122 124 136 115 128 268 213 98 162 221 291 293 127 288 169 49 6 252 229 10 277 166 78 174 307 105 43 296 ];
            
            %unikatni cisla umelych scen, v nespravnem poradi
            artificial = [231 226 247 240 251 272 273 235 248 255 228 253 232 250 233 260 245 236 266 254 264 244 241 243 249 258 275 263 230 259 239 267 227 257 269 270 237 242 234 256 265 271 274 262 238 246 261 268 252 229 ];
            iArt = ismember(jpgnum,artificial); %indexy v jpgnum, kde jsou cisla z artificial
            %unikatni cisla prirozenych scen
            natural = [303 283 317 318 306 278 289 298 284 319 287 325 276 305 294 316 312 323 324 290 300 280 309 279 322 292 285 313 304 301 297 286 308 315 311 321 320 282 299 295 302 281 314 310 291 293 288 277 307 296  ];
            iNat = ismember(jpgnum,natural);
            
            figure;
            d = squeeze(E.d(s,ch,1:325)); %jen prvnich 325 epoch
            x = 1:325; %fixne potech epoch
            plot(x,d);
            hold on;
            plot(x(iNat),d(iNat),'og');
            nmean = mean(d(iNat));
            nerr = std(d(iNat))/sqrt(sum(iNat));
            plot(x(iArt),d(iArt),'*r');
            amean = mean(d(iArt));
            aerr = std(d(iArt))/sqrt(sum(iArt));
            line([1 325],[nmean nmean],'Color',[0 1 0]);
            line([1 325],[nmean+nerr nmean+nerr],'Color',[0 0.5 0]);
            line([1 325],[amean amean],'Color',[1 0 0]);
            line([1 325],[amean+aerr amean+aerr],'Color',[0.5 0 0]);
            ylim([-2 2]);
            fprintf('Art: %f +- %f, Nat %f +- %f\n',amean,aerr,nmean,nerr);
        end

    end
    
end

