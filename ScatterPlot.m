classdef ScatterPlot < handle
    %SCATTERPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ieegdata
        dispChannels
        dispData; %vyobrazena data ve scatterplotu
        
        selCh ; %kopie ieegdata.plotRCh.selCh;
        selChNames ; %kopie ieegdata.plotRCh.selChNames;
        dispSelCh
        dispSelChName
        
        dispFilterCh
        
        connectPairs
        pairsPlot
        
        fig
        ax
        plots
        sbox
        pbox
        
        showNumbers; %zobrazit popisky hodnot, 0=nic 1=cisla 2=jmena kanalu
        numbers
        markerSize
        
        axisX
        axisY
        
        valFraction
        intFraction
        
        categories %seznam cisel kategorii od nejdulezitejsi (podle poradi ve statistice)
        categoryNames %jmena kategorii odpovidajici obj.categories
        categoriesSelectionIndex %1:numel(obj.categories)
        
        filterListener
    end
    
    methods
        function obj = ScatterPlot(ieegdata)
            %SCATTERPLOT Construct an instance of this class
            %   Detailed explanation goes here
            obj.ieegdata = ieegdata;
            obj.selCh = ieegdata.plotRCh.selCh; %vyber kanalu fghjkl * pocet kanalu
            obj.selChNames = ieegdata.plotRCh.selChNames; %vyber kanalu fghjkl * pocet kanalu
            obj.dispSelChName = [];
            obj.dispSelCh = 1:size(obj.selCh,1);  % Zobrazuji vse
            
            obj.dispFilterCh = obj.ieegdata.CH.sortorder; % Vyber podle FilterChannels
            obj.connectPairs = false;
            obj.showNumbers = 0; %defaultne se zadne labels nezobrazuji
            obj.numbers = [];
            obj.markerSize = 34;
            
            obj.setCategories();
            
            obj.drawScatterPlot(0.5, 0.5, 'tmax', 'valmax'); %TODO; zmenit na tmax, valmax, 0.5, 0.5
        end
        
        function drawScatterPlot(obj, valFraction, intFraction, axisX, axisY)
            obj.setTriggerValues(valFraction, intFraction);
            obj.initAxes(axisX, axisY);
            obj.updatePlot();
            obj.fixAxesLimits();
        end
        
        function PlotBrain(obj,katnum,xy,rangeZ)
            if ~exist('katnum','var') || isempty(katnum), katnum = 1; end
            if ~exist('xy','var'), xy = 'y'; end %defaultne osaY = valmax napriklad
            if ~exist('rangeZ','var') %pokud neni zadane
                rangeZ = iff(xy=='x',xlim(obj.ax),ylim(obj.ax));  %nastavi se podle limitu scatterplotu
            end; 
            selChFiltered = obj.selCh(obj.dispChannels,:); %chci zobrazovat jen signif odpovedi
            iData = logical(selChFiltered(:,katnum)); %kanaly se signifikantim rozdilem vuci baseline v teto kategorii
            if(xy=='x')
                data = obj.dispData(katnum).dataX(iData);
                dataName = obj.axisX;
            else
                data = obj.dispData(katnum).dataY(iData);
                dataName = obj.axisY;
            end
            obj.ieegdata.CH.plotCh3D.selch = []; %nechci mit vybrany zadny kanal z minula
            obj.ieegdata.CH.ChannelPlot([],0,data,... %param chnvals
                obj.dispChannels(iData),... %chnsel jsou cisla kanalu, pokud chci jen jejich vyber
                [],[],[obj.categoryNames{katnum} ':' dataName ', show:' obj.ieegdata.CH.plotCh2D.chshowstr], ... %popis grafu = title
                rangeZ); %rozsah hodnot - meritko barevne skaly
            %set(obj.plotAUC.Eh.CH.plotCh3D.fh, 'WindowButtonDownFcn', {@obj.hybejPlot3Dclick, selch});
        end
        
        function setXYLim(obj,xrange,yrange) %set axes limit
            if ~exist('xrange','var') || isempty(xrange), xrange = xlim(obj.ax); end
            if ~exist('yrange','var') || isempty(yrange), yrange = ylim(obj.ax); end
            xlim(xrange);
            ylim(yrange);
        end

    end

    methods(Access = private)
        
        function setCategories(obj)
            if ~isempty(obj.ieegdata.Wp) && isfield(obj.ieegdata.Wp(obj.ieegdata.WpActive), 'kats') %prvni volba je pouzit kategorie ze statistiky
                obj.categories  = flip(obj.ieegdata.Wp(obj.ieegdata.WpActive).kats);
            else
                obj.categories = obj.ieegdata.PsyData.Categories(); %pokud nejsou, pouziju vsechny kategorie
            end
            obj.categoryNames = strings(size(obj.categories));
            for k = 1 : numel(obj.categories)                
                catnum = cellval(obj.categories,k);%cislo kategorie, muze byt cell, pokud vice kategorii proti jedne
                obj.categoryNames(k) = obj.ieegdata.PsyData.CategoryName(catnum);
            end
            obj.categoriesSelectionIndex = 1:numel(obj.categories);
        end
 
        function initAxes(obj, axisX, axisY)
           switch axisX
                case 'valmax'
                    labelX = 'v_{max}';
                case 'tmax'
                    labelX = 't_{max}';
                case 'tfrac'
                    labelX = ['t_{' num2str(obj.valFraction) '}'];
                case 'tint'
                    labelX =  ['t_{int, ' num2str(obj.intFraction) '}'];
                otherwise
                    disp('X axis specification must be one of: valmax, tmax, tfrac, tint');
                    return;
            end
            
            switch axisY
                case 'valmax'
                    labelY = 'v_{max}';
                case 'tmax'
                    labelY = 't_{max}';
                case 'tfrac'
                    labelY = ['t_{' num2str(obj.valFraction) '}'];
                case 'tint'
                    labelY =  ['t_{int, ' num2str(obj.intFraction) '}'];
                otherwise
                    disp('Y axis specification must be one of: valmax, tmax, tfrac, tint');
                    return;
            end
            
            obj.axisX = axisX;
            obj.axisY = axisY;
            
            obj.fig = figure('CloseRequestFcn',@obj.tearDownFigCallback,'Name','ScatterPlot');
            
            obj.filterListener = addlistener(obj.ieegdata.CH, 'FilterChanged', @obj.filterChangedCallback);
            
            obj.ax = axes(obj.fig);
            xlabel(obj.ax, labelX);
            ylabel(obj.ax, labelY);
            set(obj.fig,'KeyPressFcn',@obj.hybejScatterPlot);
        end
        
        function setTriggerValues(obj, valFraction, intFraction)
            obj.valFraction = valFraction;
            obj.intFraction = intFraction;
        end
        
        function fixAxesLimits(obj)
            xlim(xlim(obj.ax));
            ylim(ylim(obj.ax));
        end
        
        function updatePlot(obj)
            obj.setDisplayedChannels(); % Kombinace voleb pro zobrazeni kanalu
            selChFiltered = obj.selCh(obj.dispChannels,:); %filter kanalu ve vyberu fghjkl
            delete(obj.plots);
            obj.plots = [];
            delete(obj.sbox);
            delete(obj.pbox);            
            delete(obj.pairsPlot); obj.pairsPlot = [];
            delete(obj.numbers); obj.numbers = [];            
            
            if ~isempty(obj.dispSelChName)
                obj.sbox = annotation(obj.fig, 'textbox',[0 .9 .4 .1], 'String', obj.dispSelChName, 'EdgeColor', 'none');
            end
            
            catlist = strjoin(obj.categoryNames(obj.categoriesSelectionIndex), ', ');
            obj.pbox = annotation(obj.fig, 'textbox', [0 0 .4 .1], 'String', ['C: ' catlist], 'EdgeColor', 'none');
            
            stats = struct();
            if isempty(obj.dispChannels)
                disp('No channels corresponding to the selection');
                return;
            end
            
            for k = obj.categoriesSelectionIndex %1:numel(obj.categories)
                catnum = obj.categories(k);
                [stats(k).valmax, stats(k).tmax, stats(k).tfrac, stats(k).tint] = obj.ieegdata.ResponseTriggerTime(obj.valFraction, obj.intFraction, catnum, obj.dispChannels);
            end
            
            categoryMarkers = {'o', 's', 'd','x'};
            
            hold(obj.ax, 'on');
            legend(obj.ax, 'off');
            
            if obj.connectPairs     % Nakresli linku spojujici prislusny par. Ruzne barvy musi byt samostatny plot (aby mohl scatter zustat ve stejnych osach)
                if length(obj.categoriesSelectionIndex) > 1
                    catIndex = zeros(size(obj.categoriesSelectionIndex(1)));
                    x = zeros(length(obj.categoriesSelectionIndex(1)), length(stats(1).(obj.axisX)));
                    y = zeros(length(obj.categoriesSelectionIndex(1)), length(stats(1).(obj.axisX)));
                    for cat = 1:length(obj.categoriesSelectionIndex)
                        catIndex(cat) = obj.categoriesSelectionIndex(cat);
                        x(cat,:) = stats(cat).(obj.axisX);
                        y(cat,:) = stats(cat).(obj.axisY);
                    end
                    l = length(stats(catIndex(1)).(obj.axisX));
                    for c1 = 1:length(obj.categoriesSelectionIndex)
                        for c2 = 1:c1-1
                            for k = 1:l
                                obj.pairsPlot(end+1) = plot([x(c1,k) x(c2,k)], [y(c1,k) y(c2,k)], 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
                            end
                        end
                    end
                else
                    disp('No categories to connect');
                    obj.connectPairs = false;
                end
            end
            
            baseColors = [0 1 0; 0 0 1; 1 0 0; 1 1 0; 1 0 1; 0 1 0];
            for k = obj.categoriesSelectionIndex %1:numel(obj.categories)
                dataX = stats(k).(obj.axisX);
                dataY = stats(k).(obj.axisY);
                iData = logical(selChFiltered(:,k)); %kanaly se signifikantim rozdilem vuci baseline v teto kategorii
                if any(iData) 
                    obj.plots(k,1) = scatter(obj.ax, dataX(iData), dataY(iData), obj.markerSize, repmat(baseColors(k,:), sum(iData), 1), categoryMarkers{k}, 'MarkerFaceColor', 'flat', 'DisplayName', obj.categoryNames{k});
                end
                iData = ~iData; %kanaly bez signif rozdilu vuci baseline v teto kategorii
                if any(iData) 
                    obj.plots(k,2) = scatter(obj.ax, dataX(iData), dataY(iData), obj.markerSize, repmat(baseColors(k,:), sum(iData), 1), categoryMarkers{k}, 'MarkerFaceColor', 'none', 'DisplayName', obj.categoryNames{k},...
                        'HandleVisibility','off'); %nebude v legende
                end
                if obj.showNumbers > 0
                    if obj.showNumbers == 1
                        labels = cellstr(num2str(obj.dispChannels')); %cisla kanalu
                    else
                        labels = {obj.ieegdata.CH.H.channels(obj.dispChannels).name}'; %jmena kanalu
                    end
                    dx = diff(xlim)/100;
                    th = text(dataX+dx, dataY, labels, 'FontSize', 8);
                    set(th, 'Clipping', 'on');
                    obj.numbers = [obj.numbers th];
                end
                obj.dispData(k).dataX = dataX; %zalohuju vyobrazena data pro jine pouziti
                obj.dispData(k).dataY = dataY;
            end
            legend(obj.ax, 'show');
            hold(obj.ax, 'off'); 
            title(['show:' obj.ieegdata.CH.plotCh2D.chshowstr]);
        end
        
        function hybejScatterPlot(obj,~,eventDat)
            switch eventDat.Key
                case {'u','i','o','p'}
                    ik = find('uiop'==eventDat.Key); % index 1-4
                    if ik <= numel(obj.categories)
                        if(ismember(ik, obj.categoriesSelectionIndex))
                            obj.categoriesSelectionIndex = setdiff(obj.categoriesSelectionIndex, ik);
                        else
                            obj.categoriesSelectionIndex = union(obj.categoriesSelectionIndex, ik);
                        end
                    end    
                    obj.updatePlot();
                case {'a'}  % Resset do zakladniho vyberu
                    obj.dispFilterCh = obj.ieegdata.CH.sortorder; % Vyber podle FilterChannels
                    obj.dispSelCh = 1:size(obj.selCh,1);   % Zrusit vyber dle SelCh
                    obj.dispSelChName = [];
                    obj.updatePlot();
                case {'f','g','h','j','k','l'}
                    obj.dispSelCh = find(obj.selCh(:,'fghjkl'==eventDat.Key)');
                    obj.dispSelChName = obj.selChNames{'fghjkl'==eventDat.Key};
                    obj.updatePlot();
                case {'s'}
                    obj.connectPairs = ~obj.connectPairs;
                    obj.updatePlot();
                case {'n'}
                    obj.showNumbers = iff(obj.showNumbers==0, 1, iff(obj.showNumbers==1,2,0)); %0->1->2->0
                    obj.updatePlot();
                case {'add'}
                    obj.markerSize = obj.markerSize + 8;
                    obj.updatePlot();
                case {'subtract'}
                    obj.markerSize = max(2, obj.markerSize - 8);
                    obj.updatePlot();
            end
        end
        
        function setDisplayedChannels(obj)
            obj.dispChannels = intersect(obj.dispFilterCh, obj.dispSelCh);
        end
        
        function filterChangedCallback(obj,~,~)
            obj.dispFilterCh = obj.ieegdata.CH.sortorder;   % Zmena vyberu dle filtru
            obj.updatePlot();
        end
        
        function tearDownFigCallback(obj,src,~)
            delete(obj.filterListener);
            delete(src)
        end      
        
    end
    
end

