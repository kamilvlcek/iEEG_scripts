classdef ScatterPlot < handle
    %SCATTERPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ieegdata
        dispChannels
        
        selCh
        selChNames
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
        
        showNumbers
        numbers
        markerSize
        
        axisX
        axisY
        
        valFraction
        intFraction
        
        categories
        categoryNames
        categoriesSelectionIndex
        
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
            obj.showNumbers = false;
            obj.numbers = [];
            obj.markerSize = 34;
            
            obj.setCategories();
            
            obj.drawScatterPlot(0.5, 0.5, 'tmax', 'valmax');
        end
        
        function drawScatterPlot(obj, valFraction, intFraction, axisX, axisY)
            obj.setTriggerValues(valFraction, intFraction);
            obj.initAxes(axisX, axisY);
            obj.updatePlot();
            obj.fixAxesLimits();
        end

    end

    methods(Access = private)
        
        function setCategories(obj)
            obj.categories = obj.ieegdata.PsyData.Categories();
            obj.categoryNames = strings(size(obj.categories));
            for k = 1 : numel(obj.categories)
                if iscell(obj.categories) %iff tady nefunguje, to by bylo samozrejme lepsi85858
                    catnum = obj.categories{k}; %cislo kategorie, muze byt cell, pokud vice kategorii proti jedne
                else
                    catnum = obj.categories(k); %cislo kategorie, muze byt cell, pokud vice kategorii proti jedne
                end
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
            
            obj.fig = figure('CloseRequestFcn',@obj.tearDownFigCallback);
            
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
            
            delete(obj.plots);
            obj.plots = [];
            
            delete(obj.sbox);
            if ~isempty(obj.dispSelChName)
                obj.sbox = annotation(obj.fig, 'textbox',[0 .9 .4 .1], 'String', obj.dispSelChName, 'EdgeColor', 'none');
            end
            
            delete(obj.pbox);
            catlist = strjoin(obj.categoryNames(obj.categoriesSelectionIndex), ', ');
            obj.pbox = annotation(obj.fig, 'textbox', [0 0 .4 .1], 'String', ['C: ' catlist], 'EdgeColor', 'none');
            
            stats = struct();
            for k = obj.categoriesSelectionIndex
                catnum = obj.categories(k);
                [stats(k).valmax, stats(k).tmax, stats(k).tfrac, stats(k).tint] = obj.ieegdata.ResponseTriggerTime(obj.valFraction, obj.intFraction, catnum, obj.dispChannels);
            end
            
            categoryMarkers = {'x', 'o', 's', 'd'};
            colors = obj.dispChannels/max(obj.dispChannels);    % Normalizace na 0,1
            
            hold(obj.ax, 'on');
            delete(obj.pairsPlot); obj.pairsPlot = [];
            if obj.connectPairs     % Nakresli linku spojujici prislusny par. Ruzne barvy musi byt samostatny plot (aby mohl scatter zustat ve stejnych osach)
                legend(obj.ax, 'off');  % V tomhle pripade by mela kazda linka vlastni legendu
                if length(obj.categoriesSelectionIndex) == 2
                    k1 = obj.categoriesSelectionIndex(1);
                    k2 = obj.categoriesSelectionIndex(2);
                    l = length(stats(k1).(obj.axisX));
                    x1 = stats(k1).(obj.axisX);
                    x2 = stats(k2).(obj.axisX);
                    y1 = stats(k1).(obj.axisY);
                    y2 = stats(k2).(obj.axisY);
                    
                    % Mapovani barev
                    cmap = colormap(obj.ax);
                    colorIndex = floor(colors * size(cmap,1));
                    colorIndex(colorIndex == 0) = 1;
                    colorIndex(colorIndex > size(cmap,1)) = size(cmap,1);
                    
                    for k = 1:l-1
                        obj.pairsPlot(k) = plot([x1(k) x2(k)], [y1(k) y2(k)], 'Color', cmap(colorIndex(k), :));
                    end
                else
                    disp('Pair connection must include exactly 2 categories');
                    obj.connectPairs = false;
                end
            else
                legend(obj.ax);  % Zobrazit legendu, pokud nejsou zobrazene linky
            end
            
            delete(obj.numbers); obj.numbers = [];
            for k = obj.categoriesSelectionIndex
                dataX = stats(k).(obj.axisX);
                dataY = stats(k).(obj.axisY);
                obj.plots(k) = scatter(obj.ax, dataX, dataY, obj.markerSize, colors, categoryMarkers{k}, 'MarkerFaceColor', 'flat', 'DisplayName', obj.categoryNames(k));
                if obj.showNumbers
                    dx = diff(xlim)/100;
                    th = text(dataX+dx, dataY, cellstr(num2str(obj.dispChannels')), 'FontSize', 8);
                    set(th, 'Clipping', 'on');
                    obj.numbers = [obj.numbers th];
                end
            end
            
            hold(obj.ax, 'off');
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
                    obj.showNumbers = ~obj.showNumbers;
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

