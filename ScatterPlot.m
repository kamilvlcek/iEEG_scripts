classdef ScatterPlot < handle
    %SCATTERPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ieegdata
        dispChannels
        
        selCh
        selChNames
        dispSelChName
        
        connectPairs
        pairsPlot
        
        fig
        ax
        plots
        sbox
        
        axisX
        axisY
        
        valFraction
        intFraction
        
        categories
        categoryNames
        categoriesSelectionIndex
    end
    
    methods
        function obj = ScatterPlot(ieegdata)
            %SCATTERPLOT Construct an instance of this class
            %   Detailed explanation goes here
            obj.ieegdata = ieegdata;
            obj.selCh = ieegdata.plotRCh.selCh; %vyber kanalu fghjkl * pocet kanalu
            obj.selChNames = ieegdata.plotRCh.selChNames; %vyber kanalu fghjkl * pocet kanalu
            obj.dispSelChName = [];
            %TODO: Hodilo by se zachytit zmenu FilterChannles pres event?
            obj.dispChannels = obj.ieegdata.CH.sortorder; % Vyber podle FilterChannels
            obj.connectPairs = false;
            
            obj.setCategories();
            
            obj.drawScatterPlot(0.5, 0.5, 'tmax', 'valmax');
        end
        
        function drawScatterPlot(obj, valFraction, intFraction, axisX, axisY)
            obj.setTriggerValues(valFraction, intFraction);
            obj.initAxes(axisX, axisY);
            obj.updateScatterPlot();
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
            
            obj.fig = figure;
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
        
        function updateScatterPlot(obj)
            delete(obj.plots);
            obj.plots = [];
            
            delete(obj.sbox);
            if ~isempty(obj.dispSelChName)
                obj.sbox = annotation(obj.fig, 'textbox',[0 .9 .4 .1],'String',obj.dispSelChName,'EdgeColor','none');
            end
            
            stats = struct();
            for k = obj.categoriesSelectionIndex
                catnum = obj.categories(k);
                [stats(k).valmax, stats(k).tmax, stats(k).tfrac, stats(k).tint] = obj.ieegdata.ResponseTriggerTime(obj.valFraction, obj.intFraction, catnum, obj.dispChannels);
            end
            
            hold(obj.ax, 'on');
            for k = obj.categoriesSelectionIndex
                dataX = stats(k).(obj.axisX);
                dataY = stats(k).(obj.axisY);
                        
                obj.plots(k) = scatter(obj.ax, dataX, dataY, 8, 'filled', 'DisplayName', obj.categoryNames(k));
                legend(obj.ax);
            end
            
            delete(obj.pairsPlot); obj.pairsPlot = [];
            if obj.connectPairs     % Nakresli linku spojujici prislusny par. Pro rychlejsi vykreslovani je pouzit jeden plot, zdrojova data jsou dvojice souradnic za sebou, oddelene NaNem
                if length(obj.categoriesSelectionIndex) == 2
                    k1 = obj.categoriesSelectionIndex(1);
                    k2 = obj.categoriesSelectionIndex(2);
                    l = length(stats(k1).(obj.axisX));
                    conx = NaN(1, 3*l);
                    conx(1:3:3*l) = stats(k1).(obj.axisX);
                    conx(2:3:3*l) = stats(k2).(obj.axisX);
                    cony = NaN(1, 3*l);
                    cony(1:3:3*l) = stats(k1).(obj.axisY);
                    cony(2:3:3*l) = stats(k2).(obj.axisY);
                    obj.pairsPlot = plot(obj.ax, conx, cony, 'Color', [0.8, 0.8, 0.8]);
                else
                    disp('Pair connection must include exactly 2 categories');
                    obj.connectPairs = false;
                end
            end
            hold(obj.ax, 'off');
        end
        
        function hybejScatterPlot(obj,~,eventDat)
            switch eventDat.Key
                case {'u','i','o','p'}
                    disp(obj.categoriesSelectionIndex);
                    ik = find('uiop'==eventDat.Key); % index 1-4
                    if ik <= numel(obj.categories)
                        if(ismember(ik, obj.categoriesSelectionIndex))
                            obj.categoriesSelectionIndex = setdiff(obj.categoriesSelectionIndex, ik);
                        else
                            obj.categoriesSelectionIndex = union(obj.categoriesSelectionIndex, ik);
                        end
                    end    
                    obj.updateScatterPlot();
                case {'a'}
                    obj.dispChannels = obj.ieegdata.CH.sortorder; % Vyber podle FilterChannels
                    obj.dispSelChName = [];
                    obj.updateScatterPlot();
                case {'f','g','h','j','k','l'}
                    obj.dispChannels = find(obj.selCh(:,'fghjkl'==eventDat.Key))';
                    obj.dispSelChName = obj.selChNames{'fghjkl'==eventDat.Key};
                    obj.updateScatterPlot();
                case {'s'}
                    obj.connectPairs = ~obj.connectPairs;
                    obj.updateScatterPlot();
            end            
        end
        
    end
    
end

