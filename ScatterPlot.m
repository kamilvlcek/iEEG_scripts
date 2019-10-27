classdef ScatterPlot < handle
    %SCATTERPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ieegdata
        header; %handle na kopii ieegdata.CH, kvuli PlotBrain
        dispChannels %zobrazene kanaly v grafu
        dispData; %vyobrazena data ve scatterplotu
        
        stats; % spocitane statistiky
        
        is3D
        
        selCh ; %kopie ieegdata.plotRCh.selCh;
        selChNames ; %kopie ieegdata.plotRCh.selChNames;
        dispSelCh % vyber klavesami fghjkl 
        dispSelChName
        
        dispFilterCh % CH.sortorder
        
        connectChannels
        connectionsPlot
        
        fig
        ax
        plots
        sbox
        pbox
        texthandles ; %handle na texty v obrazku, abych je mohl smazat
        
        highlights
        
        showNumbers; %zobrazit popisky hodnot, 0=nic 1=cisla 2=jmena kanalu
        numbers
        markerSize
        
        axisX
        axisY
        axisZ
        
        valFraction
        intFraction
        
        categories %seznam cisel kategorii od nejdulezitejsi (podle poradi ve statistice)
        categoryNames %jmena kategorii odpovidajici obj.categories
        categoriesSelectionIndex %1:numel(obj.categories)
        
        filterListener
        channelListener
        
        baseColors = [0 1 0; 0 0 1; 1 0 0; 1 1 0; 1 0 1; 0 1 0];
        categoryMarkers = {'o', 's', 'd', 'x'};
        transparent; % jestli maji byt markers nakreslene s pruhlednosti
    end
    
    methods
        function obj = ScatterPlot(ieegdata, is3D)
            %SCATTERPLOT Construct an instance of this class
            %   Detailed explanation goes here
            obj.ieegdata = ieegdata;
            obj.header = copy(obj.ieegdata.CH); %vytvorim nezavislou kopii kvuli PlotBrain, aby ne kolize z AUCPlotM grafem
            obj.selCh = ieegdata.plotRCh.selCh; %vyber kanalu fghjkl * pocet kanalu
            obj.selChNames = ieegdata.plotRCh.selChNames; %vyber kanalu fghjkl * pocet kanalu
            obj.dispSelChName = [];
            obj.dispSelCh = 1:size(obj.selCh,1);  % Zobrazuji vse
            
            obj.dispFilterCh = obj.ieegdata.CH.sortorder; % Vyber podle FilterChannels
            obj.connectChannels = 0;
            obj.showNumbers = 0; %defaultne se zadne labels nezobrazuji
            obj.numbers = [];
            obj.markerSize = 34;
            obj.transparent = false;
            
            if ~exist('is3D','var')
                obj.is3D = false; 
            else
                obj.is3D = is3D;
            end
            
            obj.setCategories();
            
            obj.drawScatterPlot(0.5, 0.5, 'tmax', 'valmax', 'category');
        end
        
        function drawScatterPlot(obj, valFraction, intFraction, axisX, axisY, axisZ)
            obj.setTriggerValues(valFraction, intFraction);
            obj.initAxes(axisX, axisY, axisZ);
            obj.updatePlot();
            obj.fixAxesLimits();
        end
        
        function PlotBrain(obj,katnum,xy,rangeZ)
            if ~exist('katnum','var') || isempty(katnum), katnum = 1; end
            if ~exist('xy','var'), xy = 'y'; end %defaultne osaY = valmax napriklad
            if ~exist('rangeZ','var') %pokud neni zadane
                rangeZ = iff(xy=='x',xlim(obj.ax),ylim(obj.ax));  %nastavi se podle limitu scatterplotu
            end 
            selChFiltered = obj.selCh(obj.dispChannels,:); %chci zobrazovat jen signif odpovedi
            iData = logical(selChFiltered(:,katnum)); %kanaly se signifikantim rozdilem vuci baseline v teto kategorii
            if(xy=='x')
                data = obj.dispData(katnum).dataX(iData);
                dataName = obj.axisX;
            else
                data = obj.dispData(katnum).dataY(iData);
                dataName = obj.axisY;
            end
            obj.header.plotCh3D.selch = []; %nechci mit vybrany zadny kanal z minula
            chshowstr = obj.header.plotCh2D.chshowstr; %musim udelat kopii jinak se do grafu preda odkaz
            obj.header.ChannelPlot([],0,data,... %param chnvals
                obj.dispChannels(iData),... %chnsel jsou cisla kanalu, pokud chci jen jejich vyber
                [],[],{[dataName '(' obj.categoryNames{katnum} '), SelCh: ' cell2str(obj.dispSelChName) ], ... %popis grafu = title - prvni radek
                ['show:' chshowstr]}, ... %popis grafu title, druhy radek
                rangeZ); %rozsah hodnot - meritko barevne skaly
            %set(obj.plotAUC.Eh.CH.plotCh3D.fh, 'WindowButtonDownFcn', {@obj.hybejPlot3Dclick, selch});
        end
        
        function setXYLim(obj,xrange,yrange) %set axes limit
            if ~exist('xrange','var') || isempty(xrange), xrange = xlim(obj.ax); end
            if ~exist('yrange','var') || isempty(yrange), yrange = ylim(obj.ax); end
            xlim(obj.ax,xrange);
            ylim(obj.ax,yrange);
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
 
        function initAxes(obj, axisX, axisY, axisZ)
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

            switch axisZ
                case 'valmax'
                    labelZ = 'v_{max}';
                case 'tmax'
                    labelZ = 't_{max}';
                case 'tfrac'
                    labelZ = ['t_{' num2str(obj.valFraction) '}'];
                case 'tint'
                    labelZ =  ['t_{int, ' num2str(obj.intFraction) '}'];
                case 'category'
                    labelZ = 'category';
                case 'channel'
                    labelZ = 'channel';
                otherwise
                    disp('Z axis specification must be one of: valmax, tmax, tfrac, tint, channel, category');
                    return;
            end
            
            obj.axisX = axisX;
            obj.axisY = axisY;
            obj.axisZ = axisZ;

            obj.fig = figure('CloseRequestFcn', @obj.tearDownFigCallback, 'Name', 'ScatterPlot');
            
            obj.ax = axes(obj.fig);
            xlabel(obj.ax, labelX);
            ylabel(obj.ax, labelY);
            zlabel(obj.ax, labelZ);

            obj.filterListener = addlistener(obj.ieegdata.CH, 'FilterChanged', @obj.filterChangedCallback);
            obj.channelListener = addlistener(obj.ieegdata, 'SelectedChannel', 'PostSet', @obj.channelChangedCallback);
            
            set(obj.fig, 'KeyPressFcn', @obj.hybejScatterPlot);
            set(obj.fig, 'WindowButtonDownFcn', @obj.hybejScatterPlotClick);
        end
        
        function setTriggerValues(obj, valFraction, intFraction)
            obj.valFraction = valFraction;
            obj.intFraction = intFraction;
        end
        
        function fixAxesLimits(obj)
            xlim(xlim(obj.ax));
            ylim(ylim(obj.ax));
        end
        
        function updatePlot(obj,recompute)
            if ~exist('recompute','var'), recompute = 1; end
            obj.setDisplayedChannels(); % Kombinace voleb pro zobrazeni kanalu
            selChFiltered = obj.selCh(obj.dispChannels,:); %filter kanalu ve vyberu fghjkl
            RjCh = intersect(intersect(obj.dispFilterCh, obj.dispSelCh),obj.ieegdata.RjCh); %vyrazene kanaly z tech nyni zobrazenych - cisla kanalu
            selChRj = obj.selCh(RjCh,:); %filter vyrazenych kanalu ve vyberu fghjkl
            delete(obj.plots); obj.plots = [];
            delete(obj.sbox);
            delete(obj.pbox);            
            delete(obj.connectionsPlot); obj.connectionsPlot = [];
            delete(obj.numbers); obj.numbers = [];  
            for j = 1:numel(obj.texthandles)
                if isgraphics(obj.texthandles(j)) && obj.texthandles(j) ~= 0
                    delete(obj.texthandles(j));  
                end
            end            
            obj.texthandles = [];
            
            if ~isempty(obj.dispSelChName)
                obj.sbox = annotation(obj.fig, 'textbox',[0 .9 .4 .1], 'String', obj.dispSelChName, 'EdgeColor', 'none');
            end
            
            catlist = strjoin(obj.categoryNames(obj.categoriesSelectionIndex), ', ');
            obj.pbox = annotation(obj.fig, 'textbox', [0 0 .4 .1], 'String', ['C: ' catlist], 'EdgeColor', 'none');
            
           
            if isempty(obj.dispChannels)
                disp('No channels corresponding to the selection');
                return;
            end
            if recompute
                obj.stats = struct();
                for k = obj.categoriesSelectionIndex %1:numel(obj.categories)
                    catnum = obj.categories(k);
                    [obj.stats(k).valmax, obj.stats(k).tmax, obj.stats(k).tfrac, obj.stats(k).tint] = obj.ieegdata.ResponseTriggerTime(obj.valFraction, obj.intFraction, catnum, obj.dispChannels);
                end
            end
            
            hold(obj.ax, 'on');
            legend(obj.ax, 'off');
            
            if obj.connectChannels > 0
                obj.drawConnectChannels();
            end
            obj.drawPlot(selChFiltered, selChRj);
            
            legend(obj.ax, 'show');
            hold(obj.ax, 'off');
            if isfield(obj.ieegdata.CH.plotCh2D, 'chshowstr') && ~isempty(obj.ieegdata.CH.plotCh2D.chshowstr)
                title(['show:' obj.ieegdata.CH.plotCh2D.chshowstr]);
            else
                title('show: all');
            end
            
            if obj.is3D
                grid(obj.ax, 'on');
            else
                grid(obj.ax, 'off');
            end
        end
        
        function drawPlot(obj, selChFiltered, selChRj)
            pocty = zeros(numel(obj.categoriesSelectionIndex),3); %pocty zobrazenych kanalu  
            for k = flip(obj.categoriesSelectionIndex) %nejpozdeji, cili nejvic na vrchu, chci mit prvni kategorii %1:numel(obj.categories)
                dataX = obj.stats(k).(obj.axisX);
                dataY = obj.stats(k).(obj.axisY);
                if obj.is3D %TODO: Zobecnit a ukladat nekde na zacatku (u vypoctu stats?) - opakuje se to napr. v drawConnectChannels, highlighSelected, nebo hybejPlotClick
                    if strcmp(obj.axisZ, 'channel')
                        dataZ = obj.dispChannels;
                    elseif strcmp(obj.axisZ, 'category')
                        dataZ = k * ones(size(dataX));
                    else
                        dataZ = obj.stats(k).(obj.axisZ);
                    end
                end
                iData = logical(selChFiltered(:,k)); %kanaly se signifikantim rozdilem vuci baseline v teto kategorii
                if any(iData)
                    if obj.is3D
                        obj.plots(k,1) = scatter3(obj.ax, dataX(iData), dataY(iData), dataZ(iData), obj.markerSize, repmat(obj.baseColors(k,:), sum(iData), 1), obj.categoryMarkers{k}, 'MarkerFaceColor', 'flat', 'DisplayName', obj.categoryNames{k});
                        if obj.transparent, alpha(obj.plots(k,1),.5); end %volitelne pridani pruhlednosti
                    else
                        obj.plots(k,1) = scatter(obj.ax, dataX(iData), dataY(iData), obj.markerSize, repmat(obj.baseColors(k,:), sum(iData), 1), obj.categoryMarkers{k}, 'MarkerFaceColor', 'flat', 'DisplayName', obj.categoryNames{k});
                        if obj.transparent, alpha(obj.plots(k,1),.5); end %volitelne pridani pruhlednosti
                    end
                    pocty(k,1) = sum(iData); %pocet signifikantnich v teto kategorii
                end
                iData = ~iData; %kanaly bez signif rozdilu vuci baseline v teto kategorii
                if any(iData) && obj.connectChannels >= 0
                    if obj.is3D
                        obj.plots(k,2) = scatter3(obj.ax, dataX(iData), dataY(iData), dataZ(iData), obj.markerSize, repmat(obj.baseColors(k,:), sum(iData), 1), obj.categoryMarkers{k}, 'MarkerFaceColor', 'none', 'DisplayName', obj.categoryNames{k},...
                            'HandleVisibility','off'); %nebude v legende
                    else
                        obj.plots(k,2) = scatter(obj.ax, dataX(iData), dataY(iData), obj.markerSize, repmat(obj.baseColors(k,:), sum(iData), 1), obj.categoryMarkers{k}, 'MarkerFaceColor', 'none', 'DisplayName', obj.categoryNames{k},...
                            'HandleVisibility','off'); %nebude v legende
                    end
                    pocty(k,2) = sum(iData); %pocet ne signifikantnich v teto kategorii
                end
                if any(logical(selChRj(:,k))) %pocet rejectovanych pro tuto kategorii
                    pocty(k,3) = sum(logical(selChRj(:,k))); 
                end
                if obj.showNumbers > 0
                    if obj.showNumbers == 1
                        labels = cellstr(num2str(obj.dispChannels')); %cisla zobrazenych kanalu 
                    elseif obj.showNumbers==2
                        labels = {obj.ieegdata.CH.H.channels(obj.dispChannels).name}'; %jmena kanalu
                    else
                        labels = {obj.ieegdata.CH.H.channels(obj.dispChannels).neurologyLabel}'; %anatomicke oznaceni kanalu
                    end
                    dx = diff(xlim)/100;
                    iData = iff( obj.connectChannels >= 0 , true(size(selChFiltered,1),1), logical(selChFiltered(:,k)) );   
                    if obj.is3D
                        th = text(dataX(iData)+dx, dataY(iData), obj.dispChannels(iData), labels(iData), 'FontSize', 8);
                    else
                        th = text(dataX(iData)+dx, dataY(iData), labels(iData), 'FontSize', 8);
                    end
                    set(th, 'Clipping', 'on');
                    obj.numbers = [obj.numbers ; th]; % do dlouheho sloupce, kazda kategorie ma ruzny pocet zobrazenych kanalu a textu
                end
                obj.dispData(k).dataX = dataX; %zalohuju vyobrazena data pro jine pouziti
                obj.dispData(k).dataY = dataY;
                
                xsize = xlim(obj.ax); ysize = ylim(obj.ax);
                xtext = xsize(2)-diff(xsize)/5;
                ytext = ysize(2)-diff(ysize)/5;
                
                obj.texthandles(k) = text(xtext, ytext - (k-1)*diff(ysize)/20 ,[obj.categoryNames{k} ':' num2str(pocty(k,1)) '+' num2str(pocty(k,2)) ' (' num2str(pocty(k,3)) ')' ]);
            end
        end
        
        function drawConnectChannels(obj)
        % Nakresli linku spojujici stejne kanaly. Ruzne barvy musi byt samostatny plot (aby mohl scatter zustat ve stejnych osach)
            if length(obj.categoriesSelectionIndex) > 1
                validCategoryIndex = obj.categoriesSelectionIndex(1);
                catIndex = zeros(size(obj.categoriesSelectionIndex));
                x = zeros(length(obj.categoriesSelectionIndex), length(obj.stats(validCategoryIndex).(obj.axisX)));
                y = zeros(length(obj.categoriesSelectionIndex), length(obj.stats(validCategoryIndex).(obj.axisX)));
                if ~strcmp(obj.axisZ, 'channel') && ~strcmp(obj.axisZ, 'category')
                    z = zeros(length(obj.categoriesSelectionIndex), length(obj.stats(validCategoryIndex).(obj.axisX)));
                end
                for cat = 1:length(obj.categoriesSelectionIndex)
                    catIndex(cat) = obj.categoriesSelectionIndex(cat);
                    x(cat,:) = obj.stats(catIndex(cat)).(obj.axisX);
                    y(cat,:) = obj.stats(catIndex(cat)).(obj.axisY);
                    if ~strcmp(obj.axisZ, 'channel') && ~strcmp(obj.axisZ, 'category')
                        z(cat,:) = obj.stats(catIndex(cat)).(obj.axisZ);
                    end
                end
                l = length(obj.stats(catIndex(1)).(obj.axisX));
                for k = 1:l
                    xx = [x(:,k); x(1,k)]; %tri body za vsechny kategorie + pridam prvni na konec znovu, aby se spojily
                    yy = [y(:,k); y(1,k)]; %kvuli zrychleni u velkeho mnozstvi bodu
                    if obj.is3D
                        if strcmp(obj.axisZ, 'channel')
                            zz = obj.dispChannels(k) * ones(size(xx));
                        elseif strcmp(obj.axisZ, 'category')
                            zz = [obj.categoriesSelectionIndex obj.categoriesSelectionIndex(1)] * ones(1, size(xx,2));
                        else
                            zz = [z(:,k); z(1,k)];
                        end
                                
                        obj.connectionsPlot(end+1) = plot3(xx, yy, zz, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
                    else
                        obj.connectionsPlot(end+1) = plot(xx, yy, 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
                    end
                end
            else
                disp('No categories to connect');
                obj.connectChannels = -1;
            end
        end
                
        function highlightSelected(obj, ch)
            if ~isempty(obj.highlights)
                delete(obj.highlights);
            end
            obj.highlights = [];
            idx = find(obj.dispChannels == ch);
            hold(obj.ax, 'on');
            if idx
                for k = obj.categoriesSelectionIndex
                    dataX = obj.stats(k).(obj.axisX);
                    dataY = obj.stats(k).(obj.axisY);
                    if obj.is3D
                        if strcmp(obj.axisZ, 'channel')
                            dataZ = obj.dispChannels;
                        elseif strcmp(obj.axisZ, 'category')
                            dataZ = k * ones(size(dataX));
                        else
                            dataZ = obj.stats(k).(obj.axisZ);
                        end
                        obj.highlights(k) = scatter3(obj.ax, dataX(idx), dataY(idx), dataZ(idx), 3*obj.markerSize, 0.75*obj.baseColors(k,:), 'o', 'MarkerFaceColor', 'none', 'LineWidth', 2, 'HandleVisibility','off');
                    else
                        obj.highlights(k) = scatter(obj.ax, dataX(idx), dataY(idx), 3*obj.markerSize, 0.75*obj.baseColors(k,:), 'o', 'MarkerFaceColor', 'none', 'LineWidth', 2, 'HandleVisibility','off');
                    end
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
                    obj.connectChannels = obj.connectChannels + 1;
                    if obj.connectChannels >1,  obj.connectChannels = -1; end
                    obj.updatePlot(0);
                case {'t'}
                    obj.transparent =  ~obj.transparent;
                    if obj.connectChannels >1,  obj.connectChannels = -1; end
                    obj.updatePlot(0); %neprepocitat hodnoty
                case {'n'}
                    obj.showNumbers =  obj.showNumbers + 1;
                    if obj.showNumbers > 3, obj.showNumbers=0; end %0->1->2->3->0
                    obj.updatePlot(0); %neprepocitat hodnoty
                case {'add'}
                    obj.markerSize = obj.markerSize + 8;
                    obj.updatePlot();
                case {'subtract'}
                    obj.markerSize = max(2, obj.markerSize - 8);
                    obj.updatePlot();
            end
        end
        
        function hybejScatterPlotClick(obj,h,~)
            if isempty(obj.dispChannels) % Pokud nejsou zobrazene zadne kanaly, nedelam nic
              return;
            end

            mousept = get(gca,'currentPoint');
            p1 = mousept(1,:); p2 = mousept(2,:); % souradnice kliknuti v grafu - predni a zadni bod
            chs  = zeros(size(obj.categoriesSelectionIndex));
            dist = zeros(size(obj.categoriesSelectionIndex));
            for k = 1:length(obj.categoriesSelectionIndex) % vsechny zobrazene kategorie
              categoryIndex = obj.categoriesSelectionIndex(k);
              dataX = obj.stats(categoryIndex).(obj.axisX);
              dataY = obj.stats(categoryIndex).(obj.axisY);
              if obj.is3D
                  if strcmp(obj.axisZ, 'channel')
                      dataZ = obj.dispChannels;
                  elseif strcmp(obj.axisZ, 'category')
                      dataZ = k * ones(size(dataX));
                  else
                      dataZ = obj.stats(k).(obj.axisZ);
                  end
                  coordinates = [dataX; dataY; dataZ]; % souradnice zobrazenych kanalu
                  [chs(k), dist(k)] = findClosestPoint(p1, p2, coordinates, 0.05);    % najdu kanal nejblize mistu kliknuti
              else
                  x = p1(1); y = p1(2); % souradnice v grafu (ve 2D pouze "predni" bod)
                  [chs(k), dist(k)] = dsearchn([dataX' dataY'], [x y]); %najde nejblizsi kanal a vzdalenost k nemu
                  if dist(k) > mean([diff(ylim(obj.ax)), diff(xlim(obj.ax))])/20 % kdyz kliknu moc daleko od kanalu, nechci nic vybrat - nastavim [0 inf] stejne jako to dela funkce findClosestPoint
                      chs(k) = 0;
                      dist(k) = inf;
                  end
              end
            end

            [mindist, k_min] = min(dist); % vyberu skutecne nejblizsi kanal ze vsech kategorii

            if mindist < inf
              ch = obj.dispChannels(chs(k_min));
              %TODO: Pokud neni otevreny PlotResponseCh, nebude po otevreni znat cislo vybraneho kanalu. Lepsi by bylo pouzit proxy objekt, ktery drzi informaci o vybranem kanalu a v pripade zmeny vyberu posle signal, ktery se tak zpropaguje do vsech plotu, ktere ho potrebuji.
              if isfield(obj.ieegdata.plotRCh, 'fh') && isvalid(obj.ieegdata.plotRCh.fh)  % Zjistim, jeslti je otevreny PlotResponseCh
                  sortChannel = find(obj.ieegdata.CH.sortorder == ch);
                  obj.ieegdata.PlotResponseCh(sortChannel);    % Pokud mam PlotResponseCh, updatuju zobrezene kanaly
                  % Nevolam highlightSelected, protoze ten se zavola diky eventu
                  figure(obj.fig); %kamil - dam do popredi scatter plot
              else
                  obj.highlightSelected(ch);
                  %TODO: Pokud ted manualne otevru PlotResponseCh bez parametru, neuvidim v nem spranvy kanal
              end
            else
                obj.highlightSelected(0);   % zrusi vyber
            end
        end

        
        function setDisplayedChannels(obj)
            obj.dispChannels = intersect(obj.dispFilterCh, obj.dispSelCh); % CH.sortorder (vysledek CH.FilterChannels) & vyber klavesami fghjkl 
            obj.dispChannels = setdiff(obj.dispChannels,obj.ieegdata.RjCh); %kamil 15.10 - vyradim ze zobrazeni vyrazene kanaly
        end
        
        function filterChangedCallback(obj,~,~)
            obj.dispFilterCh = obj.ieegdata.CH.sortorder;   % Zmena vyberu dle filtru
            obj.updatePlot();
        end

        function channelChangedCallback(obj, ~, eventData)
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %TODO: Pri kliknuti do ScatterPlotu se tohle zavola dvakrat!!! Nejspis problem se skoro-zacyklenim z PlotResponseCh. Sice to navenek funguje spravne, ale dvoji volani je nesmysl.
            obj.highlightSelected(eventData.AffectedObject.SelectedChannel);
            %disp(['change in SP: ' num2str(eventData.AffectedObject.SelectedChannel)]);
        end
        
        function tearDownFigCallback(obj,src,~)
            delete(obj.filterListener);
            delete(obj.channelListener);
            delete(src);
        end      
        
    end
    
end

