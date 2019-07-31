classdef ScatterPlot < handle
    %SCATTERPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ieegdata
        gui
        vizTypes
        figs
        
        categories
        categoryNames
        neurologyLabels
        
        categoriesSelection
        neurologyLabelsSelection
    end
    
    methods
        function obj = ScatterPlot(ieegdata)
            %SCATTERPLOT Construct an instance of this class
            %   Detailed explanation goes here
            obj.ieegdata = ieegdata;
            obj.setVizTypes();
            obj.setCategories();
            obj.setNeurologyLabels();
            obj.CreateInterface();
        end
    end
    
    methods(Access = private)
        function setVizTypes(obj)
            vizList = {
                'Scatter plot'          'drawScatterPlot'
            };
            settings = {
                struct('val_fraction', 0.5, 'int_fraction', 0.5, ...
                    'x_valmax', 0, 'x_tmax', 1, 'x_tfrac', 0, 'x_tint', 0, ...
                    'y_valmax', 1, 'y_tmax', 0, 'y_tfrac', 0, 'y_tint', 0) ...
            };
            obj.vizTypes = struct( ...
                'names', {vizList(:,1)'}, ...
                'functions', {vizList(:,2)'}, ...
                'settings', {settings}, ...
                'selected', 1 ...
            );
        end
        
        function updateVizSettings(obj)
            if ishandle(obj.gui.vizSettings)
                delete(obj.gui.vizSettings);
            end
            obj.gui.vizSettings = uix.VBox( 'Parent', obj.gui.vizSettingsLayout);
            vset = obj.gui.vset;
            st = obj.vizTypes.settings{obj.vizTypes.selected};
            hs = [];
            obj.addNumval('val_fraction', 'Value onset fraction: ', vset, st); hs = [hs 20];
            obj.addNumval('int_fraction', 'Integral fraction: ',    vset, st); hs = [hs 20];
            
            layoutAxes = uix.HBox('Parent', obj.gui.vizSettings);
            
            layoutX =  uix.VBox('Parent', layoutAxes);

            uicontrol('Parent', layoutX, 'Style', 'Text', 'String', 'x axis');
            
            obj.addCheckbox('x_valmax', 'maximum value', vset, st, layoutX);
            obj.addCheckbox('x_tmax',   'maximum time', vset, st, layoutX);
            obj.addCheckbox('x_tfrac',  'fraction time', vset, st, layoutX);
            obj.addCheckbox('x_tint',   'integral fraction time', vset, st, layoutX);
            
            layoutY =  uix.VBox('Parent', layoutAxes);
            
            uicontrol('Parent', layoutY, 'Style', 'Text', 'String', 'y axis');
            obj.addCheckbox('y_valmax', 'maximum value', vset, st, layoutY);
            obj.addCheckbox('y_tmax',   'maximum time', vset, st, layoutY);
            obj.addCheckbox('y_tfrac',  'fraction time', vset, st, layoutY);
            obj.addCheckbox('y_tint',   'integral fraction time', vset, st, layoutY);
            
            hs = [hs 80];
            
            set(obj.gui.vizSettings, 'Heights', hs);
        end
        
        function CreateInterface(obj)
            obj.gui = struct();
            obj.gui.vset = struct();
            obj.gui.window = figure( 'Name', 'Scatter Plot Controller', ...
                'NumberTitle', 'off', ...
                'MenuBar', 'none', ...
                'Toolbar', 'none', ...
                'HandleVisibility', 'off', ...
                'Position', [800 300 1200 1000]);

            mainLayout = uix.VBoxFlex( 'Parent', obj.gui.window, 'Spacing', 3);

            wrapperLayout = uix.HBox( 'Parent', mainLayout );

            obj.gui.vizSettingsLayout = uix.VBox('Parent', wrapperLayout );
            obj.gui.vizSettingsPanel = uiextras.BoxPanel( ...
                'Parent', obj.gui.vizSettingsLayout, ...
                'Title', 'Visualization settings:' );
            
            obj.gui.vizSettings = uix.VBox( 'Parent', obj.gui.vizSettingsLayout);
            obj.updateVizSettings();
            
            set(obj.gui.vizSettingsLayout, 'Heights', [28 -1]);

            obj.gui.categoriesLayout = uix.VBox('Parent', wrapperLayout );
            obj.gui.categoriesPanel = uiextras.BoxPanel( ...
                'Parent', obj.gui.categoriesLayout, ...
                'Title', 'Categories:' );
            obj.gui.categoriesList = uicontrol( 'Style', 'list', ...
                'Parent', obj.gui.categoriesLayout, ...
                'Callback', @obj.onCategoriesSelection );
            set(obj.gui.categoriesLayout, 'Heights', [28 -1]);

            obj.gui.neurologyLabelsLayout = uix.VBox('Parent', wrapperLayout );
            obj.gui.neurologyLabelsPanel = uiextras.BoxPanel( ...
                'Parent', obj.gui.neurologyLabelsLayout, ...
                'Title', 'Neurology Labels:' );
            obj.gui.neurologyLabelsList = uicontrol( 'Style', 'list', ...
                'Parent', obj.gui.neurologyLabelsLayout, ...
                'Callback', @obj.onNeurologyLabelsSelection );
            set(obj.gui.neurologyLabelsLayout, 'Heights', [28 -1]);
            
            controls_layout = uix.HBox( 'Parent', mainLayout );
            obj.gui.draw_button = uicontrol( 'Style', 'PushButton',...
                'Parent', controls_layout, ...
                'String', 'Draw!', ...
                'Callback', @obj.onDrawButton);

            set(mainLayout, 'Heights', [-1 28]);

            obj.resetCategories();
            obj.resetNeurologyLabels();
            
            obj.updateInterface();
        end
      
        function updateInterface(obj)
            obj.updateVizSettings();
        end
        
        function addCheckbox(obj, varname, title, vset, st, parent)
            if nargin == 4
                parent = obj.gui.vizSettings;
            end
            vset.(varname) = uicontrol('Parent', parent, 'Style', 'Checkbox', 'String', title, ...
                        'Value', st.(varname), 'Callback', {@obj.onVsetCheckbox, obj.vizTypes.selected, varname});
        end

        function addNumval(obj, varname, title, vset, st)
            layout_loc = uix.HBox('Parent', obj.gui.vizSettings);
            uicontrol('Parent', layout_loc, 'Style', 'Text', 'String', title);
            vset.(varname) = uicontrol('Parent', layout_loc, 'Style', 'Edit', 'String', st.(varname), ...
                    'Value', st.(varname), 'Callback', {@obj.onVsetNumval, obj.vizTypes.selected, varname});
        end

        function onVsetCheckbox(obj, src, ~, vid, varname)
            obj.vizTypes.settings{vid}.(varname) = get(src, 'Value');
            obj.updateInterface();
        end

        function onVsetNumval(obj, src, ~, vid, varname)
            obj.vizTypes.settings{vid}.(varname) = str2double(get(src, 'String'));
            obj.updateInterface();
        end
        
        
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
        end
        
        function setNeurologyLabels(obj)
            obj.neurologyLabels = unique({obj.ieegdata.CH.H.channels.neurologyLabel}, 'sorted');
        end
        
        function onCategoriesSelection(obj, src, ~)
            obj.categoriesSelection = get(src, 'Value');
            obj.updateInterface();
        end
        
        function onNeurologyLabelsSelection(obj, src, ~)
            obj.neurologyLabelsSelection = get(src, 'Value');
            obj.updateInterface();
        end
        
        function onDrawButton(obj, ~, ~)
            func = obj.vizTypes.functions(obj.vizTypes.selected);
            obj.updateInterface();
            obj.(func{1});
        end
        
        function resetCategories(obj)
            if isempty(obj.categories)
                set(obj.gui.categoriesList, 'String', 'No categories loaded.');
            else
                len = max(2, length(obj.categories));   % Listbox musi mit alespon 2 hodnoty
                set(obj.gui.categoriesList, 'Max', len);
                set(obj.gui.categoriesList, 'String', obj.categoryNames);
            end
        end
        
        function resetNeurologyLabels(obj)
            if isempty(obj.neurologyLabels)
                set(obj.gui.neurologyLabelsList, 'String', 'No neurology labels loaded.');
            else
                len = max(2, length(obj.neurologyLabels));   % Listbox musi mit alespon 2 hodnoty
                set(obj.gui.neurologyLabelsList, 'Max', len);
                set(obj.gui.neurologyLabelsList, 'String', obj.neurologyLabels);
            end
        end
        
        function drawScatterPlot(obj)
            graphSettings = obj.vizTypes.settings{obj.vizTypes.selected};
            
            labels = {'v_{max}' 't_{max}' ['t_{' num2str(graphSettings.val_fraction) '}'] ['t_{int, ' num2str(graphSettings.int_fraction) '}']};
            
            axesX = {};
            labelsX = {};
            if graphSettings.x_valmax
                axesX(end+1) = {'valmax'};
                labelsX(end+1) = labels(1);
            end
            if graphSettings.x_tmax
                axesX(end+1) = {'tmax'};
                labelsX(end+1) = labels(2);
            end
            if graphSettings.x_tfrac
                axesX(end+1) = {'tfrac'};
                labelsX(end+1) = labels(3);
            end
            if graphSettings.x_tint
                axesX(end+1) = {'tint'};
                labelsX(end+1) = labels(4);
            end
            
            axesY = {};
            labelsY = {};
            if graphSettings.y_valmax
                axesY(end+1) = {'valmax'};
                labelsY(end+1) = labels(1);
            end
            if graphSettings.y_tmax
                axesY(end+1) = {'tmax'};
                labelsY(end+1) = labels(2);
            end
            if graphSettings.y_tfrac
                axesY(end+1) = {'tfrac'};
                labelsY(end+1) = labels(3);
            end
            if graphSettings.y_tint
                axesY(end+1) = {'tint'};
                labelsY(end+1) = labels(4);
            end

            obj.figs = struct('axisX', {}, 'axisY', {}, 'figure', {}, 'axes', {}, 'plots', {}, 'highlights', {});
            
            for k = obj.categoriesSelection
                catnum = obj.categories(k);
                stats = struct();
                [stats.valmax, stats.tmax, stats.tfrac, stats.tint] = obj.ieegdata.ResponseTriggerTime(graphSettings.val_fraction, graphSettings.int_fraction, catnum);
                for xx = 1:length(axesX)
                    axisX = stats.(axesX{xx});
                    for yy = 1:length(axesY)
                        axisY = stats.(axesY{yy});
                        
                        figFilter = strcmp({obj.figs.axisX}, axesX{xx}) & strcmp({obj.figs.axisY}, axesY{yy});
                        fig = obj.figs(figFilter);
                        if isempty(fig)
                            fig = struct('axisX', axesX{xx}, 'axisY', axesY{yy}, 'figure', figure, 'axes', '', 'plots', [], 'highlights', []);
                            fig.axes = axes(fig.figure);
                            xlabel(fig.axes, labelsX{xx});
                            ylabel(fig.axes, labelsY{yy});
                            obj.figs(end+1) = fig;
                        elseif length(fig) > 1
                            disp('Error: mulitple graphs retrieved for single query');
                        end
                        
                        ax = fig.axes;
                        hold(ax, 'on');
                        fig.plots(end+1) = scatter(ax, axisX, axisY, '.', 'DisplayName', obj.categoryNames(k));
                        legend(ax);
                        hold off;
                        
                        obj.figs(figFilter) = fig;
                    end
                end
            end
        end
        
    end
    
end

