classdef CHilbertL < CHilbert
    % HILBERT.CLASS extension for CHilbert class
    %   Lukas Hejtmanek
    
    properties (Access = public)
        % HFreq; % hilberova obalka pro kazde frekvenci pasmo - time x channel x freq (x kategorie)
        % HFreqEpochs; %Hf bez priemerovania cez epochy - time x channel x frequency x epoch
        % fphaseEpochs; % epochovane fazy fphase ~HFreqEpochs
        % Hf; %frekvencni pasma pro ktere jsou pocitany obalky - okraje pasem, pocet je tedy vetsi o 1 nez pocet spocitanych pasem
        % Hfmean; %stredni hodnoty pasem  - pocet = pocet spocitanych pasem
        % hfilename; %jmeno souboru CHilbert  
        % plotF = struct; %udaje o stavu plotu PlotResponseFreq
        % fphase; %faze vsech zpracovavanych frekvenci - premiestnene z CMorlet pre vykreslenie a porovnanie faz z MW a Hilberta do buducna        
        % frealEpochs; % epochovane filtrovane eeg
        %normalization; %typ normalizace
    end
    
    % -------------- public instance methods -------------------------
    methods (Access = public)
        %% Constructor description??
        % Just calls CHilbert class constructor
        function obj = CHilbertL(varargin)
            obj@CHilbert(varargin{:});
        end
        
        %% Plot Response Frequency
        % ch: max number of channels to plot (from 1 to ch)
        % categories: which category to plot - if not defined, plots all.
        % Takes array of category indices beginning with 1? e.g. [1,3]
        % TODO - DEPRECATE the possibility to use cells in categories {}
        function obj = PlotResponseFreq(obj, ch, categories)
            if ~exist('ch', 'var'), ch = []; end
            if ~exist('categories', 'var'), categories = []; end
           
            obj.prepareplotcategoriespowertime(ch, categories);
            % Shouldn't the Ch be obj.plotF.ch? - WHAT IS THE
            % DIFFERENCE? the CHHeader has it as a field not afunction
            % and thus it is very difficult to unravel :(
            % The original function has the ch in plotF.ch, but uses the
            % obj.CH.sortorder(ch) to plot things - don't fully understand
            % if that is correct
            obj.plotcategoriespowertime(obj.CH.sortorder(ch), obj.plotF.kategories);
            obj.plotlabels(obj.CH.sortorder(ch), obj.plotF.kategories);
            set(obj.plotF.fh, 'KeyPressFcn', @obj.hybejPlotF);
        end
        
        % Buffering of the plot or taking settings from saved structs
        % if we are only redrawing existing plots
        function obj = prepareplotcategoriespowertime(obj, ch, categories)
            if ~numel(ch) == 0 && ~isfield(obj.plotF, 'ch'), obj.plotF.ch = 1;
            else, obj.plotF.ch = ch; end
            
            % Only rewrites categorties if not already set
            if numel(categories) > 0, obj.plotF.kategories = categories;
            elseif ~isfield(obj.plotF, 'kategories') && numel(obj.plotF.kategories) == 0
                if isfield(obj.Wp(obj.WpActive), 'kats'), categories = obj.Wp(obj.WpActive).kats;
                else, categories = obj.PsyData.Categories(); end
                obj.plotF.kategories = categories;
            end
            
            % redraws the plot or creates one
            if isfield(obj.plotF, 'fh') && ishandle(obj.plotF.fh)
                figure(obj.plotF.fh);
            else
                obj.plotF.fh = figure('Name', 'ResponseFreq', 'Position', [20, 500, 1200, 300]);
                colormap jet;
            end
        end
        
        % TODO - deprecate option to pass cells into categories
        function obj = plotcategoriespowertime(obj, ch, categories)
             [miny, maxy] = obj.getplotylimit();
             for k = 1:numel(categories)
                % QUESTION - If I pass it a cell I'd expect the
                % cell to contain NAMES of categories, but that is not what
                % is happening. The categories are not taken as strings,
                % and are not correlated to obj.PsyData.Category name.
                % Neither here, nor in PlotLabels. They are merely taken in
                % their order - this should be DEPRECATED
                if iscell(categories(k))                  
                    dd = zeros(size(obj.HFreq, 1), size(obj.HFreq, 3), numel(categories{k}));
                    for ikat = 1:numel(categories{k})
                        dd(:, :, ikat) = obj.getaverageenvelopes(ch, categories{k}(ikat));
                    end
                    % QUESTION - WHY IS THIS AVERAGING?
                    D = mean(dd, 3);
                else
                    D = obj.getaverageenvelopes(ch, categories(k));
                end
                subplot(1, numel(categories), k);
                T = obj.epochtime(1):0.1:obj.epochtime(2);
                imagesc(T, obj.Hfmean, D');
                maxy = max([maxy max(D(:))]); miny = min([miny min(D(:))]);
                axis xy;
                xlabel('time [s]');   
            end
            obj.plotF.ylim = [miny maxy];
        end
        
        function obj = plotlabels(obj, ch, categories)
            [miny, maxy] = obj.getplotylimit();
            for k = 1:numel(categories)
                subplot(1, numel(categories), k);
                caxis([miny, maxy]);
                title(obj.PsyData.CategoryName(cellval(categories, k)));
                if k == 1
                    chstr = iff(isempty(obj.CH.sortedby), num2str(ch), ...
                        [num2str(ch) '(' obj.CH.sortedby  num2str(obj.plotF.ch) ')']);
                    ylabel(['channel ' chstr ' - freq [Hz]']);
                    % TODO - Temporary fix for the multiple channel
                    % selection original:`any(obj.plotRCh.selCh(ch, :), 2) == 1`
                    % QUESTION - not sure what the plotRch.selCh does
                    if isprop(obj, 'plotRCh') && isfield(obj.plotRCh, 'selCh') && ...
                            any(any(obj.plotRCh.selCh(ch, :)))
                        klavesy = 'fghjkl';
                        text(0, obj.Hf(1), ['*' klavesy(logical(obj.plotRCh.selCh(ch, :)))], ...
                            'FontSize', 15, 'Color', 'red');
                    end
                end
                if k == numel(categories), colorbar('Position', [0.92 0.1 0.02 0.82]); end
            end
        end
                
        function [miny, maxy] = getplotylimit(obj)
             if isfield(obj.plotF, 'ylim') && numel(obj.plotF.ylim) >= 2
                miny = obj.plotF.ylim(1); maxy = obj.plotF.ylim(2);
             else
                miny = 0; maxy = 0;
             end
        end
        
        %% Getters
        % QUESTION - not sure why we add 1 to the categories 
        %
        % returns HFreqEpochs for given set of channels and categories. For
        % average call getaverageenvelopes
        % 
        % channels: array of channels. eg. [1, 5, 20]. If empty, returns
        % all channels. default []
        % categories: array of categories. eg. [1,3]. If empty, returns all
        % categories. Adds +1 to category number because of reasons - so if
        % you want category 1, you need to pass 0. default []
        % RETURN: matrix [time x channel x frequency x category]
        function envelopes = getenvelopes(obj, channels, categories)
            if ~exist('channels', 'var') || numel(channels) == 0
                channels = 1:size(obj.HFreqEpochs, 2);
            end
            if ~exist('categories', 'var') || numel(categories) == 0
                categories = 1:size(obj.HFreqEpochs, 4);
            else
                categories = obj.getcategoryindices(categories);
            end
            envelopes = obj.HFreqEpochs(:, channels, :, categories);
        end
        
        % averages envelopes per channels and categories
        % 
        % channels: vector(numeric) of channels to average across. If empty, returns
        % data only averaged across categories
        % categories: vector(numeric) of categories to average across. If
        % empty, returns
        % RETURN: matrix [time x frequency]
        function envelopes = getaverageenvelopes(obj, channels, categories)
            if ~exist('channels', 'var') || numel(channels) == 0
                channels = 1:size(obj.HFreq, 2);
            end
            if ~exist('categories', 'var') || numel(categories) == 0
                categories = 1:size(obj.HFreq, 4);
            else
                categories = obj.getaveragecategoryindices(categories);
            end
            envelopes = obj.HFreq(:, channels, :, categories);
            if numel(channels) >= 2, envelopes = mean(envelopes, 2); end
            if numel(categories) >= 2, envelopes = mean(envelopes, 4);end
            envelopes = squeeze(envelopes);
        end
        
        % Returns indices of given categories in the obj.HFreqEpochs
        % categories: vector of either characters cells or numbers defining
        % categories to be found int the obj.epochData. Zero based category
        % numbering
        % RETURN: indices of fitting categories
        % example: 
        %   obj.getcategoryindices([0 3])
        %   obj.getcategoryindices([{'Scene'} {'Object'}])
        function indices = getcategoryindices(obj, categories)
            switch class(categories)
                case 'double'
                    comparing = cellfun(@(x)x, obj.epochData(:, 2));
                case 'cell'
                    comparing = obj.epochData(:, 1);
                otherwise
                    return
            end
            indices = find(ismember(comparing, categories));
        end
        
        % Returns indices of categories as are in the obj.HFreq
        % categories: vector of either characters cells or numbers defining
        % categories to be found int the obj.epochData. Zero based category
        % numbering
        % example: 
        %   obj.getaveragecategoryindices([0 3])
        %   obj.getaveragecategoryindices([{'Scene'} {'Object'}])
        function indices = getaveragecategoryindices(obj, categories)
            conditions = obj.PsyData.P.strings.podminka;
            switch class(categories)
                case 'double'
                    indices = categories + 1;
                case 'cell'
                    iCategory = ismember(conditions(:,1), categories);
                    indices = cellfun(@(x)x + 1, conditions(iCategory, 2));
                otherwise
                    return
            end
        end
        
        % Returns indices in the envelope for given timewindow
        %
        % timewindow: numeric(2) with time in seconds
        % RETURNS: numeric(2) defining indices or [] if failed
        function indices = gettimewindowindices(obj, timewindow)
            indices = [];
            if numel(timewindow) ~= 2, return; end
            
            tick = (obj.epochtime(2) - obj.epochtime(1))/obj.fs;
            time = obj.epochtime(1) + (0:(size(obj.HFreq, 1) - 1)) * tick;
            % check if timewindow is in the epochtime boundaries
            if timewindow(1) >= max(time), return; end
            if timewindow(2) <= min(time), return; end
            
            indices = [find(time >= timewindow(1), 1) find(time <= timewindow(2), 1, 'last')];
        end
        
        %% Statistics
        % 
        % baselinetime: numeric(2) in seconds defining baseline timewindow
        % responsetime: numeric(2) in seconds defining response timewindow
        % RETURNS: calculated p values by CStat.Wilcox2D
        % example: obj.wilcoxbaseline([-0.1 0], [0.2 0.5])
        function wp = wilcoxbaseline(obj, baselinetime, responsetime, channels, categories)
            if ~exist('channels', 'var'), channels = []; end
            if ~exist('categories', 'var'), categories = []; end
            
            iResponse = obj.gettimewindowindices(responsetime);
            iBaseline = obj.gettimewindowindices(baselinetime);
            if any([numel(iBaseline) ~= 2, numel(iResponse) ~= 2]), return; end
            
            envelopes = obj.getenvelopes(channels, categories);
            response = envelopes(iResponse(1):iResponse(2), :, :, :);
            baseline = envelopes(iBaseline(1):iBaseline(2), :, :, :);
            
            wp = CStat.Wilcox2D(response, baseline, 1, [], 'mean vs baseline');
        end
        
        % Compares two category response in given time
        % responsetime: numeric(2) in seconds defining response timewindow.
        % defaults to the obj.epochtime
        % categories: numeric(2) or cell(2){character} defining categories.
        % Zero based. e.g [0 2]
        % RETURNS: calculated p values by CStat.Wilcox2D
        % example:
        %   obj.wilcoxcategories([0 1])
        %   obj.wilcoxcategories([{'Ovoce'} {'Scene'}], [0 0.5], 5:6)
        function wp = wilcoxcategories(obj, categories, responsetime, channels)
            if ~exist('responsetime', 'var') || numel(responsetime) == 0
                responsetime = obj.epochtime(1:2);
            end
            if ~exist('channels', 'var'), channels = []; end
            
            iResponse = obj.gettimewindowindices(responsetime);
            if any([numel(iResponse) ~= 2, numel(categories) ~= 2]), return; end
            
            responseA = obj.getenvelopes(channels, categories(1));
            responseA = responseA(iResponse(1):iResponse(2), :, :, :);
            responseB = obj.getenvelopes(channels, categories(2));
            responseB = responseB(iResponse(1):iResponse(2),:, :, :);
            
            wp = CStat.Wilcox2D(responseA, responseB, 1, [], 'mean vs baseline');
        end
        
        
        function obj = RemoveChannels(obj,channels)       
            %smaze se souboru vybrane kanaly. Kvuli redukci velikost aj
            keepch = setdiff(1:obj.channels,channels); %channels to keep
            obj.channels = obj.channels - numel(channels);
            obj.d = obj.d(:,keepch,:);
            if isprop(obj,'mults'), obj.mults = obj.mults(:,keepch); end
            obj.HFreq = obj.HFreq(:,keepch,:,:);
            if isprop(obj,'HFreqEpochs'), obj.HFreqEpochs = obj.HFreqEpochs(:,keepch,:,:); end
            if isprop(obj,'fphaseEpochs'), obj.fphaseEpochs = obj.fphaseEpochs(:,keepch,:,:); end
            if isprop(obj,'frealEpochs'), obj.frealEpochs = obj.frealEpochs(:,keepch,:,:); end            
            if isprop(obj,'plotRCh') && isfield(obj.plotRCh,'selCh')
               obj.plotRCh.selCh = obj.plotRCh.selCh(keepch,:); 
            end
            if isprop(obj,'RjEpochCh'), obj.RjEpochCh = obj.RjEpochCh(keepch,:); end
            obj.CH.RemoveChannels(channels);
            obj.els = obj.CH.els; %ty uz se redukuji v CHHeader
            obj.RjCh = obj.CH.RjCh;      
            
            for j = 1:numel(obj.Wp)
                obj.Wp(j).D2 = obj.Wp(j).D2(:,keepch);
                obj.Wp(j).DiEpCh = obj.Wp(j).DiEpCh(keepch,:);
                for k = 1:numel(obj.Wp(j).WpKat)
                    if numel(obj.Wp(j).WpKat{k}) > 0
                        obj.Wp(j).WpKat{k} = obj.Wp(j).WpKat{k}(:,keepch);
                    end
                end
                for k = 1:numel(obj.Wp(j).WpKatBaseline)
                    if numel(obj.Wp(j).WpKatBaseline{k}) > 0
                        obj.Wp(j).WpKatBaseline{k} = obj.Wp(j).WpKatBaseline{k}(:,keepch);
                    end
                end
            end
        end    
    end
    
    %% Visualisations
    methods  (Access = private)
        function obj = hybejPlotF(obj, ~, eventDat)
            if numel(obj.plotF.ch) == 1
                switch eventDat.Key
                   case 'rightarrow'
                       obj.PlotResponseFreq(min([obj.plotF.ch + 1 , obj.channels]));
                   case 'pagedown'
                       obj.PlotResponseFreq(min([obj.plotF.ch + 10 , obj.channels]));
                   case 'leftarrow'
                       obj.PlotResponseFreq(max([obj.plotF.ch - 1 , 1]));
                   case 'pageup'
                       obj.PlotResponseFreq(max([obj.plotF.ch - 10 , 1]));
                   case 'space' % zobrazi i prumerne krivky
                       obj.PlotResponseCh(obj.plotF.ch);
                       obj.PlotEpochs(obj.plotRCh.ch, obj.Wp(obj.WpActive).kats);
                       figure(obj.plotF.fh);
                end
            end
            switch eventDat.Key
                case {'multiply', '8'} % dialog na vlozeni minima a maxima osy y
                    answ = inputdlg('Enter ymax and min:', 'Yaxis limits', [1 50], {num2str(obj.plotF.ylim)});
                    if numel(answ) == 0, return; end
                    % answ cell 1x1 - cancel is cell 0x0
                    if isempty(answ{1}) || any(answ{1} == '*'), obj.plotF.ylim = [];
                    else
                        data = str2num(answ{:}); %#ok<ST2NM>
                        if numel(data) >= 2, obj.plotF.ylim = [data(1) data(2)]; end
                    end
                    obj.PlotResponseFreq(obj.plotF.ch);
                case {'divide', 'slash'} % automaticke meritko na ose z - power
                    obj.plotF.ylim = [];
                    obj.PlotResponseFreq(obj.plotF.ch);
                case {'add', 'equal','s'} % + oznaceni kanalu
                    obj.SelChannel(obj.plotF.ch);
                    obj.PlotResponseFreq(obj.plotF.ch);
                case {'numpad6', 'd'} % + oznaceni kanalu
                    ch2 = obj.plotRCh.selCh(find(obj.plotRCh.selCh > obj.plotF.ch, 1));
                    ch = iff(isempty(ch2), obj.plotF.ch, ch2);
                    obj.PlotResponseFreq(ch);
               % QUESTION - this keeps returning 0
               case {'numpad4', 'a'} % + oznaceni kanalu
                    ch2 = obj.plotRCh.selCh(find(obj.plotRCh.selCh < obj.plotF.ch, 1, 'last'));
                    ch = iff(isempty(ch2), obj.plotF.ch, ch2);
                    obj.PlotResponseFreq(ch);
            end
        end       
     end
end