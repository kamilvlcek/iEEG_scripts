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
        plotFrequency; % replaces plotF so both can be called seaprately
    end
    
    % -------------- public instance methods -------------------------
    methods (Access = public)
        %% Constructor description??
        % Just calls CHilbert class constructor
        function obj = CHilbertL(varargin)
            obj@CHilbert(varargin{:});
        end
        
        %% Plot Response Frequency
        
        function obj = plotresponsefrequency(obj, channels, categories)
        % Restructured PlotResponseFreq function which allows plotting of
        % either single or multiple channels averaged together.
        %
        % channels: numeric array definich which channels to plot. If more than
        %   one channel is passed, data are averaged across all given
        %   channels. e.g. 
        % categories: which category to plot - if not defined, plots all.
        %   Takes array of category indices beginning with 0. e.g. [0,3]. 
        %   See getenvelopes for better explanation
        %
        % Example:
        %   obj.plotresponsefrequency(1, [0:2])
        %   obj.plotresponsefrequency([1, 2:5, 8], [0:2])
        %   obj.plotresponsefrequency([1:5], 2)
            if ~exist('channels', 'var'), channels = []; end
            if ~exist('categories', 'var'), categories = []; end
           
            obj.prepareplotcategoriespowertime(channels, categories);
            % Shouldn't the Ch be obj.plotF.ch? The CHHeader has it as a 
            % field not a function and thus it is very difficult to unravel
            % The original function has the ch in plotF.ch, but uses the
            % obj.CH.sortorder(ch) to plot things - don't fully understand
            % if that is correct
            obj.plotcategoriespowertime(obj.CH.sortorder(channels), obj.plotFrequency.kategories);
            obj.plotlabels(obj.CH.sortorder(channels), obj.plotFrequency.kategories);
            set(obj.plotFrequency.fh, 'KeyPressFcn', @obj.frequencyplothandle);
        end
        
        function obj = PlotResponseFreqMean(obj, ch, categories)
        % Wrapper around plotresponsefrequency for plotting a response
        % frequencies for multiple channels
        % ch: numeric array definich which channels to plot. If more than
        %   one channel is passed, data are averaged across all given
        %   channels. e.g. 
        % categories: which category to plot - if not defined, plots all.
        %   Takes array of category indices beginning with 0. e.g. [0,3]. 
        %   See getenvelopes for better explanation
        %
        % Example:
        %   obj.PlotResponseFreqMean(1, [0:2])
        %   obj.PlotResponseFreqMean([1, 2:5, 8], [0:2])
        %   obj.PlotResponseFreqMean([1:5], 2)
            obj.plotresponsefrequency(ch, categories)
        end
        
        %% Getters
        
        % TODO! - rewrite to allow using RjEpCh a add Param parser probably
        function envelopes = getenvelopes(obj, channels, frequencies, categories)
        % returns HFreqEpochs for given set of channels and categories. For
        %   average call getaverageenvelopes
        % 
        % channels: array of channels. eg. [1, 5, 20]. If empty, returns
        %   all channels. default []
        % frequencies: array of frequencies or their indices. If floats are
        %   passed, then it is assumed its frequnecies. If integers are
        %   passed, then it evaluates as indices
        % categories: array of categories. eg. [1,3]. If empty, returns all
        % categories. Adds +1 to category number because of reasons - so if
        %   you want category 1, you need to pass 0. default []
        % RETURN: matrix [time x channel x frequency x category]
        % EXAMPLES:
        %   
            if ~exist('channels', 'var') || numel(channels) == 0
                channels = 1:size(obj.HFreqEpochs, 2);
            end
            if ~exist('frequencies', 'var') || numel(frequencies) == 0
                frequencies = 1:size(obj.HFreqEpochs, 3);
            else
                frequencies = obj.getfrequencyindices(frequencies);
            end
            if ~exist('categories', 'var') || numel(categories) == 0
                categories = 1:size(obj.HFreqEpochs, 4);
            else
                categories = obj.getcategoryindices(categories);
            end
            envelopes = obj.HFreqEpochs(:, channels, frequencies, categories);
        end
        
        function envelopes = getaverageenvelopes(obj, channels, frequencies, categories)
        % averages envelopes per channels and categories and 
        %   for particular frequencies. Uses data in the obj.HFreq, 
        %   selects based on channels and category indices and returns 
        %   a single averaged matrix
        % 
        % channels: vector(numeric) of channels to average across. If empty, 
        %   returns all channels averaged across categories
        % frequencies: array of frequencies or their indices. If floats are
        %   passed, then it is assumed its frequnecies. If integers are
        %   passed, then it evaluates as indices. See
        %   obj.getfrequencyindices
        % categories: vector(numeric) of categories to average across. If
        %   empty, averages across all categories
        % RETURN: matrix [time x frequency]
        % examples:
        %   % averages first 4 channels in the first category
        %   hilbert.getaverageenvelopes([0:3], [], 0)
        %   % averages all channels in the first category
        %   hilbert.getaverageenvelopes([], [], 1)
        %   % averages all channels across all categories
        %   hilbert.getaverageenvelopes([], [], [])
            if ~exist('channels', 'var') || numel(channels) == 0
                channels = 1:size(obj.HFreq, 2);
            end
            if ~exist('frequencies', 'var') || numel(frequencies) == 0
                frequencies = 1:size(obj.HFreq, 3);
            else
                frequencies = obj.getfrequencyindices(frequencies);
            end
            if ~exist('categories', 'var') || numel(categories) == 0
                categories = 1:size(obj.HFreq, 4);
            else
                categories = obj.getaveragecategoryindices(categories);
            end
            envelopes = obj.HFreq(:, channels, frequencies, categories);
            if numel(channels) >= 2, envelopes = mean(envelopes, 2); end
            if numel(categories) >= 2, envelopes = mean(envelopes, 4);end
            envelopes = squeeze(envelopes);
        end
        
        function indices = getcategoryindices(obj, categories)
        % Returns indices of given categories in the obj.HFreqEpochs
        % categories: vector of either characters cells or numbers defining
        % categories to be found int the obj.epochData. Zero based category
        %   numbering
        % RETURN: indices of fitting categories
        % example: 
        %   obj.getcategoryindices([0 3])
        %   obj.getcategoryindices([{'Scene'} {'Object'}])
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

        function indices = getaveragecategoryindices(obj, categories)
        % Returns indices of categories as are in the obj.HFreq
        % categories: vector of either characters cells or numbers defining
        % categories to be found int the obj.epochData. Zero based category
        %   numbering
        % example: 
        %   obj.getaveragecategoryindices([0 3])
        %   obj.getaveragecategoryindices([{'Scene'} {'Object'}])
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
        
        function indices = gettimewindowindices(obj, timewindow)
        % Returns indices in the envelope for given timewindow
        % 
        % timewindow: numeric(2) with time in seconds
        % RETURNS: numeric(2) defining indices or [] if failed
            indices = [];
            if numel(timewindow) ~= 2, return; end
            
            tick = (obj.epochtime(2) - obj.epochtime(1))/obj.fs;
            time = obj.epochtime(1) + (0:(size(obj.HFreq, 1) - 1)) * tick;
            % check if timewindow is in the epochtime boundaries
            if timewindow(1) >= max(time), return; end
            if timewindow(2) <= min(time), return; end
            
            indices = [find(time >= timewindow(1), 1) find(time <= timewindow(2), 1, 'last')];
        end
          
        function indices = getfrequencyindices(obj, frequencies)
        % Returns indices of freqneuncies oif passed as float array. Or
        %   returns indices back if passed as integer array
        %
        % frequencies: array of frequencies or their indices. If doubles are
        %   passed, then it is assumed its frequnecies. If integers/whole numbers
        %   are passed, then they are evaluated as indices
            if(~all(frequencies == floor(frequencies))) %checks if doubles
                % selects only those frequencies which actually are in the
                % Hfmean frequency vector
                frequencies = frequencies(ismember(frequencies, obj.Hfmean));
                indices = arrayfun(@(x) find(obj.Hfmean == x), frequencies);
            else % returns indices
                indices = frequencies;
            end
        end
        %% Statistics

        function wp = wilcoxbaseline(obj, baselinetime, responsetime, frequencies, categories)
        % Calculates rank test for passed categories and frequencies against a baseline.
        %
        % Compares epochs during the response time to the epochs during to
        % baselineteime in each channel, frequency and category. Returns a
        % matrix with time x channel x frequency x category
        % 
        % baselinetime: numeric(2) in seconds defining baseline timewindow
        % responsetime: numeric(2) in seconds defining response timewindow
        % frequencies: array of frequencies to analyse. if [], runs for all
        %   frequencies. See obj.getenvelopes for description
        % categories: list of cells of wanted categories. e.g. {'Ovoce'} or
        %   category numbers as they are in the psychodata (starting with 0)
        % RETURNS: a 4d matrix calculated p values by CStat.Wilcox2D.  
        %   returned matrix is p value in time x channel x frequency x category
        %   Non calculated comparisons are left empty with NaNs.    
        %   e.g. result(:, 1, 1, 2) is a array of p values for all times
        %   for first channel, first frequency, and second category
        %   if the results doesn't have a category 2 calculated, retunrs an
        %   array of NaNs
        % example: 
        %   obj.wilcoxbaseline([-0.1 0], [0.2 0.5])
            iResponse = obj.gettimewindowindices(responsetime);
            iBaseline = obj.gettimewindowindices(baselinetime);
            if any([numel(iBaseline) ~= 2, numel(iResponse) ~= 2]), return; end
            
            % returned matrix is time x channel x frequency x category comparison
            wpSize = [diff(iResponse) + 1, ...
                size(obj.HFreqEpochs, 2), ...
                size(obj.HFreqEpochs, 3), ...
                size(unique(obj.epochData(:,1)), 1)];
            wp = NaN(wpSize);
            
            if ~exist('frequencies', 'var') || numel(frequencies) == 0
                frequencies = 1:wpSize(3);
            end
            if ~exist('categories', 'var') || numel(categories) == 0
                categories = 1:wpSize(4);
            end
            
            % separate comparison for each channel, frequency and category
            for iFrequency = 1:numel(frequencies)
                fprintf('Comparing for frequency %f\n', obj.Hfmean(iFrequency));
                for category = categories
                    if(iscell(category))
                        iCategory = obj.PsyData.CategoryNum(category);
                    else
                        iCategory = category;
                    end
                    [~, ~, RjEpCh] = obj.CategoryData(iCategory);
                    
                    envelopes = obj.getenvelopes([], frequencies(iFrequency), iCategory);
                    response = squeeze(envelopes(iResponse(1):iResponse(2), :, :, :));
                    baseline = squeeze(envelopes(iBaseline(1):iBaseline(2), :, :, :));
                    % As we are passing the same set of epochs, we can
                    % basically "clone" the epochs?
                    tempWp = CStat.Wilcox2D(response, baseline, 0, [], '', ...
                        RjEpCh, RjEpCh);
                    wp(:, :, iFrequency, iCategory + 1) = tempWp;
                end
            end
        end

        %%% comparing sets of channels instead of epcohs within a channel
        % against each other
        
        function wp = wilcoxcategories(obj, categories, responsetime, frequencies)
        % Compares two category responses in given time
        % 
        % categories: numeric(2) or cell(2){character} defining categories.
        %   Zero based. e.g [0 2]
        % responsetime: numeric(2) in seconds defining response timewindow.
        %   defaults to the obj.epochtime
        % RETURNS: a 3d matrix with calculated p values by CStat.Wilcox2D
        %   matrix is time x channel x frequency
        % example:
        %   obj.wilcoxcategories([0 1])
        %   obj.wilcoxcategories([{'Ovoce'} {'Scene'}], [0 0.5], 5:6)
            if ~exist('responsetime', 'var') || numel(responsetime) == 0
                responsetime = obj.epochtime(1:2);
            end
            iResponse = obj.gettimewindowindices(responsetime);
            if any([numel(iResponse) ~= 2, numel(categories) ~= 2]), return; end
            
            wpSize = [diff(iResponse) + 1, ...
                size(obj.HFreqEpochs, 2), ...
                size(obj.HFreqEpochs, 3)];
            wp = NaN(wpSize);
            
            if ~exist('frequencies', 'var') || numel(frequencies) == 0
                frequencies = 1:wpSize(3);
            end
            
            rejectedEpochs = {2};
            for iCategory = 1:2
                [~, ~, RjEpCh] = obj.CategoryData(categories{(iCategory)});
                rejectedEpochs(iCategory) = {RjEpCh};
            end
            
            for iFrequency = 1:numel(frequencies)
                fprintf('Comparing for frequency %f\n', obj.Hfmean(iFrequency));
                responseA = obj.getenvelopes([], frequencies(iFrequency), categories(1));
                responseA = squeeze(responseA(iResponse(1):iResponse(2), :, :, :));
                responseB = obj.getenvelopes([], frequencies(iFrequency), categories(2));
                responseB = squeeze(responseB(iResponse(1):iResponse(2), :, :, :));
                % As we are passing the same set of epochs, we can
                % basically "clone" the epochs?
                tempWp = CStat.Wilcox2D(responseA, responseB, 0, [], '', ...
                    rejectedEpochs{(1)}, rejectedEpochs{(2)});
                wp(:, :, iFrequency) = tempWp;
            end
        end
        
        function wp = wilcoxaveragebaseline(obj, baselinetime, responsetime, channels1, channels2, frequencies, categories)
        % Compares sets of channels against each other. for given
        % frequencies and categories
        % Runs a wilcox test for an average across all non rejected epochs
        % 
        % categories: 
        % channels: two dimensional matrix of channels to average and run 
        % responsetime:
        % frequencies
        %
        
        % Get the epoch data for given parameters and only non rejected
        
        % [~, ~, RjEpCh] = obj.CategoryData(iCategory);
        % Average non rejected epochs across channels
        % - resp: time x channel x frequency - cpomaring sets of channels
        % against each other
            
        
            wp = [];
        end
        
        function wp = wilcoxaveragecategories(obj, categories, responsetime, channels1, channels2, frequencies)
        % Runs a wilcox test for an average across epochs and
        % compares them betweeen passed categories
        % 
        % categories: two categories to be compared
        % channels: two dimensional matrix of channels to average and run 
        % responsetime:
        % frequencies
        % Notes: Basically runs like wilcoxchannelcategories but replaces
        % Validate only two categories
        if iscell(categories)
            categories = obj.PsyData.CategoryNum(categories);
        end
        [~, ~, ~, iEpochs1] = obj.CategoryData(categories(1));
        [~, ~, ~, iEpochs2] = obj.CategoryData(categories(2));
        % Filter epochs for given data
        category1 = obj.HFreqEpochs;
        % Filter given channels
            wp = [];
        end
    end  
    
    %% Visualisations
    methods  (Access = private)
        % Buffering of the plot, loading settings from saved �n the 
        % obj.plotFrequency. Returns existing plot if it exists
        function obj = prepareplotcategoriespowertime(obj, ch, categories)
            if ~numel(ch) == 0 && ~isfield(obj.plotFrequency, 'ch'), obj.plotFrequency.ch = 1;
            else, obj.plotFrequency.ch = ch; end
            
            % Only rewrites categorties if not already set
            if numel(categories) > 0, obj.plotFrequency.kategories = categories;
            elseif ~isfield(obj.plotFrequency, 'kategories') && numel(obj.plotFrequency.kategories) == 0
                if isfield(obj.Wp(obj.WpActive), 'kats'), categories = obj.Wp(obj.WpActive).kats;
                else, categories = obj.PsyData.Categories(); end
                obj.plotFrequency.kategories = categories;
            end
            
            % redraws the plot or creates one
            if isfield(obj.plotFrequency, 'fh') && ishandle(obj.plotFrequency.fh)
                figure(obj.plotFrequency.fh);
            else
                obj.plotFrequency.fh = figure('Name', 'ResponseFreq', 'Position', [20, 500, 1200, 300]);
                colormap jet;
            end
        end
        
        % TODO - deprecate option to pass cells into categories
        % Actual plotting function. Takes channels and categories and
        % creates imagesc based on averaged data (averaged across channels
        % and categories passed)
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
                        dd(:, :, ikat) = obj.getaverageenvelopes(ch, [], categories{k}(ikat));
                    end
                    % QUESTION - WHY IS THIS AVERAGING?
                    D = mean(dd, 3);
                else
                    D = obj.getaverageenvelopes(ch, [], categories(k));
                end
                subplot(1, numel(categories), k);
                T = obj.epochtime(1):0.1:obj.epochtime(2);
                imagesc(T, obj.Hfmean, D');
                axis xy;
                xlabel('time [s]');
                % Buffers the x and y limits to use in future plots
                maxy = max([maxy max(D(:))]); miny = min([miny min(D(:))]);
            end
            obj.plotFrequency.ylim = [miny maxy];
        end
        
        %Adds labels to the plotfrequency plot
        function obj = plotlabels(obj, ch, categories)
            [miny, maxy] = obj.getplotylimit();
            for k = 1:numel(categories)
                subplot(1, numel(categories), k);
                caxis([miny, maxy]);
                %TODO - make this a function, this can change in the future
                %or error out
                title(obj.PsyData.CategoryName(cellval(categories, k)));
                if k == 1
                    chstr = iff(isempty(obj.CH.sortedby), num2str(ch), ...
                        [num2str(ch) '(' obj.CH.sortedby  num2str(obj.plotFrequency.ch) ')']);
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
          
        % Gets plot y limits, either existing or defaults to [0,0]
        function [miny, maxy] = getplotylimit(obj)
             if isfield(obj.plotFrequency, 'ylim') && numel(obj.plotFrequency.ylim) >= 2
                miny = obj.plotFrequency.ylim(1); maxy = obj.plotFrequency.ylim(2);
             else
                miny = 0; maxy = 0;
             end
        end
        
        % Replaces hybejPlotF to separate functionality of parent CHilbert
        % and CHilbertL
        function obj = frequencyplothandle(obj, ~, eventDat)
            if numel(obj.plotFrequency.ch) == 1
                switch eventDat.Key
                   case 'rightarrow'
                       obj.plotresponsefrequency(min([obj.plotFrequency.ch + 1 , obj.channels]));
                   case 'pagedown'
                       obj.plotresponsefrequency(min([obj.plotFrequency.ch + 10 , obj.channels]));
                   case 'leftarrow'
                       obj.plotresponsefrequency(max([obj.plotFrequency.ch - 1 , 1]));
                   case 'pageup'
                       obj.plotresponsefrequency(max([obj.plotFrequency.ch - 10 , 1]));
                   case 'space' % zobrazi i prumerne krivky
                       obj.PlotResponseCh(obj.plotFrequency.ch);
                       obj.PlotEpochs(obj.plotRCh.ch, obj.Wp(obj.WpActive).kats);
                       figure(obj.plotFrequency.fh);
                end
            end
            switch eventDat.Key
                case {'multiply', '8'} % dialog na vlozeni minima a maxima osy y
                    answ = inputdlg('Enter ymax and min:', 'Yaxis limits', [1 50], {num2str(obj.plotFrequency.ylim)});
                    if numel(answ) == 0, return; end
                    % answ cell 1x1 - cancel is cell 0x0
                    if isempty(answ{1}) || any(answ{1} == '*'), obj.plotFrequency.ylim = [];
                    else
                        data = str2num(answ{:}); %#ok<ST2NM>
                        if numel(data) >= 2, obj.plotFrequency.ylim = [data(1) data(2)]; end
                    end
                    obj.plotresponsefrequency(obj.plotFrequency.ch);
                case {'divide', 'slash'} % automaticke meritko na ose z - power
                    obj.plotFrequency.ylim = [];
                    obj.plotresponsefrequency(obj.plotFrequency.ch);
                case {'add', 'equal','s'} % + oznaceni kanalu
                    obj.SelChannel(obj.plotFrequency.ch);
                    obj.plotresponsefrequency(obj.plotFrequency.ch);
                case {'numpad6', 'd'} % + oznaceni kanalu
                    ch2 = obj.plotRCh.selCh(find(obj.plotRCh.selCh > obj.plotFrequency.ch, 1));
                    ch = iff(isempty(ch2), obj.plotFrequency.ch, ch2);
                    obj.plotresponsefrequency(ch);
               % QUESTION - this keeps returning 0
               case {'numpad4', 'a'} % + oznaceni kanalu
                    ch2 = obj.plotRCh.selCh(find(obj.plotRCh.selCh < obj.plotFrequency.ch, 1, 'last'));
                    ch = iff(isempty(ch2), obj.plotFrequency.ch, ch2);
                    obj.plotresponsefrequency(ch);
            end
        end       
     end
end