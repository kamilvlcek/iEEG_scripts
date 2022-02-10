classdef CHilbertL < CHilbert
    % HILBERT.CLASS extension for CHilbert class
    % Original author: Lukáš Hejtmánek
    
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
    
    %% -------------- public instance methods -------------------------
    methods (Access = public)
        %% Constructor description??
        % Just calls CHilbert class constructor
        function obj = CHilbertL(varargin)
            obj@CHilbert(varargin{:});
        end
        
        %% Plot Response Frequency
        
        function obj = plotresponsefrequency(obj, channels, varargin)
        % Restructured PlotResponseFreq function which allows plotting of
        % either single or multiple channels averaged together.
        %
        % channels: numeric array definich which channels to plot. If more than
        %   one channel is passed, data are averaged across all given
        %   channels. e.g. 
        % varargin
        % categories: which category to plot - if not defined, plots all.
        %   Takes array of category indices beginning with 0. e.g. [0,3]. 
        %   See getenvelopes for better explanation
        %
        % Example:
        %   obj.plotresponsefrequency(1, [0:2])
        %   obj.plotresponsefrequency([1, 2:5, 8], 'categories', [0:2])
        %   obj.plotresponsefrequency([1:5], 'categories', 2)
            p = inputParser;
            p = obj.addcategoriesparameter(p);
            addParameter(p, 'ylim', obj.getplotylimit, ...
                @(x)(numel(x) == 2 && isnumeric(x)));
            parse(p, varargin{:}); %p.Parameters contains names of arguments and p.Results contains their values
            
            if ~exist('channels', 'var'), channels = []; end
            obj.plotFrequency = struct;
            obj.prepareplotcategoriespowertime(channels, p.Results.categories); %creates figure, stores args in obj.plotFrequency
            obj.plotFrequency.ylim = p.Results.ylim;
            obj.updateplotcategoriespowertime;
        end
        
        function obj = PlotResponseFreqMean(obj, categories, zlim, ch)
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
            if ~exist('ch','var') || isempty(ch)
                ch = obj.CH.sortorder; %currently filtered channels
            end %ch is later used as index in obj.CH.sortorder
            if ~exist('categories','var') || isempty(categories)                
                categories = obj.Wp(obj.WpActive).kats;
            elseif ~isempty(obj.Wp) && isfield(obj.Wp,'kats')
                katnum = obj.Wp(obj.WpActive).kats(categories); %use numbers from currect contrast if available                
                categories = katnum;                
            end
            if exist('zlim','var') && numel(zlim) == 2                
                obj.plotresponsefrequency(ch, 'categories', categories,'ylim', zlim);
            else
                obj.plotresponsefrequency(ch, 'categories', categories);
            end
            
        end
        
        %% Getters
        
        function envelopes = getenvelopes(obj, varargin)
        % returns HFreqEpochs for given set of channels, frequencies, and
        %   categories. For average call getaverageenvelopes. Contains
        % 
        % channels: array of channels. eg. [1, 5, 20]. If empty, returns
        %   all channels. default []
        % frequencies: array of frequencies or their indices. If floats are
        %   passed, then it is assumed its frequnecies. If integers are
        %   passed, then it evaluates as indices
        % categories: array of categories. eg. [1,3]. If empty, returns all
        %   categories. Adds +1 to category number because of reasons - so 
        %   if you want category 1, you need to pass 0. default []
        % time: time in seconds to extract
        % reject: logical defining if the rejected epochs should be
        %   excluded. You can ONLY reject epochs if you are selecting a
        %   single category , as epochs are rejected per channel x category 
        %   and trying to return multiple would result in an uneven matrix. 
        %   If multiple channels are selected in the channels parameter,
        %   then can only return if all have the same number of final
        %   epochs. Otherwise errors out.
        %   Defaults to false.
        % RETURN: 4D matrix [time x channel x frequency x epoch]
        % EXAMPLES: 
        %   hilbert.getenvelopes('channels', [1 3])
        %   hilbert.getenvelopes('channels', 1, 'frequencies', [52.75])
        %   hilbert.getenvelopes('frequencies', [1 3], 'categories', {'Ovoce'})
            p = inputParser;
            p = obj.addchannelsparameter(p);
            p = obj.addfrequenciesparameter(p);
            p = obj.addcategoriesparameter(p);
            addParameter(p, 'time', obj.epochtime(1:2),@(x)(numel(x) == 2 && ...
                x(1) >= obj.epochtime(1) && x(1) <= obj.epochtime(2)));
            addParameter(p, 'reject', false, @islogical);
            parse(p, varargin{:});
            
            frequencies = obj.getfrequencyindices(p.Results.frequencies);
            iCategories = obj.getcategoryindices(p.Results.categories);
            iTime = obj.gettimewindowindices(p.Results.time);
 
            if p.Results.reject
            % NOT WORKING because I am not sure how to implement it
            % This only works if there are no per channel specific epoch
            % rejections - e.g. all channels have the same epoch removed.
            % Otherwise the returned matrix would be uneven - e.g. 51
            % epochs in 4th dimention for channel 1 but only 50 for channel
            % 2. This could be solved by inserting NaNs but if feels as odd
            % behaviour
                if numel(p.Results.categories) ~= 1
                    error(['Can only reject epochs for a single category' ...
                        ' and a single channel']);
                end
                category = obj.categorytonum(p.Results.categories);
                [~, ~, ~, iEpochs] = obj.CategoryData(category(1));
                % we are only interested in the subset of channels
                iEpochs = iEpochs(:, p.Results.channels); 
                if ~all(sum(iEpochs) == sum(iEpochs(:, 1)))
                    error(['Not all channels have the same number of' ...
                        ' epochs rejected, cannot return.']);
                else
                    iCategories = iEpochs(:,1);
                end
            end
            envelopes = obj.HFreqEpochs(iTime(1):iTime(2),p.Results.channels,...
                frequencies,iCategories);
        end
        
        function envelopes = getaverageenvelopes(obj, varargin)
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
            p = inputParser;
            p = obj.addchannelsparameter(p);
            p = obj.addfrequenciesparameter(p);
            p = obj.addcategoriesparameter(p);
            parse(p, varargin{:});
            
            frequencies = obj.getfrequencyindices(p.Results.frequencies);
            categories = obj.getaveragecategoryindices(p.Results.categories);

            envelopes = obj.HFreq(:, p.Results.channels, frequencies, categories); %all samples, selected channels, frequencies and categories
            
            if numel(p.Results.channels) >= 2, envelopes = mean(envelopes, 2); end %mean over channels
            if numel(categories) >= 2, envelopes = mean(envelopes, 4);end %mean over categories
            envelopes = squeeze(envelopes);
        end
        
        % TODO - use existing psyData functions
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
        % TODO - the fdr correction switch
        function wp = wilcoxbaseline(obj, varargin)
        % Calculates rank test for power against a baseline. 
        %
        % Compares epochs during the response time to the epochs during to
        %   baselineteime in each channel, frequency and category separately. 
        %   Returns a matrix with time x channel x frequency x category
        % OPTIONAL PARAMETERS
        %   baseline: numeric(2) in seconds defining baseline timewindow
        %   response: numeric(2) in seconds defining response timewindow
        %   frequencies: array of frequencies to analyse. See 
        %       obj.getenvelopes for description
        %   categories: array of wanted categories. See obj.getenvelopes for description
        %   channels: array of wanted channels. See obj.getenvelopes for description
        %   fdrMethod: if the fdr Correction shhould be applied. Correction
        %       is applied for each channel and category separately (across all
        %       calculated frequencies and time). Values can be either 
        %       [] (no correction - default), pdep or dep. 
        %   squeeze: if true, only selects calculated values, not returning NaNs
        % RETURNS: a 4d matrix calculated p values by CStat.Wilcox2D. 
        %   returned matrix is p value in time x channel x frequency x category
        %   Non calculated comparisons are left empty with NaNs, or
        %   squeezed.
        %   e.g. result(:, 1, 1, 2) is a array of p values for all times
        %   for first channel, first frequency, and second category
        %   if the results doesn't have a category 2 calculated, retunrs an
        %   array of NaNs
        % example: 
        %   obj.wilcoxbaseline('baselinetime', [-0.1 0], 'responsetime', [0.2 0.5])
            p = inputParser;
            p = obj.addchannelsparameter(p);
            p = obj.addfrequenciesparameter(p);
            p = obj.addcategoriesparameter(p);
            p = obj.addfdrparameter(p);
            addParameter(p, 'baseline', obj.baseline, @(x)(numel(x) == 2));
            addParameter(p, 'response', [obj.baseline(2) obj.epochtime(2)],...
                @(x)(numel(x) == 2));
            addParameter(p, 'squeeze', false, @islogical);
            parse(p, varargin{:});
            
            iResponse = obj.gettimewindowindices(p.Results.response);
            iBaseline = obj.gettimewindowindices(p.Results.baseline);
            if any([numel(iBaseline) ~= 2, numel(iResponse) ~= 2]), return; end
            
            % returned matrix is time x channel x frequency x category comparison
            wp = NaN(diff(iResponse) + 1,...
                size(obj.HFreqEpochs, 2),...
                size(obj.HFreqEpochs, 3),...
                size(obj.PsyData.Categories(false), 2));
            
            iCategories = obj.categorytonum(p.Results.categories);
            iFrequencies = obj.getfrequencyindices(p.Results.frequencies);
            % separate comparison for each channel, frequency and category
            for iFrequency = iFrequencies
                fprintf('Comparing for frequency %f\n', obj.Hfmean(iFrequency));
                for iCategory = iCategories
                    [~, ~, RjEpCh] = obj.CategoryData(iCategory);
                    envelopes = obj.getenvelopes('channels', p.Results.channels,...
                        'frequencies', iFrequency, 'categories', iCategory);
                    % selects times and squeezes the categories and frequency
                    response = squeeze(envelopes(iResponse(1):iResponse(2), :, :, :));
                    
                    %% This should be probably averaged in the first dim
                    baseline = squeeze(envelopes(iBaseline(1):iBaseline(2), :, :, :));
                    baseline = mean(baseline, 1);
                    % NO FDR CORRECTION 
                    tempWp = CStat.Wilcox2D(response, baseline, 0, 0, '',...
-                        RjEpCh, RjEpCh);
                    wp(:, p.Results.channels, iFrequency, iCategory + 1) = tempWp;
                end
            end
            % FDR correction for each channel and category separately
            if ~isempty(p.Results.fdrMethod)
                message(['Performing fdr correction ' p.Results.fdrMethod]);
                for iCategory = iCategories
                    for channel = p.Results.channels
                        tempWp = wp(:, channel, :, iCategory + 1);
                        [~, ~, adj_p] = fdr_bh(tempWp, 0.05, p.Results.fdrMethod, 'no');
                        wp(:, channel, :, iCategory + 1) = adj_p;
                    end
                end
            end
            % squeezing
            if p.Results.squeeze
                wp = wp(:, p.Results.channels, iFrequencies, iCategories + 1);
            end
        end

        % TODO - the fdr correction switch - currently run per frequnecy
        function wp = wilcoxcategories(obj, categories, varargin)
        % Compares two category responses in given time
        % 
        % REQUIRED:
        %   categories: two categories to compare. see getenvelopes on how
        %   to pass the list
        % OPTIONAL PARAMETERS:
        %   channels: array of channels to select. see getenvelopes
        %   frequencies: array of frequencies to select. see getenvelopes
        %   response: numeric(2) in seconds defining response timewindow.
        %   fdrMethod: if the fdr correction shhould be applied. Correction
        %       is applied for each channel separately (across all
        %       calculated frequencies and time). Values can be either 
        %       [] (no correction - default), pdep or dep. 
        %   squeeze: should the non calculated values be dropped?
        % RETURNS: a 3d matrix with calculated p values by CStat.Wilcox2D
        %   matrix is time x channel x frequency
        % example:
        %   obj.wilcoxcategories([0 1])
        %   obj.wilcoxcategories([{'Ovoce'} {'Scene'}], [0 0.5], 5:6)
            p = inputParser;
            addRequired(p, 'categories', @(x)(numel(x) == 2));
            p = obj.addchannelsparameter(p);
            p = obj.addfrequenciesparameter(p);
            p = obj.addresponseparameter(p);
            p = obj.addfdrparameter(p);
            addParameter(p, 'squeeze', false, @islogical);
            parse(p, categories, varargin{:});
            
            iResponse = obj.gettimewindowindices(p.Results.response);
            wp = NaN([diff(iResponse) + 1, size(obj.HFreqEpochs, 2),...
                size(obj.HFreqEpochs, 3)]);
            
            rejectedEpochs = cell(2);
            for iCategory = 1:2
                [~, ~, RjEpCh] = obj.CategoryData(categories{(iCategory)});
                rejectedEpochs(iCategory) = {RjEpCh};
            end
            
            iFrequencies = obj.getfrequencyindices(p.Results.frequencies); 
            for iFrequency = iFrequencies
                fprintf('Comparing for frequency %f\n', obj.Hfmean(iFrequency));
                responseA = obj.getenvelopes('frequencies', iFrequency,...
                    'categories', categories(1));
                responseA = squeeze(responseA(iResponse(1):iResponse(2), :, :, :));
                responseB = obj.getenvelopes('frequencies', iFrequency,...
                    'categories', categories(2));
                responseB = squeeze(responseB(iResponse(1):iResponse(2), :, :, :));
                % As we are passing the same set of epochs, we can
                % basically "clone" the epochs?
                % NO FDR at this point
                tempWp = CStat.Wilcox2D(responseA, responseB, 0, 0, '',...
                    rejectedEpochs{(1)}, rejectedEpochs{(2)});
                wp(:, :, iFrequency) = tempWp;
            end
            if ~isempty(p.Results.fdrMethod)
                message(['Performing fdr correction ' p.Results.fdrMethod]);
                for channel = p.Results.channels
                    tempWp = wp(:, channel, :);
                    [~, ~, adj_p] = fdr_bh(tempWp, 0.05, p.Results.fdrMethod, 'no');
                    wp(:, channel, :) = adj_p;
                end
            end
            % squeezing
            if p.Results.squeeze
                wp = wp(:, p.Results.channels, iFrequencies);
            end
        end
        
        %%% comparing sets of channels instead of epochs within a channel
        % against each other
        
        function wp = wilcoxaveragebaseline(obj, varargin)
        % Compares averages of channels against each other
        %   Runs a wilcox test for an average across all non rejected epochs
        % OPTIONAL PARAMETERS:
        %   channels: array of channels to select. see getenvelopes
        %   frequencies: array of frequencies to select. see getenvelopes
        %   categories: array of categories to select. see getenvelopes
        %   baseline: numeric(2) in seconds defining baseline timewindow.
        %   response: numeric(2) in seconds defining response timewindow.
        %   squeeze: should the non calculated values be dropped?
        % RETURNS:
        %   returns 3d matrix time x frequency x category with p values and
        %   NaNs where data was not required. If squeeze is true, then only
        %   returns squeezed amtrix with only selected values should be
        %   returned
            p = inputParser;
            p = obj.addchannelsparameter(p);
            p = obj.addfrequenciesparameter(p);
            p = obj.addcategoriesparameter(p);
            p = obj.addresponseparameter(p);
            addParameter(p, 'baseline', obj.baseline, @(x)(numel(x) == 2));
            addParameter(p, 'fdrMethod', [],...
                @(x)isstring(x) && any(strcmp({'pdep' 'dep'},x)));
            addParameter(p, 'squeeze', false, @(x)(numel(x) == 1 && islogical(x)));
            parse(p, varargin{:});
            
            channels = p.Results.channels;
            
            iCategories = obj.categorytonum(p.Results.categories);
            iFrequencies = obj.getfrequencyindices(p.Results.frequencies);
            
            % creates a buffer so we don't change the matrix size in the
            % loop -  3d matrix time x frequency x category
            wp = NaN(diff(obj.gettimewindowindices(p.Results.response)) + 1,...
                obj.nfrequencies, obj.ncategories);
            for iCategory = iCategories
                disp(['Calculating for category ', obj.PsyData.CategoryName(iCategory)]);
                % Prepares empty matrices. Each is time x frequency x channel
                baseline = NaN(1, obj.nfrequencies, numel(channels));
                response = NaN(size(wp, 1), obj.nfrequencies, numel(channels));
                for iFrequency = iFrequencies
                    % because it is possible that each channel has a
                    % different number of rejected epochs, we have to
                    % select and average per channel :( TODO - when the
                    % getenvelopes returns in case of unequal rejections,
                    % then this can be heavilly streamlined
                    for channel = channels
                        tempBaseline = obj.getenvelopes('channels',channel,...
                            'frequencies',iFrequency,'categories',iCategory,...
                            'time',p.Results.baseline,'reject',true);
                        tempResponse = obj.getenvelopes('channels',channel,...
                            'frequencies',iFrequency,'categories',iCategory,...
                            'time',p.Results.response, 'reject',true);
                        % averaging across epochs
                        response(:, iFrequency, channel) = mean(tempResponse, 4);
                        tempBaseline = mean(tempBaseline, 4);
                        %averaging baseline time to a single measurement
                        baseline(1, iFrequency, channel) = mean(tempBaseline, 1);
                    end
                end
                % removes the NaNs from the baseline and response
                baseline = baseline(:, iFrequencies, channels);
                response = response(:, iFrequencies, channels);
                tempWp = wilcox3d(response, baseline);
                wp(:, iFrequencies, iCategory + 1) = tempWp;
            end
           
            if ~p.Results.squeeze
                wp = wp(:, iFrequencies, iCategories + 1);
            end
        end
        
        function wp = wilcoxaveragecategories(obj, categories, varargin)
        % Runs a wilcox test for an average resposne across epochs and 
        %   compares them betweeen passed categories
        % REQUIRED PARAMETERS:
        %   categories: two categories to be compared
        % OPTIONAL PARAMETERS  
        %   channels: array of channels to select. see getenvelopes
        %   frequencies: array of frequencies to select. see getenvelopes
        %   response: numeric(2) in seconds defining response timewindow.
        %   squeeze: should the non calculated values be dropped?
        % RETURNS
        %   2d matrix with p values with time x frequency. NaNs at positions 
        %   where the calculation was not required. If squeeze is true, then only
        %   returns squeezed amtrix with only selected values should be returned
            p = inputParser;
            addRequired(p, 'categories', @(x)(numel(x) > 0));
            p = obj.addchannelsparameter(p);
            p = obj.addfrequenciesparameter(p);
            p = obj.addresponseparameter(p);
            addParameter(p, 'squeeze', false, @(x) islogical(x));
            parse(p, categories, varargin{:});

            iCategories = obj.categorytonum(p.Results.categories);
            iFrequencies = obj.getfrequencyindices(p.Results.frequencies);
            channels = p.Results.channels;

            % default matrix is time x single frequency x channel
            category1 = NaN(diff(obj.gettimewindowindices(p.Results.response)) + 1,...
                obj.nfrequencies, numel(channels));
            category2 = category1;
            for iFrequency = iFrequencies
                for channel = channels
                    % because it is possible that each channel has a
                    % different number of rejected epochs, we have to
                    % select and average per channel :( TODO - when the
                    % getenvelopes returns in case of unequal rejections,
                    % then this can be heavilly streamlined
                    tempCategory1 = obj.getenvelopes('channels',channel,...
                        'frequencies',iFrequency,'categories',iCategories(1),...
                        'time',p.Results.response,'reject',true);
                    tempCategory2 = obj.getenvelopes('channels',channel,...
                        'frequencies',iFrequency,'categories',iCategories(2),...
                        'time',p.Results.response,'reject',true);
                    category1(:, iFrequency, channel) = mean(tempCategory1, 4);
                    category2(:, iFrequency, channel) = mean(tempCategory2, 4);
                end
            end
            % removes the NaNs from the baseline and response
            category1 = category1(:, iFrequencies, :);
            category2 = category2(:, iFrequencies, :);
            wp = wilcox3d(category1, category2);
            if ~p.Results.squeeze
                % out is time x frequency
                tempWp = NaN(size(wp, 1), obj.nfrequencies);
                tempWp(:, iFrequencies) = wp;
                wp = tempWp;
            end
        end
    end
    
    %% Private helpers
    methods(Access = private)
        function iCategories = categorytonum(obj, categories)
            %Converts categories cells to their numeric counterparts
            if iscell(categories)
                iCategories = obj.PsyData.CategoryNum(categories);
            else
                iCategories = categories;
            end
        end
        
        %% Parsers
        function parser = addchannelsparameter(obj, parser)
            % Adds channels frequencies and categories parameters
            addParameter(parser, 'channels', 1:obj.nchannels,...
                @(x)(isnumeric(x) && (numel(x) > 0)));
        end
        
        function parser = addfrequenciesparameter(obj, parser)
            addParameter(parser, 'frequencies',...
                1:obj.nfrequencies,...
                @(x)(isnumeric(x) && (numel(x) > 0)));
        end
        
        function parser = addcategoriesparameter(obj, parser)
            if isfield(obj.Wp(obj.WpActive), 'kats'), categories = obj.Wp(obj.WpActive).kats;
            else, categories = obj.PsyData.Categories(false); end
            addParameter(parser, 'categories',...
                categories,...
                @(x)(numel(x) > 0));
        end
        
        function parser = addresponseparameter(obj, parser)
             addParameter(parser, 'response',...
                 [obj.baseline(2) obj.epochtime(2)],...
                 @(x)(numel(x) == 2));
        end
        
        function parser = addfdrparameter(~, parser)
             addParameter(parser, 'fdrMethod', [],...
                @(x)isstring(x) && any(strcmp({'pdep' 'dep'},x)));
        end
        %% getters
        function n = nchannels(obj)
         % returns number of channels in the obj.HFreq
            n = size(obj.HFreq, 2);
        end
        
        function n = nfrequencies(obj)
         % returns number of frequencies in the obj.HFreq
            n = size(obj.HFreq, 3);
        end
        
        function n = ncategories(obj)
        % RETURNS total number of categories, not necessarily active onse
            n = size(obj.PsyData.Categories(false), 2);
        end
        
        function n = nsortedchannels(obj)
        % returns number of channels in the obj.CH.sortorder
            n = numel(obj.CH.sortorder);
        end
        
        function i = sortedchannelorder(obj, channel)
        % returns channel order in the obj.CH.sortorder
            i = find(obj.CH.sortorder == channel, 1, 'first');
        end
    end
    
    %% Visualisations
    methods  (Access = private)
        % Buffering of the plot, loading settings from saved ïn the 
        % obj.plotFrequency. Returns existing plot if it exists
        function obj = prepareplotcategoriespowertime(obj, ch, categories)
            if numel(ch) == 0 && ~isfield(obj.plotFrequency, 'ch')
                obj.plotFrequency.ch = obj.CH.sortorder(1);
            else
                obj.plotFrequency.ch = ch;
            end
            
            if numel(obj.plotFrequency.ch) == 1
                obj.plotFrequency.iCurrentChannel = ...
                    obj.sortedchannelorder(obj.plotFrequency.ch);
            end
            
            % Only rewrites categorties if not already set
            if numel(categories) > 0, obj.plotFrequency.kategories = categories;
            elseif ~isfield(obj.plotFrequency, 'kategories') || numel(obj.plotFrequency.kategories) == 0
                if isfield(obj.Wp(obj.WpActive), 'kats'), categories = obj.Wp(obj.WpActive).kats;
                else, categories = obj.PsyData.Categories(); end
                obj.plotFrequency.kategories = categories;
            end
            
            % redraws the plot or creates one
            if isfield(obj.plotFrequency, 'fh') && ishandle(obj.plotFrequency.fh)
                figure(obj.plotFrequency.fh);
            else
                obj.plotFrequency.fh = figure('Name', 'ResponseFreq',...
                    'Position', [20, 500, 400*numel(obj.plotFrequency.kategories), 300]);
                set(obj.plotFrequency.fh, 'KeyPressFcn', @obj.frequencyplothandle);
                colormap jet;
            end
        end
        
        function obj = updateplotcategoriespowertime(obj)
            obj.plotcategoriespowertime(obj.plotFrequency.ch, obj.plotFrequency.kategories);
            obj.plotlabels(obj.plotFrequency.ch, obj.plotFrequency.kategories);
        end
        
        % TODO - deprecate option to pass cells into categories
        % Actual plotting function. Takes channels and categories and
        % creates imagesc based on averaged data (averaged across channels
        % and categories passed)
        function obj = plotcategoriespowertime(obj, channels, categories)
             [miny, maxy] = obj.getplotylimit();
             for k = 1:numel(categories)
                % QUESTION - If I pass it a cell I'd expect the
                % cell to contain NAMES of categories, but that is not what
                % is happening. The categories are not taken as strings,
                % and are not correlated to obj.PsyData.Category name.
                % Neither here, nor in PlotLabels. They are merely taken in
                % their order - this should be DEPRECATED
                if iscell(categories(k))               
                    %cell array is used here for combined categories
                    dd = zeros(size(obj.HFreq, 1), obj.nfrequencies, numel(categories{k})); %time samples x epochs x categories
                    for ikat = 1:numel(categories{k}) 
                        %for one of combined categories
                        dd(:, :, ikat) = obj.getaverageenvelopes('channels', channels, 'categories', categories{k}(ikat));
                    end
                    % QUESTION - WHY IS THIS AVERAGING?
                    D = mean(dd, 3); %mean over combined categories , result is time sample x epochs
                else                    
                    D = obj.getaverageenvelopes('channels', channels, 'categories', categories(k)); %time sample x epochs
                end
                subplot(1, numel(categories), k);
                T = obj.epochtime(1):0.1:obj.epochtime(2);
                imagesc(T, obj.Hfmean, D');
                axis xy; xlabel('time [s]');
                % Buffers the x and y limits to use in future plots
                % THIS can "RESET" the ylimits if the values overshoot the
                % given ones. Kinda "troublesome" in case of scrolling, as
                % each drawn plot redoes the z-axis
                maxy = max([maxy max(D(:))]); miny = min([miny min(D(:))]);
            end
            obj.plotFrequency.ylim = [miny maxy];
        end
        
        % Adds labels to the plotfrequency plot
        % ch: absolute channel numbers
        function obj = plotlabels(obj, channels, categories)
            [miny, maxy] = obj.getplotylimit;
            for k = 1:numel(categories)
                subplot(1, numel(categories), k);
                caxis([miny, maxy]);
                title(obj.PsyData.CategoryName(cellval(categories, k))); %works also for combined categories
                if k == 1
                    channelStrings = iff(numel(channels) > 5, [num2str(numel(channels)) ' channels'],...
                        ['channels ' num2str(channels)]);
                    ylabChannels = iff(isempty(obj.CH.sortedby), channelStrings, ... % if empty obj.CH.sortedby
                        [channelStrings, '(', obj.CH.sortedby, num2str(obj.plotFrequency.ch), ')']);
                    ylabel([ylabChannels ' - freq [Hz]']);
                end
                if k == numel(categories)
                    colorbar('Position', [0.92 0.1 0.02 0.82]);
                end
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
                iCurrentCh = obj.plotFrequency.iCurrentChannel;
                switch eventDat.Key
                   case 'rightarrow'
                       iCurrentCh = min([iCurrentCh + 1, obj.nsortedchannels]);
                   case 'pagedown'
                       iCurrentCh = min([iCurrentCh + 10, obj.nsortedchannels]);
                   case 'leftarrow'
                       iCurrentCh = max([iCurrentCh - 1, 1]);
                   case 'pageup'
                       iCurrentCh = max([iCurrentCh - 10, 1]);
                   case 'space' % zobrazi i prumerne krivky
                       obj.PlotResponseCh(obj.plotFrequency.ch);
                       obj.PlotEpochs(obj.plotRCh.ch, obj.Wp(obj.WpActive).kats);
                       figure(obj.plotFrequency.fh);
                end
                obj.plotFrequency.iCurrentChannel = iCurrentCh;
                obj.plotFrequency.ch = obj.CH.sortorder(iCurrentCh);
                obj.updateplotcategoriespowertime;
            end
            switch eventDat.Key
                case {'multiply', '8'} % dialog na vlozeni minima a maxima osy y
                    answ = inputdlg('Enter ymax and min:', 'Yaxis limits', ...
                        [1 50], {num2str(obj.plotFrequency.ylim)});
                    if numel(answ) == 0, return; end
                    % answ cell 1x1 - cancel is cell 0x0
                    if isempty(answ{1}) || any(answ{1} == '*'), obj.plotFrequency.ylim = [];
                    else
                        data = str2num(answ{:}); %#ok<ST2NM>
                        if numel(data) >= 2, obj.plotFrequency.ylim = [data(1) data(2)]; end
                    end
                    obj.updateplotcategoriespowertime;
                case {'divide', 'slash'} % automaticke meritko na ose z - power
                    obj.plotFrequency.ylim = [];
                    obj.updateplotcategoriespowertime;
                case {'add', 'equal', 's'} % + oznaceni kanalu
                    obj.SelChannel(obj.plotFrequency.ch);
                    obj.updateplotcategoriespowertime;
                % NOT SURE WHAT THIS DOES
                case {'numpad6', 'd'} % + oznaceni kanalu
                    warning('NOT IMPLEMENTED');
                    ch2 = obj.plotRCh.selCh(find(obj.plotRCh.selCh > obj.plotFrequency.ch, 1));
                    ch = iff(isempty(ch2), obj.plotFrequency.ch, ch2);
                    obj.plotFrequency.ch = ch;
                    obj.plotFrequency.iCurrentChannel = obj.sortedchannelorder(ch);
                    obj.updateplotcategoriespowertime;
               % NOT SURE WHAT THIS DOES
               case {'numpad4', 'a'} % + oznaceni kanalu
                    warning('NOT IMPLEMENTED');
                    ch2 = obj.plotRCh.selCh(find(obj.plotRCh.selCh < obj.plotFrequency.ch, 1, 'last'));
                    ch = iff(isempty(ch2), obj.plotFrequency.ch, ch2);
                    obj.plotFrequency.ch = ch;
                    obj.plotFrequency.iCurrentChannel = obj.sortedchannelorder(ch);
                    obj.updateplotcategoriespowertime;
            end
        end       
     end
end