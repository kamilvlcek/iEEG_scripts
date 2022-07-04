classdef CPlots < matlab.mixin.Copyable
    %CPLOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Eh; %handle to CiEEGData object       
        PlotRChMean; %data for  PlotResponseChMean        
        plotTimeInt = struct; %stavove udaje o grafu TimeIntervals
        plotCorrelChanData = struct; % save info about the plot PlotCorelChan %Sofiia 25.11.2020 originally in CiEEGData
    end
    methods (Access = public)
        function obj = CPlots(E) %constructor
            obj.Eh = E; %save handle to main object
        end 
        function PlotResponseChMean(obj,kategories,ylimits, labels, channels,fdrlevel)
            %description ... 
            %arguments: 
            % kategories ? 
            % ylimits ?
            % labels ?
            % channels ?
            % fdrlevel ?
            
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
                if ~isempty(obj.Eh.CH.plotCh2D.chshow)
                    channels = obj.Eh.CH.plotCh2D.chshow;                
                else
                    channels = 1:numel(obj.Eh.CH.H.channels);
                end
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
                [labels,barvy]=obj.Eh.CH.channelPlot.GetBrainlabelsSaved();
                figuretitle = [figuretitle ' - ' obj.Eh.PsyData.CategoryName(cellval(kategories,1))];
            elseif ~strcmp(labels{1},'no') %plotting individual labels
                barvy = distinguishable_colors(max(numel(labels), numel(kategories)));               
                figuretitle = [figuretitle ' - ' obj.Eh.PsyData.CategoryName(cellval(kategories,1))];
            else %plotting categories
                barvy = cell2mat(obj.Eh.colorskat(1:end)'); %from cell array of rgb to matrix
            end            
            colorsErrorBars = obj.GetColorsErrorBars(barvy); 

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
                katnum = cellval(kategories,k);     
                kk = obj.Eh.KatIndex(kategories,k); %index of color
                for l = 1:numel(labels)
                    if strcmp(labels{l},'no') % all brainlabels, plotting categories
                        chsel = channels; %all channels
                        colorkatk = [barvy(kk,:) ; colorsErrorBars(kk,:)]; %dve barvy, na caru a stderr plochu kolem                            
                        legendStr{k} = obj.Eh.PsyData.CategoryName(katnum); % array of names of categories for a legend
                    else % plotting individual brainlabels for one stimulus category
                        brailabels = lower({obj.Eh.CH.brainlabels.label});
                        chsel = intersect(find(contains(brailabels,labels{l})),channels);
                        colorkatk = [barvy(max(k,l),:) ; colorsErrorBars(max(k,l),:)]; %dve barvy, na caru a stderr plochu kolem  
                        legendStr{max(k,l)} = [labels{l} ' x' num2str(numel(chsel)) ':'  obj.Eh.PsyData.CategoryName(katnum)];
                    end
                    [katdata,~,RjEpCh] = obj.Eh.CategoryData(katnum,[],[],chsel);%katdata now time x channels x epochs                                
                    for ich = 1:numel(chsel)
                        CHM(:,k,ich) = mean(katdata(:,chsel(ich),~RjEpCh(ich,:)),3); %mean over valid epochs from katdata
                    end
                    M = nanmean(CHM(:,k,1:numel(chsel)),3);  %mean over channels, nan values can emerge if all epochs for some channel were excluded
                    nanvals = find(isnan(CHM(1,k,1:numel(chsel)))); %count of nan channels (for the first time point)
                    if numel(nanvals) > 0, disp([ legendStr{k} ':' num2str(nanvals') ' are channels with only NaN values']); end                           
                    E = nanstd(CHM(:,k,1:numel(chsel)),[],3)/sqrt(numel(chsel)-numel(nanvals)); %std err of mean over channels          
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
            
            % call separate function for computing statistic
            [adj_p,allTukey] = PlotResponseChMeanStat(obj, CHM, fdrlevel); % return adjusted pvalue after FDR corr
            iadj_p = adj_p <= 0.05;      % Wilcoxon signed rank or ANOVA - find significant time points       
            if ~isempty(allTukey)
                SignPairs = allTukey(iadj_p);  % needs for ploting stars between different cats
            end
            y = ymin*0.95;
            adj_pLims = [find(iadj_p,1) find(iadj_p,1,'last')];
            
            % plot results of Wilcoxon test or ANOVA
            if ~isempty(adj_pLims)
                if numel(kategories)==2 % results of Wilcoxon test for 2 cats
                    y = ymin*0.85;
                    plot(T(iadj_p),ones(1,sum(iadj_p))*y, '*','Color','red'); %stars for Wilcoxon test
                else % general results of ANOVA                  
                    iTimeWndw = obj.Eh.Wp(obj.Eh.WpActive).iepochtime(2,:); %indexes of from what was the statistics computed
                    Tnew = T(iTimeWndw(1):iTimeWndw(2));
                    Tgaps = Tnew; % in case there are non-significant parts in ANOVA results
                    Tgaps(~iadj_p) = NaN; % replace non-signif parts by NaN
                    plot(Tgaps(adj_pLims(1):adj_pLims(2)),ones(1,diff(adj_pLims)+1)*y,'LineWidth',5,'Color',[255 99 71]./255); % line between start and end of significance (could be with gaps)
                    
                    % plot post hoc results only for ANOVA - significance between categories %15.7.2021 Sofiia
                    if ~isempty(SignPairs)
                        yposkat = [3 0 0 0; 4 5 0 0; 6 7 8 0]; %pozice y pro kombinaci k a l - k v radcich a l-1 ve sloupcich
                        ybottom = iff(numel(kategories)>3,0.4,0.3); %odkud se maji umistovat kontrasty mezi kategoriemi - parametr vypoctu y
                        iTimeSign = find(~isnan(Tgaps)); % indexes of significant time points (after ANOVA FDR-corr)                      
                        TukeyResult = zeros(numel(SignPairs), numel(kategories)+1); % initialize matrix for all tukey results

                        for si = 1:numel(SignPairs)  % for each time point where ANOVA result<0.05 (after FDR corr)
                            isignCats = find(SignPairs{si}{:,5} < 0.05);   % find in which pairs pvalue <0.05 in post-hoc tukey table
                            TukeyResult(si,1) = Tgaps(iTimeSign(si)); % save signif time point
                            
                            for ipair = 1:numel(isignCats) % for each sign pair
                                catIndex1 = regexp(SignPairs{si}{isignCats(ipair),1},'(\d+)','match'); % split the letter and index of category, save only index
                                catIndex1 = str2double(catIndex1{1});
                                catIndex2 = regexp(SignPairs{si}{isignCats(ipair),2},'(\d+)','match');
                                catIndex2 = str2double(catIndex2{1});
                                if catIndex1==1
                                    colorstar = barvy(catIndex2,:); %green a red jsou proti kategorii 0, cerna je kat 1 vs kat 2
                                    TukeyResult(si,catIndex2) = 1;   % in aedist: ego > ctrl or allo > ctrl
                                elseif catIndex2==1
                                    colorstar = barvy(catIndex1,:);  
                                else
                                    colorstar = barvy(1,:);  % allo vs ego
                                    TukeyResult(si,4) = 1;
                                end
                                y = ymin + (ymax-ymin)*(ybottom - (yposkat(isignCats(ipair),isignCats(ipair)))*0.05)  ; %pozice na ose y                             
                                plot(Tgaps(iTimeSign(si)),y, '*','Color',colorstar); % plot the star 
                            end                                
                        end                       
                        
                        % create names for pairs of categories - columns names for the table %21.9.2021 Sofiia
                        ColNames = cell(1,numel(legendStr)+1);
                        ni = 0;
                        for k = 1:numel(legendStr)
                            for j = k+1:numel(legendStr)
                                ni=ni+1;
                                ColNames{ni+1} = [legendStr{k} 'X' legendStr{j}];
                            end
                        end
                        ColNames{1} = 'TimePoint'; 
                        % create table for easier manipulations
                        TukeyTable = array2table(TukeyResult,'VariableNames', ColNames); 
                        
                        % names of pairs of categories and tsig on the plot
                        for pi = 1 : numel(kategories)
                            itsig = find(TukeyTable{:,pi+1} == 1, 1,'first'); % first time point of signif
                            y = ymin + (ymax-ymin)*(ybottom - (yposkat(pi,pi))*0.05);
                            if pi == 3
                                colortxt = barvy(1,:);  % allo vs ego
                            else
                                colortxt = barvy(pi+1,:);
                            end
                            text(0.2+obj.Eh.epochtime(1),y, ColNames{pi+1},'FontSize',12, 'Color', colortxt); % name of categories pair
                            if ~isempty(itsig)
                                text(TukeyTable{itsig,1} - 0.07,y, [ num2str(round(TukeyTable{itsig,1}*1000)) 'ms']); % time of start of significance
                            end
                        end
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
        function [adj_p, allTukey] = PlotResponseChMeanStat(obj, CHM, fdrlevel) %Sofiia since 3.6.2021
            % computes statistic for figure PlotResponseChMean
            % if 2 categories, Wilcoxon signed rank is used; otherwise one-way rm ANOVA with FDR correction is used
            % returns adjusted pvalues after FDR corr and cell array with all Tukey post-hoc results for each time point
            % calls from PlotResponseChMean function
            
            % process arguments
            if ~exist('fdrlevel','var') || isempty(fdrlevel)
                fdrlevel = 1;  %less strict fdr correction, 2 is more strict
            end
            
            if size(CHM,2) == 2  % if 2 cats, compute Wilcoxon signed rank
                paired = 1; %paired
                adj_p = CStat.Wilcox2D(CHM(:,1,:), CHM(:,2,:),0,fdrlevel,'kat 1 vs 2',[],[],paired); %more strict dep method
                allTukey = [];
            else
                % time window in which ANOVA should be computed (from stimulus time to the end of epoch)
                iTimeWndw = obj.Eh.Wp(obj.Eh.WpActive).iepochtime(2,:); %indexes of from what was the statistics computed 
                CHMnewTime = CHM(iTimeWndw(1):iTimeWndw(2),:,:);
                % initialize matrix with all p values 
                allpvalue = zeros(size(CHMnewTime,1),1);
                % initialize cell array for all tukey tables - needs for ploting stars
                allTukey = cell(size(CHMnewTime,1),1);
                for ti = 1:size(CHMnewTime,1) % for each time point
                    data1sample =  squeeze(CHMnewTime(ti, :,:)); % cats x chan
                    [pvalue, ~, Tukey_tbl] = CStat.ANOVA1rm(data1sample'); % compute 1-way rm ANOVA for this ti; transpose - chan x cats
                    allpvalue(ti) = pvalue;
                    allTukey{ti} = Tukey_tbl;
                end
                % apply FDR-correction
                if fdrlevel == 2, method = 'dep'; else method = 'pdep'; end 
                [~, ~, adj_p]=fdr_bh(allpvalue,0.05,method,'no'); % 'dep' is more strict                
            end            
        end
        function TimeIntervals(obj,vch,intervals,ylimits, nofile,store) % Sofiia, Kamil since 14.1.2020 
            %  average time intervals for one channel or for a vector of channels (e.g. across a particular structure)
            %  numbers of channels can be obtained, for example, after applying CM.CH.FilterChannels({'lobe','precun'})
            %  or if it's empty it will use numbers of channels after current filtering (from CM.CH.sortoder)
            %  default intervals = [0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8];
            %  nofile - if 1, do not save any xls file, only plot figure
            
            %  TODO move stats to a separate function
                        
            %process arguments
            if exist('vch','var') && ~isempty(vch) && isstruct(vch) && ~exist('store','var')
                store = vch; %ve can use first argument as store struct
                vch = [];
            end
            if ~exist('vch','var') || isempty(vch) %vector of channels
                if ~isempty(obj.Eh.CH.sortorder)
                   vch = obj.Eh.CH.sortorder; % will use numbers of channels after current filtering
                else
                   vch = 1:numel(obj.Eh.CH.H.channels);
                end
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
                kats = obj.Eh.Wp(obj.Eh.WpActive).kats; % define categories
            end
            if ~isfield(obj.plotTimeInt,'outputstyle')
                obj.plotTimeInt.outputstyle = 0; %stores the default value, which can be than changed in var editor
            end
            if ~isfield(obj.plotTimeInt,'interaction')
                obj.plotTimeInt.interaction = 1; % how to perform post-hoc test: 1(default)-category by time, 2-time by category(duration); can be than changed in var editor
            end
            
            % initialize matrix for all channels
            chanMeans = zeros (size(intervals,1), numel(kats), length(vch)); % intervals x kats x channels  -  for statistica - all time intervals and category in columns
            
            %get data from obj.CategoryData for all channels at onces - to make the whole function faster
            categorydata = cell(numel(kats),2); %store precomputed category data here 
            legendStr = cell(1,numel(kats)); % strings for legend and xls table
            for k = 1:numel(kats)                
                [d,~,RjEpCh,~]= obj.Eh.CategoryData(cellval(kats,k),[],[],vch); % d:time x channels x epochs - get data for all channels
                categorydata(k,:) = {d(:,vch,:),RjEpCh}; % save d (time x channels x epochs) and RjEpCh (channels x epochs)
                legendStr{k} = obj.Eh.PsyData.CategoryName(cellval(kats,k)); % array of names of categories for a legend
            end
            
            for ch = 1:size(vch,2) % for each channel                
                % initialize matrix for an individual channel
                allmeans = cell(size(intervals,1),numel(kats)); %  time intervals x kategories (x epochs)
                categories = cell(size(obj.Eh.d,3), 1); % vector to define categories for all epochs for an individual channel - for xls table
                T = (0 : 1/obj.Eh.fs : (size(obj.Eh.d,1)-1)/obj.Eh.fs)' + obj.Eh.epochtime(1); % time in sec
                nInd = 1; % initialize index for variable categories
                
                for k = 1:numel(kats) % for each category%                   
                    RjEpCh = categorydata{k,2}(ch,:); %previously stored RjEpCh
                    d = squeeze(categorydata{k,1}(:,ch,~RjEpCh));  %time x epochs, previously stored d value                  
                    categories(nInd:(nInd+size(d,2)-1)) = repmat({cellval(kats,k)},size(d,2),1); % to establish the number of category for all epochs for an individual channel - for xls table
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
            if exist('store','var') && ~isempty(store) && isstruct(store)  %store or retrieve the to be plotted data  
                [vch,kats,chanMeans,legendStr,chshowstr,ylimits,barvy]=obj.TimeIntervalsStore(store,vch,kats,chanMeans,legendStr,chshowstr,ylimits,intervals);
            else
                store = struct();
            end
            
            
            % plotting
            if isfield(obj.plotTimeInt,'fh') && ishandle(obj.plotTimeInt.fh)
                figure(obj.plotTimeInt.fh); %pouziju uz vytvoreny graf
                clf(obj.plotTimeInt.fh); %graf vycistim
            else
                obj.plotTimeInt.fh = figure('Name','TimeIntervals','CloseRequestFcn', @obj.tearDownFigCallbackTimeIntervals);                 
            end                          
                           
            if length(vch)==1 && isempty(store) && (~isfield(obj.plotTimeInt,'nosinglechannel') || obj.plotTimeInt.nosinglechannel == 0)
                % plot mean and std across epochs for one channel
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
                if obj.plotTimeInt.outputstyle
                    M = M*100; %in percent
                    E = E*100;
                    ylimits = ylimits * 100; %also the ylimits in percents
                end
                if numel(chshowstr)==2
                    if obj.plotTimeInt.outputstyle 
                        title({' '}); %no title but empty lines
                    else
                        title(vertcat( [num2str(length(vch)) ' channels: ' chshowstr{1}] , chshowstr(2) )); %title for multiple channels
                    end
                else
                    title( [num2str(length(vch)) ' channels: ' chshowstr] );
                end
            end            
            
            ploth = zeros(1,numel(kats)); % handles of plots for legend                        
            if ~exist('barvy','var') || isempty(barvy) %if not set inside TimeIntervalsStore  
                if numel(kats)>3                                           
                    barvy = distinguishable_colors(numel(kats));                
                else
                    barvy = cell2mat(obj.Eh.colorskat(1:end)');
                end
            end
            
            if isfield(store,'style') %possibility to set style of the line, including markers
                %as a field inside the store variable. If store.kat and store.store is not set, the store.style can be the only function of store
                if ~iscell(store.style)  %if only string is given
                    store.style = {store.style};
                end %if store.style is cell array, it should be entered with double brackets. e.g. 'style',{{'-d',':d'}}
                if numel(store.style)<numel(kats)
                    style = repmat(store.style,numel(kats),1); %repeat the style for each category
                else
                    style = store.style; 
                end
            else
                style = repmat({'-'},numel(kats),1); %default style with just a line, no markers
            end
               
            for k=1:numel(kats)
                hold on
                kk = obj.Eh.KatIndex(kats,k);
                errorbar(intervals(:, 2),M(:,k),E(:,k),'.','LineWidth',2,'Color', barvy(kk,:)); % plot errors
                hold on;
                h_mean = plot(intervals(:, 2), M(:,k),style{k},'LineWidth',2,'Color',barvy(kk,:),'MarkerSize',10,'MarkerFaceColor', barvy(kk,:)); % plot means %'none'
                ploth(k) = h_mean; % plot handle                
            end
                      
            if numel(kats)==2 && (~isfield(store,'stat') || store.stat ~= 0) %for just two categoires, plot nonparametric paired fdr corrected statistics
                %if store.stat ==0 , no stars for statistics are ploted
                paired = 1; %paired
                fdrlevel= 1; %actually ANOVA would be better, with more standard results
                if ~isempty(ylimits) 
                    y = ylimits(1)*0.8;
                else
                    y = 0;
                end
                Wp = CStat.Wilcox2D(chanMeans(:,1,:), chanMeans(:,2,:),0,fdrlevel,'kat 1 vs 2',[],[],paired); %more strict dep method
                iWp = Wp <= 0.05;
                plot(intervals(iWp, 2),ones(1,sum(iWp))*y, '*','Color','red'); %stars
                
            elseif numel(kats) > 2 && (~isfield(store,'stat') || store.stat ~= 0) && length(vch)>1  % compute 2 rm ANOVA for set of channels
                %Sofiia since 12.3.2021 
                [tableX,~] = obj.TimeIntervals2XLSTable(vch,chanMeans,kats,legendStr,intervals); % get data in prepared table
                data = cell2mat(tableX(:, 7:end)); % only numerical data
                levelNames = cell(2,1);
                
                % convert to strings - double quotes - needs for function CStat.ANOVA2rm
                tempStr = string(legendStr); % names of categories
                for stri = 1:numel(tempStr)
                    levelNames{1}{stri} = tempStr(stri);
                end
                strInterv = string(num2str(intervals,'%.1f-%.1f')); % names of intervals
                for stri = 1:numel(strInterv)
                    levelNames{2}{stri} = strInterv(stri);
                end                
                factorNames = {'Category', 'Time'};
                interact = obj.plotTimeInt.interaction; % 1 - category by time; 2 - time by category (to find duration of response); %Sofiia since 20.5.2021 
                [ranova_table, Tukey_table] = CStat.ANOVA2rm(data, factorNames, levelNames, interact,1); % compute 2-way repeated measures ANOVA
                
                % define y limits for ploting stars
                if ~isempty(ylimits)
                    yrange = ylimits;
                else
                    yrange = get(gca, 'ylim');
                end
                ydiff = 0;
                ys = yrange(1)+diff(yrange)*0.02; % initial y coordinate of star
                % add pvalue and F statistic on the figure
                text(intervals(1,2),yrange(2)*0.96,['F(' num2str(ranova_table{7,2}) ', ' num2str(ranova_table{8,2}) ') = ' num2str(ranova_table{7,4}) ', p = ' num2str(ranova_table{7,5})])
                
                % find significant pairs for ploting stars              
                if ~isempty(Tukey_table)  %Sofiia since 20.5.2021  
                    % find significant pairs
                    [~, ipval, ~] = unique(Tukey_table{:,6},'stable'); % only unique pvalue for one pair
                    Tukey_tableUniq = Tukey_table(ipval,:);
                    iSignPairs = Tukey_tableUniq{:,6} <= 0.05; % pvalue <= 0.05
                    signPairs = table2cell(Tukey_tableUniq(iSignPairs, 1:4)); % signif time intervals, category1, category2 and their difference without repetition
                    
                    % plot stars of significance between kats on the same time interval (now only for 3 categories) 
                    % TODO - for any number of categories
                    if interact == 1 && numel(kats) == 3
                        for ipair=1:size(signPairs,1)
                            color = [];
                            if signPairs{ipair,4}<0 && signPairs{ipair,2} == legendStr{1} && signPairs{ipair,3} == legendStr{2}
                                color = barvy(2,:); % asterisk of 2 category's color - cat2 > cat1
                            elseif signPairs{ipair,4}<0 && signPairs{ipair,2} == legendStr{1} && signPairs{ipair,3} == legendStr{3}
                                color = barvy(3,:); % asterisk of 3 category's color - - cat3 > cat1
                            elseif (signPairs{ipair,2} == legendStr{2} && signPairs{ipair,3} == legendStr{3}) || (signPairs{ipair,2} == legendStr{3} && signPairs{ipair,3} == legendStr{2})
                                color = [1 0 1]; % asterisk of magenta color - - cat3 > cat2 or cat2 > cat3
                            end
                            
                            iInterv = strcmp(signPairs{ipair,1},strInterv); % define time interval
                            
                            if ~isempty(color) && ipair>1 && signPairs{ipair,1}==signPairs{ipair-1,1} % to plot several stars for one time interv
                                ydiff = ydiff + diff(yrange)*0.02;
                                plot(intervals(iInterv,2), ys+ydiff, '*','Color',color); % plot star
                            elseif ~isempty(color)
                                plot(intervals(iInterv,2), ys, '*','Color',color);
                                ydiff = 0;
                            end
                        end
                        
                        % legend for stars
                        text(intervals(1,2),yrange(2)*0.9,['* ' legendStr{3} ' vs ' legendStr{2}],'Color',[1 0 1])
                        text(intervals(1,2),yrange(2)*0.86,['* ' legendStr{2} ' > ' legendStr{1}],'Color',barvy(2,:))
                        text(intervals(1,2),yrange(2)*0.82,['* ' legendStr{3} ' > ' legendStr{1}],'Color',barvy(3,:))
                        
                    elseif interact == 2 % plot stars of significance between 1st time int and all others for each category separately (duration)
                        %Sofiia 20.5.2021 
                        % leave only sign pairs with the first time interval
                        iSignFirst = strcmp(string(signPairs(:,2)),strInterv(1));
                        SignFirst = signPairs(iSignFirst,:);
                        
                        % plot stars for each kat
                        for k = 1:numel(legendStr)
                            ikat = strcmp(string(SignFirst(:,1)),tempStr(k)); % find all sign time intervals for this kat
                            ikat = find(ikat);
                            for iOnekat = 1:numel(ikat) % for each time int
                                iInterv = strcmp(string(SignFirst(ikat(iOnekat),3)),strInterv); % define time interval
                                plot(intervals(iInterv,2), ys+ydiff, '*','Color',barvy(k,:)); % plot star on this time int
                            end
                            % legend for star
                            text(intervals(1,2)-0.03, yrange(2)*(0.9 - ydiff*3),['* ' legendStr{k} '(t) > ' legendStr{k} ' (' char(strInterv(1)) ')'],'Color',barvy(k,:))
                            
                            ydiff = ydiff + diff(yrange)*0.02; % new y of next star
                        end
                    end
                end
            end
            
            if ~isempty(ylimits), ylim(ylimits); end
            xlim([intervals(1,2)-.05, max(max(intervals))+.05]);
            xticks([intervals(1, 1) (intervals( : , 2))'])
            if obj.plotTimeInt.outputstyle %for output style very simple x tickmarks
                txtinterv = num2str(intervals(:,2),'%0.1f'); %%.1f-%.1f %only higher parts of intervals
                txtinterv = char(strrep(cellstr(txtinterv),'0.','.')); %delete zeros before dot
                %txtinterv = horzcat(repmat(' ',size(txtinterv,1),15),txtinterv); %move ticklabels a bit right by preceding spaces
                yticks(0 : 20 : ylimits(2)); %ticks only 
            else
                txtinterv = num2str(intervals,'%.1f-%.1f'); 
            end
            xticklabels({'0', txtinterv(1:end, :)})            
            
            if ~obj.plotTimeInt.outputstyle %for output style no legend and no xlabel
                xlabel('time intervals, s');
                legend(ploth,legendStr,'Location','best'); 
            else
                set(obj.plotTimeInt.fh,'Name',['TimeIntevals - ',cell2str(legendStr)]) %change the figure window title 
                set(gca,'linewidth',2)
                set(gca,'FontSize',24)
            end
            
            if ~isfield(obj.plotTimeInt,'filterListener') || isempty(obj.plotTimeInt.filterListener) %function to be calle when event CHHeader.FilterChanged happens
                obj.plotTimeInt.filterListener = addlistener(obj.Eh.CH, 'FilterChanged', @obj.filterChangedCallbackTimeIntervals);
            end
                        
            % export data in xls table
            if length(vch)>1 && numel(kats)>2 && obj.plotTimeInt.nofile==0
                obj.TimeIntervals2XLS(vch,kats, legendStr, intervals, chanMeans,chshowstr, categories, allmeans,ranova_table, Tukey_table); %Sofiia 12.3.2021 
            elseif obj.plotTimeInt.nofile==0 && isempty(store) 
                obj.TimeIntervals2XLS(vch,kats, legendStr, intervals, chanMeans,chshowstr, categories, allmeans);
            end
        end
        function TimeIntervalsColect(obj,marks,labels,labelsAreCluster)
            %collects data using TimeIntervals from multiple filter settings
            %marks is cell array with two columns: 1. no of marking sets 2. markings one of fghjkl. e.g. marks = {1,1,1,2,2,2;'','j','k','j','k','l'}';
            %labels is cell array of brain labels, filtered by label = '' ; if empty, uses stored labels   in channelPlot.plotCh3D.brainlabelColors         
            %the intervals and ylimis needs to be set before
            if ~exist('labels','var') || isempty(labels) 
                labels = {obj.Eh.CH.channelPlot.plotCh3D.brainlabelColors.label};
            elseif ischar(labels) && strcmp(labels,'clear')
                obj.plotTimeInt.data = struct; %clear all stored data
                disp('all TimeIntervals data cleared');
                return;
            end
            if ~exist('labelsAreCluster','var')
                labelsAreCluster = 0; %if labelsAreCluster==1, the labels var contains names of clusters
            end
            if exist('marks','var')
                assert(iscell(marks),'marks has to be a cell array');                
                if size(marks,1) == 2 && size(marks,2)>2
                    marks = marks'; %if the dimensions are swaped, make them correct
                end
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
                        labelflag = iff(labelsAreCluster,'cluster','label'); %filter labels or clusters? even clusters names are stored under labels column in data
                        obj.Eh.CH.FilterChannels({labelflag,labels{l}},{},[marks{m,2} 'n']);
                        obj.TimeIntervals([],[],[],1,struct('store',1,'label',labels{l},'mark',markname));
                    end
                end            
            end
        end
        function TimeIntervals2XLS(obj,varargin)
            if nargin > 2 %when called directly from TimeIntervals()
               %vch,kats , legendStr, intervals, chanMeans,chshowstr,categories, allmeans, (ranova_table, Tukey_table - Sofiia 12.3.2021 )
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
            if nargin>9 % Sofiia 12.3.2021 - additional sheets with ANOVA statistic in xls file
                ranova_table = varargin{9};
                Tukey_table = varargin{10};
                writetable(ranova_table,xlsfilename,'Sheet','rANOVA_result','WriteRowNames',1); % export ANOVA statistic
                if ~isempty(Tukey_table)
                    writetable(Tukey_table,xlsfilename,'Sheet','Tukey_result');
                end
            end
            disp([ 'xls tables saved: ' xlsfilename]);
            
        end
        function obj = GetCorrelChan(obj,nofile,fraction,pval_limit,rho_limit) % Sofiia since 11.2020
            % find channels in which neural response is correlated with response time of patient (to be excluded from the further analysis)
            % correlation between t of max response and patient's RT in all channels are saved in CM.CH.correlChan (struct array)
            % nofile - if 0, saves the xls file with correlation coeff, pvalue, name of channel, patient
            
            assert(obj.Eh.epochs > 1,'only for epoched data');
            if ~exist('nofile','var'), nofile = 0; end % save the xls table by default
            if ~exist('fraction','var'), fraction = 0.9; end % time of 90% of maximum
            if ~exist('pval_limit','var'), pval_limit = 0.01; end % threshold for p value
            if ~exist('rho_limit','var'), rho_limit = 0.2; end % threshold for correlation coeff
            
            % initialize cell aray: n channel,patient,name of channel,neurology label, r coeffic, pvalue, correl (0/1)
            output = cell(size(obj.Eh.d,2),7);
            
            chanSwitch = 1;
            for iPacient = 1:obj.Eh.PsyData.nS % over each patient
                obj.Eh.PsyData.SubjectChange(iPacient); % change current subject in CPsyDataMulti
                [resp,RT,~,test] = obj.Eh.PsyData.GetResponses(); % get his RT
                correctRT = RT(resp==1 & test==1); % get RT only for correct test trials
                
                for iChannel = chanSwitch : obj.Eh.els(iPacient)   % for channels of one patient
                    
                    [rho,pval] = obj.ComputeCorrelChan(iChannel,fraction,correctRT,resp==1 & test==1);  %  
                    
                    if pval < pval_limit && rho >= rho_limit, correlCh = 1; else, correlCh = 0; end % weak-to-moderate correlation
                    
                    % put everything in cell array for the futher export
                    output{iChannel,1} = iChannel;  % n of channel
                    output{iChannel,2} = obj.Eh.PsyData.P.pacientid; % name of patient
                    output{iChannel,3} = obj.Eh.CH.H.channels(iChannel).name; % name of channel
                    output{iChannel,4} = obj.Eh.CH.H.channels(iChannel).neurologyLabel; % neurology label of channel
                    output{iChannel,5} = rho;
                    output{iChannel,6} = pval;
                    output{iChannel,7} = correlCh; % correlated activity with RT
                                                            
                end
                chanSwitch = obj.Eh.els(iPacient)+1; % switch to channels of the next patient
            end
            
            % save channels correlated/uncorrelated (1/0) with RT to CHHeader
            obj.Eh.CH.StoreCorrelChan(cell2mat(output(:,7))); %store the correlCh values 1/0 in the the CH.CorrelChan field
            
            if nofile == 0  % export xls table
                variableNames = {'channel' 'pacient' 'channelName' 'neurologyLabel' 'rho' 'p-value' 'correlated'...
                    'fraction' 'pval_limit' 'rho_limit'};
                output(1,8:10) = {fraction, pval_limit, rho_limit}; % save also parameters of function
                assert(~isempty(obj.Eh.hfilename),'the object is not saved, please save the object first');
                [~,mfilename,~] = fileparts(obj.Eh.hfilename); %hfilename se naplni az pri ulozeni objektu
                mfilename = strrep(mfilename, ' ', '_');
                logfilename = ['CorrelChan_' mfilename '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') ];
                xlsfilename = fullfile('logs', [logfilename '.xls']);
                xlswrite(xlsfilename ,vertcat(variableNames,output)); % write to xls file
                disp(['xls table ' xlsfilename ' was exported']);
            end
        end        
        function obj = PlotCorrelChan(obj,ch,tmax,fraction,conditions, nofile) % Sofiia since 11.2020
            % creates scatterplot between t of max neural response (or max size of response) and
            % patient's RT for one specific channel either across all trials together or for each category independently
            % also can export xls table with all epochs data of each condition
            
            assert(obj.Eh.epochs > 1,'only for epoched data');
            if ~exist('ch','var') || isempty(ch), ch = obj.Eh.CH.sortorder(obj.Eh.plotRCh.ch); end % use current channel from PlotResponseCh plot
               % ch = obj.Eh.CH.sortorder(1)
            if ~exist('fraction','var') || isempty(fraction), fraction = 0.9; end % time of 90% of maximum
            if ~exist('tmax','var') || isempty(tmax), tmax = 1; end % by default computes time of maximum, if 0 - computes valmax (max size of response)
            if ~exist('conditions','var') || isempty(conditions), conditions = 0; end % by default computes for all epochs together, without division to conditions
            if ~exist('nofile','var') || isempty(nofile), nofile = 1; end % by default doesn't export xls table with all epochs of each condition (beh RT and max eeg value)
            
            % return the name of the patient of this channel
            % patId = obj.Eh.CH.PacientTag(ch);
            % find index of that patient in PsyData
            obj.Eh.PsyData.SubjectChange(find(obj.Eh.els >= ch,1)); %change current subjet in CPsyDataMulti (but if called from PlotResponseCh it was changed alread) or do nothing if the file contains only one patient
            
            if conditions == 1 % if we want to get data for each condition separately
                kats = obj.Eh.Wp(obj.Eh.WpActive).kats; % define categories
                maxTrialsKats = cell(numel(kats),1); %  maxtrials and RT  x categories
                corrKats = zeros(numel(kats), 2); % to store rho and pval for each category (for futher plotting)
                legendStr = cell(1,numel(kats)); % strings for legend and xls table
                
                for k=1:numel(kats)
                    katnum = cellval(kats,k);
                    [katdata,psy_rt, RjEpCh] = obj.Eh.CategoryData(katnum,[],[],ch); % katdata - time points x channel x trials
                    katdata = squeeze(katdata(:,ch,~RjEpCh));   % exclude incorrect and bad epochs, time points x trials
                    correctRT = psy_rt(~RjEpCh); % also exclude them from RT
                    [rho,pval,maxTrials] = obj.ComputeCorrelChan(ch,tmax,fraction,correctRT,[], conditions, katdata);  % get corr and max value for all epochs of this kat
                    maxTrialsKats{k} = [maxTrials correctRT]; % to put data of each kat in one cell
                    corrKats(k,:) = [rho pval];
                    legendStr{k} = obj.Eh.PsyData.CategoryName(katnum); % array of names of categories for a legend
                end
                
            else  % if we want to get data of all epochs together
                % get patient's RT for all trials
                [resp,RT,~,test] = obj.Eh.PsyData.GetResponses();    %its better to use existing functions where available
                correctRT = RT(resp==1 & test==1); % get RT only for correct test trials - in CHilbertMulti files the traning is excluded, but in CHilbert files for one patient only not.
                [rho,pval,maxTrials] = obj.ComputeCorrelChan(ch,tmax,fraction,correctRT,resp==1 & test==1,conditions);  % get corr and max value for all correct epochs
                
                % save values in the structure obj.plotCorrelChanData
                obj.plotCorrelChanData.fraction = fraction;
                obj.plotCorrelChanData.pval = pval;
                obj.plotCorrelChanData.rho = rho;
                obj.plotCorrelChanData.chan = ch;
            end
            
            % ploting
            if isfield(obj.plotCorrelChanData,'h') && ~isempty(obj.plotCorrelChanData.h) && ishandle(obj.plotCorrelChanData.h)
                figure(obj.plotCorrelChanData.h) % use the existing plot
            else
                obj.plotCorrelChanData.h = figure('Name','PlotCorrelChan scatterplot'); % or new plot
            end
            fig = gcf; % current figure handle
            
            if conditions == 1 % if we want to plot data for each condition separately
                ploth = zeros(1,numel(kats)); % handles of subplots
                if numel(kats)>3 % to get colours for all categories
                    colors = distinguishable_colors(numel(kats));
                else
                    colors = cell2mat(obj.Eh.colorskat(1:end)');
                end
                for k=1:numel(kats) % for each category
                    kk = obj.Eh.KatIndex(kats,k);
                    ploth(k) = subplot(1,numel(kats),k);
                    scatter(maxTrialsKats{k}(:,1),maxTrialsKats{k}(:,2), 25, colors(kk,:),'filled'); % chart
                    hold on %Sofiia 31.3.2022
                    plot([(min(maxTrials)-0.1) (max(maxTrials)+0.1)], [median((maxTrialsKats{k}(:,2))) median((maxTrialsKats{k}(:,2)))], ':', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.5); % plot median behav RT of each category
                    title([legendStr{k} ' (' num2str(length(maxTrialsKats{k}(:,1))) '), rho = ' num2str(corrKats(k,1),2) ', p value = ' num2str(corrKats(k,2),3)])
                end
                linkaxes(ploth, 'xy') % link axes of all subplots to specify their limits later
                fig.Position = [20, 300, 1800, 500];
            else
                scatter(maxTrials, correctRT, 25,'m','filled') % to plot data of all epochs
                title(['All epochs (' num2str(length(maxTrials)) '), rho = ' num2str(rho,2) ', p value = ' num2str(pval,3)])
            end
            
            if tmax == 1
                xname = ['t of ' num2str(fraction) ' max response [s], fraction = ' num2str(fraction)];
            else
                xname = 'max BGA response';
            end
            
            ax = findobj(fig,'Type','Axes');
            ylabel(ax(length(ax)),'patients RT [s]') % ylabel only for the first subplot
            xlabel(ax(length(ax)),xname)
            set(ax,'FontSize',12)
            set(gca, 'Xlim', [(min(maxTrials)-0.1) (max(maxTrials)+0.1)], 'Ylim', [0 (max(correctRT)+0.1)]);
            xrange = get(gca,'XLim');
            yrange = get(gca,'YLim');
            text(xrange(1)+0.2,yrange(1)+0.1,['channel ' num2str(ch) ', ' obj.Eh.CH.H.channels(ch).name ', ' obj.Eh.CH.H.channels(ch).neurologyLabel], 'FontSize', 13)
            
            % export epochs data of each condition to xls table 
            if nofile == 0 && conditions == 1                
                allData = cell2mat(maxTrialsKats); % join all cats (max value and RT) in one array
                conditionCol = strings(size(allData,1),1); % string array for categories column
                kInd = 1; % index for each cat
                
                for k=1:numel(kats) % for each category
                    conditionCol(kInd :(kInd-1+length(maxTrialsKats{k}(:,1)))) = string(repmat(legendStr{k}, length(maxTrialsKats{k}(:,1)), 1)); % repeat category name for all corresponding epochs
                    kInd = kInd + length(maxTrialsKats{k}(:,1));
                end
                
                % put all data in cell array
%                 output = num2cell([allData, conditionCol]);
%                 output{1,4} = ch;  % n of channel
%                 output{1,5} = patId; % name of patient
%                 output{1,6} = obj.Eh.CH.H.channels(ch).name; % name of channel
%                 output{1,7} = obj.Eh.CH.H.channels(ch).neurologyLabel; % neurology label of channel                
%                 columnNames = {'max_eeg' 'patients_RT' 'condition' 'channel' 'patient' 'channelName' 'NeurologyLabel'};
%                 
                outputTable = table(allData(:,1), allData(:,2), conditionCol,'VariableNames',{'max_eeg' 'patients_RT' 'condition'});
                logfilename = ['Epoch data for channel ' num2str(ch) '_' obj.Eh.CH.H.channels(ch).name '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') ];
                xlsfilename = fullfile('logs', [logfilename '.xls']);
                %xlswrite(xlsfilename, outputTable); % write to xls file
                writetable(outputTable, xlsfilename)
                disp(['xls table ' xlsfilename ' was exported']);
            end
        end
        
        function [rho,pval,maxTrials,obj] = ComputeCorrelChan(obj,ch,tmax,fraction,correctRT,correctTrials, conditions, katdata) %Sofiia since 11.2020             
            %function to compute correlation between behavioural RT and time of eeg resp maximum (or valmax) for one channel
            % to remove duplicite code in PlotCorrelChan and GetCorrelChan
            %arguments: ch - index in obj.Eh.d, tmax - 1=returns time of maximum,0=return valmax of maximum
            % fraction - fraction of maximum, only used for tmax = 1
            % correctRT - array of behavioural reaction times to correlate - for individual epochs
            % correctTrials - ?
            % conditions - ?
            % katdata - ?
            %output arguments [rho ? ,pval ?,maxTrials ?]            
            
            if ~exist('tmax','var') || isempty(tmax), tmax = 1; end % by default computes time of maximum, if 0 - computes valmax
            % initialize matrix for all trials
            maxTrials = zeros (length(correctRT),1); %time of maximum eeg response or valmax
            assert(isfield(obj.Eh.Wp,'iepochtime'),'the statistics needs to be computed');
            iintervalyData = obj.Eh.Wp(obj.Eh.WpActive).iepochtime(2,:); %indexes of from what was the statistics computed , saved in function ResponseSearch
            
            if conditions == 1
                for iTrial = 1:size(katdata,2) % over trials of one category
                    d = katdata(iintervalyData(1):iintervalyData(2), iTrial); % time points x channel x trials
                    [valmax, idxMax, idxFrac] = cMax(d, fraction);
                    if tmax == 1 % if we want to compute corr between behavioural RT and time of maximum eeg resp
                        if fraction >= 1 %if we are interested in the real maximum
                            maxTrials(iTrial) = idxMax/obj.Eh.fs;  % put it in matrix for all trials, in seconds
                        else %if we are interested in a fraction of maximum (i.e. fraction < 1)
                            maxTrials(iTrial) = idxFrac/obj.Eh.fs;  % put it in matrix for all trials, in seconds
                        end
                    else
                        maxTrials(iTrial) = valmax; % if we want to compute corr between behavioural RT and size of eeg resp
                    end
                end
            else
                imaxTrials = 1; %index in maxTrials
                for iTrial = 1:size(obj.Eh.d,3)  % over all trials without division to categories
                    if correctTrials(iTrial)
                        d = squeeze(obj.Eh.d(iintervalyData(1):iintervalyData(2), ch, iTrial)); % obj.d - time points x channels x trials (128 x 424 x 384) - without the time before stimulus
                        [valmax, idxMax, idxFrac] = cMax(d, fraction);
                        if tmax == 1 % if we want to compute corr between behavioural RT and time of maximum eeg resp
                            if fraction >= 1 %if we are interested in the real maximum
                                maxTrials(imaxTrials) = idxMax/obj.Eh.fs;  % put it in matrix for all trials, in seconds
                            else %if we are interested in a fraction of maximum (i.e. fraction < 1)
                                maxTrials(imaxTrials) = idxFrac/obj.Eh.fs;  % put it in matrix for all trials, in seconds
                            end
                        else
                            maxTrials(imaxTrials) = valmax; % if we want to compute corr between behavioural RT and size of eeg resp
                        end
                        imaxTrials = imaxTrials +1;
                    end
                end
            end
            % compute correlation between max Power and RT
            [rho,pval] = corr(maxTrials, correctRT,'Type', 'Spearman');
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
            tableX = [num2cell(vch), names, pacients, labels, num2cell(chanMeans)]; % final table
            
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
        function [vch,kats,chanMeans,legendStr,chshowstr,ylimits,barvy]=TimeIntervalsStore(obj,store,vch,kats,chanMeans,legendStr,chshowstr,ylimits,intervals)    
            %stores or retrieves the to be plotted data in obj.plotTimeInt.data
            assert(numel(store) == 1,'store can be only of size 1');
            if isfield(store,'store') && store.store==1
                %store the values into obj.plotTimeInt.data
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
                    if isempty(n), n = numel(obj.plotTimeInt.data)+1; end %new combination of label and mark
                    obj.plotTimeInt.data(n).label = store.label;
                    obj.plotTimeInt.data(n).mark = store.mark;
                    obj.plotTimeInt.data(n).chanMeans = chanMeans; % intervals x kats x channels
                    obj.plotTimeInt.data(n).chshowstr = chshowstr; %current filter label
                    obj.plotTimeInt.data(n).legendStr = legendStr; %current filter label
                    obj.plotTimeInt.data(n).chnum = numel(vch); %number of channels
                    obj.plotTimeInt.data(n).channels = vch; %channel numbers
                    obj.plotTimeInt.data(n).intervals = intervals; %number of channels
                    obj.plotTimeInt.data(n).kats = kats; %number of categories
                    disp(['data stored as "' store.label ' & ' store.mark '"']);
                    barvy = []; %only to return something
            elseif isfield(store,'kat') 
                assert(min(store.kat) > 0 && max(store.kat) <= size(obj.plotTimeInt.data(1).chanMeans,2),['wrong kat numbers: ' num2str(store.kat) , ', should be 1 to ' num2str(size(obj.plotTimeInt.data(1).chanMeans,2))]);
                %number of category to plot
                if isfield(store,'mark')
                    marks = {obj.plotTimeInt.data.mark};
                    idata_m = contains(marks,store.mark)'; %logical index  n x 1                   
                else
                    store.mark = 'NaN';
                    idata_m = true(numel(obj.plotTimeInt.data),1);                    
                end
                if isfield(store,'label')
                    labels = {obj.plotTimeInt.data.label};
                    idata_l = contains(labels,store.label)'; %logical index n x 1                
                else
                    store.label = 'NaN';
                    idata_l = true(numel(obj.plotTimeInt.data),1);                    
                end
                idata = all([idata_m, idata_l],2); %n x 1, logical index of all rows of data to retrieve
                assert(sum(idata)>0, 'no data to plot');
                idata1 = find(idata,1); %index of first row                
                datan = find(idata); %rows in plotTimeInt.data
                if isfield(store,'markstogether') && store.markstogether==1 
                    L = unique({obj.plotTimeInt.data(idata).label}); %unique labels for this selection
                    if numel(datan)==numel(L)*numel(store.mark)
                        datan = reshape(datan,numel(store.mark),[])'; % make size numel(L) x numel(store.mark) 
                    else
                        datan = datan'; %treat all rows of obj.plotTimeInt.data together
                    end
                end
                
                if ~isfield(store,'normalize') %normalize each label interval values to it maximal interval mean - to be able to compare the speed of response
                    store.normalize = 0;
                elseif store.normalize == 1
                    ylimits = [-.1 1.1]; %the values will be [0 1], temporary change 
                end
                
                if size(datan,1) == 1 %if only one row of obj.plotTimeInt.data should be plotted
                    kats = obj.plotTimeInt.data(idata1).kats(store.kat); %real number of categories as stored 
                else
                    kats = 1 : numel(store.kat) * size(datan,1); %mimic different stored (labels and kats) as categories
                end
                
                if isfield(store,'markstogether') && store.markstogether==1 
                    chanMeans = nan(size(obj.plotTimeInt.data(idata1).chanMeans,1), numel(kats) , sum([obj.plotTimeInt.data(idata).chnum])); %intervals x (labels+kats) x channels
                else
                    chanMeans = nan(size(obj.plotTimeInt.data(idata1).chanMeans,1), numel(kats) , max([obj.plotTimeInt.data(idata).chnum])); %intervals x (labels+kats) x channels
                end
                legendStr = cell(1,numel(kats)); 
                
                if isfield(store,'colors') && strcmp(store.colors, 'nrj') %if to use colors from CM.CH.channelPlot.plotCh3D.brainlabelColors 
                    [brainlabels,labels_barvy] = obj.Eh.CH.channelPlot.GetBrainlabelsSaved();
                    barvy = zeros(numel(kats),3);
                elseif isfield(obj.plotTimeInt,'colorskat')
                    barvy = obj.plotTimeInt.colorskat; %preferentially load previously stored and maybe modified colors
                elseif size(datan,1)==1
                    barvy = cell2mat(obj.Eh.colorskat(2:end)');
                    obj.plotTimeInt.colorskat = barvy; %store colors so that they can be modified in var editor
                else
                    barvy = distinguishable_colors(numel(kats));
                end
                
                vch = []; %list of all displayed channels
                for d = 1:size(datan,1) %rows of plotTimeInt.data 
                     for ikat = 1:numel(store.kat) %all kats are in single row of obj.plotTimeInt.data, this cycle is using this single row
                        k = (d-1)*numel(store.kat)+ikat; %index in chanMeans, 1:numel(kats)
%                         if numel(datan(d,:))>1
                            ichanMeans = 1;
                            for j = 1:numel(datan(d,:))
                                chM =  obj.plotTimeInt.data(datan(d,j)).chanMeans(:,store.kat(ikat),:);
                                chanMeans(:,k,ichanMeans:ichanMeans+size(chM,3)-1) = chM;
                                ichanMeans = ichanMeans + size(chM,3);
                                vch = union(vch,obj.plotTimeInt.data(datan(d,j)).channels);
                            end
                            ichanMeans = ichanMeans -1; %make it number of channels 
%                         else
%                             chanMeans(:,k,1:obj.plotTimeInt.data(datan(d)).chnum) = obj.plotTimeInt.data(datan(d)).chanMeans(:,store.kat(ikat),:);
%                             vch = union(vch,obj.plotTimeInt.data(datan(d,:)).channels);
%                         end
                        
                        if store.normalize
                            katmax = max(mean(chanMeans(:,k,1:ichanMeans),3));
                            chanMeans(:,k,1:ichanMeans) = chanMeans(:,k,1:ichanMeans)./katmax;
                        end
                        legendStr{k} = [ obj.plotTimeInt.data(datan(d,1)).label ',' obj.plotTimeInt.data(idata1).legendStr{ikat} ' x' num2str(ichanMeans)];
                        if isfield(store,'colors') && strcmp(store.colors, 'nrj')
                            ilabel = contains(brainlabels,obj.plotTimeInt.data(datan(d,1)).label);
                            barvy(d,:) = labels_barvy(ilabel,:);                        
                        end
                    end
                end
                chshowstr = {cell2str(obj.plotTimeInt.data(idata1).legendStr(store.kat)) ['mark=' cell2str(store.mark) ] }; %assumes same category accross data fields                   
            end            
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
        function colorsErrorBars = GetColorsErrorBars(barvy)
            HSV = rgb2hsv(barvy);
            HSV(:,2) = 0.15; %decrease saturation for error bands
            colorsErrorBars = hsv2rgb(HSV); % barvycellfun(@(a) min(a+hue, 1), barvy, 'UniformOutput', false);
        end
    end
    
end

