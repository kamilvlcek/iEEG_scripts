function plot_PLV_significance(PLV_data1, PLV_data2, significant, direction)
% visualization of channel pairs with significant PLV difference between 2 periods (or 2 conditions),
% plots all significant chan pairs and saves them to the patient's folder
% PLV_data1 - filename of PLV data which we want to plot
% PLV_data2 - filename of other PLV data (if e.g. we want to plot PLV of delay, bs, and encoding in the same figure)
% significant - if to plot only significant ch pairs (1), or all(0)
% direction - if 1, significant ch pairs where cond 1 > cond 2; if 2 significant ch pairs where cond 2 > cond 1
if(~exist('PLV_data1','var')) || isempty(PLV_data1), PLV_data1 = 'PLV_VTC-IPL_last 2s delay_vs_bs_all_trials_2024-02.mat'; end
if(~exist('significant','var')) || isempty(significant), significant = 1; end % by default, plot only significant ch pairs
if(~exist('direction','var')) || isempty(direction), direction = 1; end % by default, significant ch pairs where cond 1 > cond 2

PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact, some are marked not to analyze as they don't have ROIs we want

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if exist('PLV_data2','var')
            if isfile([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data2])
                load([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data2],...
                    'dataCond2', 'PLVCond2', 'plv_signif_allPairs_clustcorr');
                dataEncod = dataCond2;
                PLVEncod = PLVCond2;
                Encod_plv_signif_allPairs_cc = plv_signif_allPairs_clustcorr;
                chanPairsToPlot2 = 1:size(Encod_plv_signif_allPairs_cc,1);
                clear dataCond2 PLVCond2 plv_signif_allPairs_clustcorr;
            end
        end
        if isfile([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data1])
            % load the PLV data with statistics
            load([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data1]);
            
            if significant == 1    
                if direction == 1 % cond 1 > cond 2
                    % find all rows with at least one positive value
                    positives = plv_signif_allPairs_clustcorr > 0;
                    % Sum across columns to find rows with at least one positive
                    rowsWithPositives = sum(positives, 2) > 0;
                    % rows with negative values
                    negatives = plv_signif_allPairs_clustcorr < 0;
                    % Sum across columns to ensure no negatives in the row
                    rowsWithoutNegatives = sum(negatives, 2) == 0;
                    % Combine conditions: rows with at least one positive and no negatives
                    desiredRows = rowsWithPositives & rowsWithoutNegatives;
                    
                elseif direction == 2 % cond 2 > cond 1
                    % find all rows with at least one negative value
                    negatives = plv_signif_allPairs_clustcorr < 0;
                    % Sum across columns to find rows with at least one negative
                    rowsWithNegatives = sum(negatives, 2) > 0;
                    % rows with positive values
                    positives = plv_signif_allPairs_clustcorr > 0;
                    % Sum across columns to ensure no positives in the row
                    rowsWithoutPositives = sum(positives, 2) == 0;
                    % Combine conditions: rows with at least one positive and no negatives
                    desiredRows = rowsWithNegatives & rowsWithoutPositives;
                end

                % Find indices of these rows
                chanPairsToPlot = find(desiredRows);
            else
                chanPairsToPlot = 1:size(ROI_chanpairs,1);
            end
            
            % plot separate chan pairs with significant difference between 2 periods
            fig_name_part = split(PLV_data1, '.mat');
            fig_filename = [basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' fig_name_part{1}];
            
            % find the max plv value for setting ylim
            maxplv1 = max(max(max(squeeze(PLVCond1.plvspctrm(ROI_chanpairs(chanPairsToPlot,1),ROI_chanpairs(chanPairsToPlot,2), :)))));
            maxplv2 = max(max(max(squeeze(PLVCond2.plvspctrm(ROI_chanpairs(chanPairsToPlot,1),ROI_chanpairs(chanPairsToPlot,2), :)))));
            maxplv = max(maxplv1, maxplv2);
            
            for i = 1:numel(chanPairsToPlot)
                
                ipair = ROI_chanpairs(chanPairsToPlot(i), :);
                significancePos = plv_signif_allPairs_clustcorr(chanPairsToPlot(i),:)>0;
                significanceNeg = plv_signif_allPairs_clustcorr(chanPairsToPlot(i),:)<0;
                
                figure(i);
                h1 = plot(PLVCond1.freq, squeeze(PLVCond1.plvspctrm(ipair(1),ipair(2),:)),'b', 'LineWidth', 1);
                hold on
                h2 = plot(PLVCond2.freq, squeeze(PLVCond2.plvspctrm(ipair(1),ipair(2),:)),'r', 'LineWidth', 1);
                hold on
                ylim([0 maxplv*1.05])
                yLimits = ylim;
                y = yLimits(2) * 0.1;
                
                % plot significance
                ySignifPosValues = ones(1, numel(PLVCond2.freq)) * NaN; % Initialize with NaN
                ySignifPosValues(significancePos) = y; % Assign y values only at significant indices
                plot(PLVCond2.freq, ySignifPosValues, '.-b', 'LineWidth', 3, 'MarkerSize', 10);
                hold on
                
                ySignifNegValues = ones(1, numel(PLVCond2.freq)) * NaN; % Initialize with NaN
                ySignifNegValues(significanceNeg) = y; % Assign y values only at significant indices
                plot(PLVCond2.freq, ySignifNegValues, '.-r', 'LineWidth', 3, 'MarkerSize', 10);
                
                xlim([PLVCond2.freq(1) PLVCond2.freq(end)])
                xlabel('Frequency, Hz')
                ylabel('PLV')
                
                if exist('PLV_data2','var') && ~isempty(PLV_data2) % if we want to plot also encoding
                    significancePos = Encod_plv_signif_allPairs_cc(chanPairsToPlot2(i),:)>0;
                    significanceNeg = Encod_plv_signif_allPairs_cc(chanPairsToPlot2(i),:)<0;
                    y = yLimits(2) * 0.05;
                    hold on
                    h3 = plot(PLVEncod.freq, squeeze(PLVEncod.plvspctrm(ipair(1),ipair(2),:)),'g', 'LineWidth', 1);
                    hold on
                    ySignifPosValues = ones(1, numel(PLVEncod.freq)) * NaN; % Initialize with NaN
                    ySignifPosValues(significancePos) = y; % Assign y values only at significant indices
                    plot(PLVEncod.freq, ySignifPosValues, '.-', 'Color', [0, 0.75, 0.75],'LineWidth', 3, 'MarkerSize', 10);
                    hold on
                    ySignifNegValues = ones(1, numel(PLVEncod.freq)) * NaN; % Initialize with NaN
                    ySignifNegValues(significanceNeg) = y; % Assign y values only at significant indices
                    plot(PLVEncod.freq, ySignifNegValues, '.-g', 'LineWidth', 3, 'MarkerSize', 10);
                    legend([h1, h2, h3],{[dataCond1.condition ': ' num2str(length(dataCond1.trial)) ' trials'], [dataCond2.condition ': ' num2str(length(dataCond2.trial)) ' trials'],...
                        [dataEncod.condition ': ' num2str(length(dataEncod.trial)) ' trials']});
                else
                    legend([h1, h2], {[dataCond1.condition ': ' num2str(length(dataCond1.trial)) ' trials'], [dataCond2.condition ': ' num2str(length(dataCond2.trial)) ' trials']});
                end
                title(['PLV between ' PLVCond1.label{ipair(1), 1} ' ' dataCond1.channelInfo(ipair(1)).ROI ' and ' PLVCond1.label{ipair(2), 1} ' ' dataCond1.channelInfo(ipair(2)).ROI]);
                
                figi_filename = [fig_filename '_chnpair_' PLVCond1.label{ipair(1), 1} ' ' PLVCond1.label{ipair(2), 1} '_' num2str(i) '.fig'];
                savefig(figi_filename);
                close all
            end
            
        end
    end
end
end