function plot_PLV_significance(PLV_data)
% visualization of channel pairs with significant PLV difference between 2 periods, 
% plots all significant chan pairs and saves them to the patient's folder
% PLV_data - filename of PLV data which we want to plot

if(~exist('PLV_data','var')) || isempty(PLV_data), PLV_data = 'PLV_VTC-IPL_last 1.5s delay_vs_bs_all_trials_2024-02.mat'; end
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact, some are marked not to analyze as they don't have ROIs we want

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if isfile([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data])
            % load the PLV data with statistics
            load([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data]);

            % first find indexes of chan pairs with significant PLV difference
            significant_chanPairs = find(sum(plv_signif_allPairs_clustcorr,2) ~= 0);

            % plot separate chan pairs with significant difference between 2 periods
            % (LOC - IPL) ROI1-ROI2
            fig_name_part = split(PLV_data, '.mat');
            fig_filename = [basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' fig_name_part{1}];

            % find the max plv value for setting ylim
            maxplv1 = max(max(max(squeeze(PLVCond1.plvspctrm(ROI_chanpairs(significant_chanPairs,1),ROI_chanpairs(significant_chanPairs,2), :)))));
            maxplv2 = max(max(max(squeeze(PLVCond2.plvspctrm(ROI_chanpairs(significant_chanPairs,1),ROI_chanpairs(significant_chanPairs,2), :)))));
            maxplv = max(maxplv1, maxplv2);
            
            for i = 1:numel(significant_chanPairs)

                ipair = ROI_chanpairs(significant_chanPairs(i), :);
                significancePos = plv_signif_allPairs_clustcorr(significant_chanPairs(i),:)>0;
                significanceNeg = plv_signif_allPairs_clustcorr(significant_chanPairs(i),:)<0;

                figure(i); plot(PLVCond1.freq, squeeze(PLVCond1.plvspctrm(ipair(1),ipair(2),:)),'b', 'LineWidth', 1)
                hold on
                plot(PLVCond2.freq, squeeze(PLVCond2.plvspctrm(ipair(1),ipair(2),:)),'r', 'LineWidth', 1)
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
%                 legend('delay last 1.5 sec', 'baseline 1.5 sec');
                legend(dataCond1.condition, dataCond2.condition);
                title(['PLV between ' PLVCond1.label{ipair(1), 1} ' ' dataCond1.channelInfo(ipair(1)).ROI ' and ' PLVCond1.label{ipair(2), 1} ' ' dataCond1.channelInfo(ipair(2)).ROI]);

                figi_filename = [fig_filename '_chnpair_' PLVCond1.label{ipair(1), 1} ' ' PLVCond1.label{ipair(2), 1} '_' num2str(i) '.fig'];
                savefig(figi_filename);
                close all
            end

        end
    end
end
end