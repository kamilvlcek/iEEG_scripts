%% set up a path for data
PLV_data =  'PLV_VTC-IPL_last 2s delay_vs_bs_all_trials_2024-02.mat';
% PLV_data =  'PLV_VTC-IPL_last 2s delay_vs_2s encoding_all_trials_2024-02.mat';
PLV_folder = 'PLV_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact();

%%
% initialize 2 cells, one for PLV of delay, and one for PLV of bs
PLVdelay_allSubj = [];
PLVdelay_allSubj.plvspctrm = {};

PLVbs_allSubj = [];
PLVbs_allSubj.plvspctrm = {};

ip = 1; % index of patients with signif ch pairs

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if isfile([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data])
            % load the PLV data with statistics
            load([basedir pacienti(p).folder '\' subfolder '\' PLV_folder '\' PLV_data]);
            
%             % Find rows (ch pairs) with at least 2 positive values (2 freq bins)
%             positive_values_count = sum(plv_signif_allPairs_clustcorr > 0, 2);
%             significant_chanPairs = find(positive_values_count >= 2);

            % Find rows with at least 2 consecutive positive values
            rows_with_consecutive_positives = [];
            
            for i = 1:size(plv_signif_allPairs_clustcorr, 1)
                row = plv_signif_allPairs_clustcorr(i, :);                
                % Use a sliding window to check for consecutive positive values
                for j = 1:(length(row)-1)
                    if row(j) > 0 && row(j+1) > 0
                        rows_with_consecutive_positives = [rows_with_consecutive_positives, i];
                        break; % Break once you find the first occurrence
                    end
                end
            end
            
            significant_chanPairs = rows_with_consecutive_positives';
            
            if isempty(significant_chanPairs) % if no significant ch pairs were found, skip this patient
                continue
            end
            
            % significant PLV diff
            PLVdiff_chpairs_signif = plv_signif_allPairs_clustcorr(significant_chanPairs, : );
            [max_PLVdiff, ichPair] = max(sum(PLVdiff_chpairs_signif,2)); % max PLV difference
            
            ichPairMax = significant_chanPairs(ichPair);
            idxChanMax = ROI_chanpairs(ichPairMax, :);  % channel indices in this pair in the original data
            
            PLVdelay_allSubj.plvspctrm{ip} = squeeze(PLVCond1.plvspctrm(idxChanMax(1), idxChanMax(2), :)); % ch x ch x freq
            PLVbs_allSubj.plvspctrm{ip} = squeeze(PLVCond2.plvspctrm(idxChanMax(1), idxChanMax(2), :)); % ch x ch x freq
            ip = ip+1;
            
            % also plot this chn pair
            significancePos = plv_signif_allPairs_clustcorr(ichPairMax,:)>0;
            figure; plot(PLVCond1.freq, squeeze(PLVCond1.plvspctrm(idxChanMax(1), idxChanMax(2), :)),'b', 'LineWidth', 1)
            hold on
            plot(PLVCond2.freq, squeeze(PLVCond2.plvspctrm(idxChanMax(1), idxChanMax(2), :)),'r', 'LineWidth', 1)
            hold on
            ylim([0 max(squeeze(PLVCond1.plvspctrm(idxChanMax(1), idxChanMax(2), :)))*1.05])
            yLimits = ylim;
            y = yLimits(2) * 0.1;
            
            % plot significance
            ySignifPosValues = ones(1, numel(PLVCond2.freq)) * NaN; % Initialize with NaN
            ySignifPosValues(significancePos) = y; % Assign y values only at significant indices
            plot(PLVCond2.freq, ySignifPosValues, '.-b', 'LineWidth', 3, 'MarkerSize', 10);
            xlim([PLVCond2.freq(1) PLVCond2.freq(end)])
            xlabel('Frequency, Hz')
            ylabel('PLV')
            legend(dataCond1.condition, dataCond2.condition);
            title(['PLV between ' PLVCond1.label{idxChanMax(1), 1} ' ' dataCond1.channelInfo(idxChanMax(1)).ROI ' and ' PLVCond1.label{idxChanMax(2), 1} ' ' dataCond1.channelInfo(idxChanMax(2)).ROI]);
            
        end
    end
end
PLVdelay_allSubj.freq = PLVCond1.freq;

%% now transform to the format suitable for ft_freqstatistics
data_arr{1}.powspctrm = PLVdelay_allSubj.plvspctrm{1}'; % delay
data_arr{2}.powspctrm = PLVbs_allSubj.plvspctrm{1}'; % baseline

for i = 2:(ip-1)
    data_arr{1}.powspctrm = [data_arr{1}.powspctrm; PLVdelay_allSubj.plvspctrm{i}'];
    data_arr{1}.label{1} = 'VTC-IPL delay';
    data_arr{1}.dimord = 'subj_freq';
    data_arr{1}.freq = PLVdelay_allSubj.freq;
    data_arr{2}.powspctrm = [data_arr{2}.powspctrm; PLVbs_allSubj.plvspctrm{i}'];
    data_arr{2}.label{1} = 'VTC-IPL baseline';
    data_arr{2}.freq = PLVdelay_allSubj.freq;
    data_arr{2}.dimord = 'subj_freq';
end

%% permutation statistic
cfg = [];
cfg.event_comparisons = {[1 2]};
stat_maint = statistics_fr(cfg, data_arr);

%%
sign_freq_bins = stat_maint{1}.VTC_IPL_delay.freq(stat_maint{1}.VTC_IPL_delay.mask);
if isempty(sign_freq_bins)
    disp('No significance frequency bands for delay > bs')
end

% plot mean, std and significance
n_subjects = ip-1;

figure;
ciplot(mean(data_arr{1}.powspctrm,1)-std(data_arr{1}.powspctrm,[],1)/sqrt(n_subjects), ...
    mean(data_arr{1}.powspctrm,1)+std(data_arr{1}.powspctrm,[],1)/sqrt(n_subjects), data_arr{1}.freq, [0.8, 0.8, 1]) ; % std err of mean in delay
hold on
ciplot(mean(data_arr{2}.powspctrm,1)-std(data_arr{2}.powspctrm,[],1)/sqrt(n_subjects), ...
    mean(data_arr{2}.powspctrm,1)+std(data_arr{2}.powspctrm,[],1)/sqrt(n_subjects), data_arr{1}.freq, [1, 0.8, 0.8]) ; % std err of mean in bs
hold on
plot(data_arr{1}.freq, mean(data_arr{1}.powspctrm,1),'b', 'LineWidth', 1) % mean of delay over subj
hold on
plot(data_arr{2}.freq, mean(data_arr{2}.powspctrm,1),'r', 'LineWidth', 1) % mean of bs over subj


% ylim([0 maxplv*1.05])
% yLimits = ylim;
% y = yLimits(2) * 0.1;
%
% % plot significance
% ySignifPosValues = ones(1, numel(PLVCond2.freq)) * NaN; % Initialize with NaN
% ySignifPosValues(significancePos) = y; % Assign y values only at significant indices
% plot(PLVCond2.freq, ySignifPosValues, '.-b', 'LineWidth', 3, 'MarkerSize', 10);
% hold on
%
% ySignifNegValues = ones(1, numel(PLVCond2.freq)) * NaN; % Initialize with NaN
% ySignifNegValues(significanceNeg) = y; % Assign y values only at significant indices
% plot(PLVCond2.freq, ySignifNegValues, '.-r', 'LineWidth', 3, 'MarkerSize', 10);

% xlim([PLVCond2.freq(1) PLVCond2.freq(end)])
xlabel('Frequency, Hz')
ylabel('PLV')
%                 legend('delay last 1.5 sec', 'baseline 1.5 sec');
legend(dataCond1.condition, dataCond2.condition);
% title(['PLV between ' PLVCond1.label{ipair(1), 1} ' ' dataCond1.channelInfo(ipair(1)).ROI ' and ' PLVCond1.label{ipair(2), 1} ' ' dataCond1.channelInfo(ipair(2)).ROI]);
title('PLV between VTC and IPL, mean over 5 patients')
