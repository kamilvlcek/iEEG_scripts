%% set up a path for data
granger_data = 'Granger_VTC-IPL last 2s delay all_trials_2024-02.mat'; 
granger_folder = 'Granger_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact();

freq_start = 2; % take freq from 2 Hz to 15 Hz
freq_end = 20; 
%%
% initialize 2 cells, one - Granger for one direction and for opposite direction
Granger_allSubj1 = [];
Granger_allSubj1.grangerspctrm = {};

Granger_allSubj2 = [];
Granger_allSubj2.grangerspctrm = {};

ip = 1; % index of patients with signif ch pairs

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if isfile([basedir pacienti(p).folder '\' subfolder '\' granger_folder '\' granger_data])
            % load the PLV data with statistics
            load([basedir pacienti(p).folder '\' subfolder '\' granger_folder '\' granger_data]);
            
            % first find indexes of chan pairs with significant direction VTC -> IPL (negative difference)
%             significant_chanPairs = find(sum(Granger_signif_allPairs_clustcorr,2) < 0);
            significant_chanPairs = find(any(Granger_signif_allPairs_clustcorr < 0,2) &...
                ~any(Granger_signif_allPairs_clustcorr > 0, 2)); % at least one negative value and no pisitive values

            iPairsGranger = 1:2:2*size(ROI_chanpairs, 1); % as we have 2 row of values for each ch pair
            
            % significant Granger diff
            Granger_diff_chpairs_signif = Granger_signif_allPairs_clustcorr(significant_chanPairs, : );
            [max_Granger_diff, ichPair] = max(abs(sum(Granger_diff_chpairs_signif,2))); % max Granger difference
            
            ichPairMax = significant_chanPairs(ichPair);
            idxChanMax = ROI_chanpairs(ichPairMax, :);  % channel indices in this pair in the original data
            
            % indexes of freq we want to take
            ifreq = find(grangerDelay.freq == freq_start);
            ifreq2 = find(grangerDelay.freq == freq_end);
            
            Granger_allSubj1.grangerspctrm{ip} = grangerDelay.grangerspctrm(iPairsGranger(ichPairMax),ifreq:ifreq2); % chpair x freq
            Granger_allSubj2.grangerspctrm{ip} = grangerDelay.grangerspctrm(iPairsGranger(ichPairMax)+1,ifreq:ifreq2); % chpair x freq
            ip = ip+1;
            
        end
    end
end
            
Granger_allSubj1.freq = grangerDelay.freq(ifreq : ifreq2);

%% now transform to the format suitable for ft_freqstatistics
data_arr{1}.powspctrm = Granger_allSubj1.grangerspctrm{1}; % one direction
data_arr{2}.powspctrm = Granger_allSubj2.grangerspctrm{1}; % opposite direction

for i = 2:(ip-1)
    data_arr{1}.powspctrm = [data_arr{1}.powspctrm; Granger_allSubj1.grangerspctrm{i}];
    data_arr{1}.label{1} = 'IPL->VTC';
    data_arr{1}.dimord = 'subj_freq';
    data_arr{1}.freq = Granger_allSubj1.freq;
    data_arr{2}.powspctrm = [data_arr{2}.powspctrm; Granger_allSubj2.grangerspctrm{i}];
    data_arr{2}.label{1} = 'VTC->IPL';
    data_arr{2}.freq = Granger_allSubj1.freq;
    data_arr{2}.dimord = 'subj_freq';
end
            
            
%% permutation statistic
cfg = [];
cfg.event_comparisons = {[1 2]};
stat_delay = statistics_fr(cfg, data_arr);

sign_freq_bins = stat_delay{1}.VTC_IPL_delay.freq(stat_delay{1}.VTC_IPL_delay.mask);
if isempty(sign_freq_bins)
    disp('No significance frequency bands for difference between 2 directions')
end            

%% plot mean, std and significance
n_subjects = ip-1;

figure;
ciplot(mean(data_arr{1}.powspctrm,1)-std(data_arr{1}.powspctrm,[],1)/sqrt(n_subjects), ...
    mean(data_arr{1}.powspctrm,1)+std(data_arr{1}.powspctrm,[],1)/sqrt(n_subjects), data_arr{1}.freq, [0.8, 0.8, 1]) ; % std err of mean in delay
hold on
ciplot(mean(data_arr{2}.powspctrm,1)-std(data_arr{2}.powspctrm,[],1)/sqrt(n_subjects), ...
    mean(data_arr{2}.powspctrm,1)+std(data_arr{2}.powspctrm,[],1)/sqrt(n_subjects), data_arr{1}.freq, [1, 0.8, 0.8]) ; % std err of mean in bs
hold on
h1 = plot(data_arr{1}.freq, mean(data_arr{1}.powspctrm,1),'b', 'LineWidth', 1); % mean of IPL->VTC over subj
hold on
h2 = plot(data_arr{2}.freq, mean(data_arr{2}.powspctrm,1),'r', 'LineWidth', 1); % mean of VTC->IPL over subj
hold on
% plot significance
yLimits = ylim;
y = yLimits(2) * 0.05;
ySignifValues = ones(1, numel(data_arr{1}.freq)) * NaN; % Initialize with NaN
ySignifValues(stat_delay{1}.VTC_IPL_delay.mask) = y; % Assign y values only at significant indices
plot(data_arr{1}.freq, ySignifValues, '.-r', 'LineWidth', 3, 'MarkerSize', 10);

xlabel('Frequency, Hz')
ylabel('Granger')
legend([h1, h2], {'IPL->VTC', 'VTC->IPL'});
title('Granger between VTC and IPL, mean over 5 patients')
