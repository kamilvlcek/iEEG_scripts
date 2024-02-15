function PLV_permutation_stat(ROI1, ROI2, condition, period, significant, freq, n_permutes, threshold)
%%% computes PLV between ROI1 and ROI2 in each patient for channels that showed alpha power increase during the delay
%%% computes permutation statistics by shuffling conditions, either periods - comparing delay vs baseline, or comparing real conditions - same and diff during the delay
%%% Code adapted from M. X. Cohen, 2014 (https://github.com/mikexcohen/AnalyzingNeuralTimeSeries)

% first, set up pacienti_memact by selecting which patients to analyze
%%% these parametrs should be configured each time:
% ROI1 - e.g. 'LOC'
% ROI2 - e.g. 'IPL'
% condition - can be 'same' , 'diff' or 'all' (trials of both conditions) or empty for the case of period = 0
% period - which period to compare with which:
% 0 - compare 2 conditions during the delay
% 1 - last 1.5 sec of delay (4.4-5.9) vs 1.5 sec baseline
% 2 - middle 1.5 sec of delay (2.9-4.4) vs 1.5 sec baseline
% 3 last 1.5 sec of delay (4.4-5.9) vs 1.5 sec of encoding
% 4 encoding vs baseline (not in the code yet)
% significant - only for period =0; if 1, select only significant chan pairs that were obtained by PLV_permutation_stat for all trials
%%% optional:
% freq - freq range for which compute PLV
% n_permutes - number of permutations, default = 200
% threshold - alpha level for p-value, default = 0.05

tic
if(~exist('freq','var')) || isempty(freq), freq = [2:40]; end % default 2-40 Hz for PLV
if(~exist('n_permutes','var')) || isempty(n_permutes), n_permutes = 200; end
if(~exist('threshold','var')) || isempty(threshold), threshold = 0.05; end

ft_defaults % set up fieldtrip

data_preproc = 'MemAct_refBipo -1.5-5.9 bdel fieldtrip_2024-02.mat';
% data_preproc ='MemAct_refBipo -2.0-5.9 bdel fieldtrip_2024-02.mat';

setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact, some are marked not to analyze as they don't have ROIs we want

for p = 1:numel(pacienti)
    if pacienti(p).todo
        % load the prepared file in fieldtrip format
        patient_path = [basedir pacienti(p).folder '\' subfolder];
        load([patient_path '\' data_preproc]);
        
        patient_id = regexp(pacienti(p).folder,' ','split');
        patient_id = patient_id{1};
        
        % select channels, trials and period for the analysis
        if period == 0 % if we want to compare two conditions for one period (delay)
            [chan_labels, ROI_labels, ROI_chanpairs] = select_chan_ROI_pair([], ROI1, ROI2, patient_id, significant, patient_path); % select only significant chan pairs (that were significant for all trials)
            
            itrials_same = get_correct_trials('same', TrialInformationTable); % good trials for same condition
            cfg = [];
            cfg.channel = chan_labels; % select channels in 2 ROIs
            cfg.trials = itrials_same; % for now, we don't reject individual epochs for each channel with spikes
            cfg.latency = [4.4, 5.9+1/data.fsample];    % last 1.5 sec of delay, if fsample=512, 768 time points
            dataCond1 = ft_selectdata(cfg,data);
            dataCond1.condition = 'same';
            
            itrials_diff = get_correct_trials('diff', TrialInformationTable); % good trials for diff condition
            cfg = [];
            cfg.channel = chan_labels; % select channels in 2 ROIs
            cfg.trials = itrials_diff;
            cfg.latency = [4.4, 5.9+1/data.fsample];    % last 1.5 sec of delay, if fsample=512, 768 time points
            dataCond2 = ft_selectdata(cfg,data);
            dataCond2.condition = 'diff';
            str2save = '_diff_vs_same_last 1.5s delay_';
            
        else  % if we want to compare two periods for one condition
            % find channels for this patient and ROIs in the table if there are some
            [chan_labels, ROI_labels, ROI_chanpairs] = select_chan_ROI_pair([], ROI1, ROI2, patient_id); % select channels for this patient and ROIs
            
            if numel(unique(ROI_labels)) < 2
                continue; % if for this patient no channels in 2 ROIs were found, switch to the next
            end
            
            itrials_good = get_correct_trials(condition, TrialInformationTable); % good trials for this condition or for both conditions
            cfg = [];
            cfg.channel = chan_labels;
            cfg.trials = itrials_good;
            dataROI = ft_selectdata(cfg,data);
            
            %% select delay and baseline period (or encoding)
            if period == 1
                delaycfg = [];
                delaycfg.latency = [4.4, 5.9+1/data.fsample];    % last 1.5 sec of delay, if fsample=512, 768 time points
%                 delaycfg.latency = [3.9+1/dataROI.fsample, 5.9]; % last 2 sec of delay
                bscfg = []; 
                bscfg.latency = [-1.5, -1/dataROI.fsample]; % baseline [-1.5 0], 768 time points
%                 bscfg.latency = [-2, -1/dataROI.fsample]; % baseline 2 sec
                str2save = '_last 1.5s delay_vs_bs_';
                str_period = 'baseline';
            elseif period == 2
                delaycfg = [];
                delaycfg.latency = [2.9+1/dataROI.fsample, 4.4];    % middle 1.5 sec of delay, if fsample=512, 768 time points
%                 delaycfg.latency = [1/dataROI.fsample, 2];   % first 2 sec of delay
                bscfg = [];
                bscfg.latency = [-1.5, -1/dataROI.fsample];
%                 bscfg.latency = [-2, -1/dataROI.fsample]; % baseline 2 sec
                str2save = '_middle 1.5s delay_vs_bs_';
%                 str2save = '_first 2s delay_vs_2s_bs_';
                str_period = 'baseline';
            elseif period == 3
                delaycfg = [];
                delaycfg.latency = [4.4, 5.9+1/data.fsample];    % last 1.5 sec of delay
                bscfg = [];
                bscfg.latency = [1/dataROI.fsample, 1.5]; % encoding 1.5 sec
                str2save = '_last 1.5s delay_vs_1.5s encoding_';
                str_period = 'encoding';
            end
            
            dataCond1 = ft_selectdata(delaycfg,dataROI); % delay
            dataCond1.condition = 'delay';
            dataCond2 = ft_selectdata(bscfg,dataROI); % baseline
            dataCond2.condition = str_period;
        end
        % save ROI labels for each condition
        for i = 1:numel(dataCond1.channelInfo)
            dataCond1.channelInfo(i).ROI = ROI_labels{i};
            dataCond2.channelInfo(i).ROI = ROI_labels{i};
        end
        
        %% Frequency and PLV analysis for true delay and baseline period (or for two conditions)
        freqcfg             = [];
        freqcfg.method      = 'mtmfft';
        freqcfg.taper       = 'dpss';
        freqcfg.output      = 'fourier';
        freqcfg.foi         = freq;
        freqcfg.tapsmofrq   = 2; % the amount of spectral smoothing through multi-tapering, Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
        freqcfg.pad         = 2;
        freqcfg.keeptrials = 'yes';
        FreqCond1   = ft_freqanalysis(freqcfg, dataCond1);
        
        plvcfg             = [];
        plvcfg.method      = 'plv';
        PLVCond1  = ft_connectivityanalysis(plvcfg, FreqCond1);
        
        % Frequency and PLV analysis for baseline period
        FreqCond2   = ft_freqanalysis(freqcfg, dataCond2);
        PLVCond2  = ft_connectivityanalysis(plvcfg, FreqCond2);
        
        %% permutation stat
        
        % initialize to store the output for this patient
        plv_signif_allPairs = zeros(size(ROI_chanpairs,1), numel(PLVCond1.freq)); % chn pairs x significant PLV diff
        p_values_allPairs = zeros(size(ROI_chanpairs,1), numel(PLVCond1.freq)); % chn pairs x p values for PLV diff
        plv_signif_allPairs_clustcorr = zeros(size(ROI_chanpairs,1), numel(PLVCond1.freq)); % chn pairs x significant PLV diff after cluster correction
        
        % n of trials in real delay and baseline data (or in each condition)
        ntrialsCond1 = length(dataCond1.trial);     % for VT62 - 132 trials
        ntrialsCond2 = length(dataCond2.trial); % for VT62 - 132 trials
        
        % concatenate trials from two conditions
        cfg = [];
        data2periods = ft_appenddata(cfg,dataCond1,dataCond2); % for VT62 - 264 trials
        
        fprintf('Computing permutat stat for patient %s \n', patient_id);
        for iPair = 1:size(ROI_chanpairs,1)
            
            % Channel indexes
            iChannel_1 = ROI_chanpairs(iPair,1);
            iChannel_2 = ROI_chanpairs(iPair,2);
            
            fprintf('Calculating %s with  %s right now. Elapsed pairs for this patient are %d..... \n',PLVCond1.label{iChannel_1, 1},PLVCond1.label{iChannel_2, 1}, size(ROI_chanpairs,1)-iPair);
            
            % select data for one pair of chan
            cfg = [];
            cfg.channel = [iChannel_1, iChannel_2];
            data2periodsChPair = ft_selectdata(cfg, data2periods);
            PLVCond1ChPair = ft_selectdata(cfg, PLVCond1);
            PLVCond2ChPair = ft_selectdata(cfg, PLVCond2);
            
            % initialize for PLV distribution of null hypothesis values
            permuted_vals_diff = zeros(n_permutes, numel(PLVCond1.freq)); % n_permutes x freq; PLV differences between 2 periods (conditions)
            
            % permute
            for permi=1:n_permutes
                
                % random permutation
                fakeconds = randperm(length(data2periods.trial));
                
                % shuffled condition labels
                fakeconds(fakeconds<ntrialsCond1+1) = 1;
                fakeconds(fakeconds>1) = 2;
                
                % Frequency and PLV analysis for fake condtition 1
                cfg = [];
                cfg.trials = fakeconds==1;
                dataFakeCond1 = ft_selectdata(cfg,data2periodsChPair);
                FreqFakeCond1   = ft_freqanalysis(freqcfg, dataFakeCond1);
                PLVFakeCond1  = ft_connectivityanalysis(plvcfg, FreqFakeCond1);
                
                % Frequency and PLV analysis for fake condtition 2
                cfg = [];
                cfg.trials = fakeconds==2;
                dataFakeCond2 = ft_selectdata(cfg,data2periodsChPair);
                FreqFakeCond2   = ft_freqanalysis(freqcfg, dataFakeCond2);
                PLVFakeCond2  = ft_connectivityanalysis(plvcfg, FreqFakeCond2);
                
                % PLV difference of 2 spectrums of shuffled conditions
                permuted_vals_diff(permi,:) = squeeze(PLVFakeCond1.plvspctrm(1,2,:)) - squeeze(PLVFakeCond2.plvspctrm(1,2,:));
                
            end
            
            % % plot the distribution of H0 values for one freq (just for checking purposes)
            % figure(1), clf
            % histogram(permuted_vals_diff(:,10),50)
            
            % real PLV difference between 2 periods
            realPLV_diff = squeeze(PLVCond1ChPair.plvspctrm(1,2,:))-squeeze(PLVCond2ChPair.plvspctrm(1,2,:)); % frequency x 1
            
            % transform to z values and set to zero values below thershold
            zmap = (realPLV_diff'-mean(permuted_vals_diff,1))./std(permuted_vals_diff);
            plv_diff_signif = realPLV_diff;
            p_value_difference = 2 * (1 - normcdf(abs(zmap))); % compute two-tailed p-values using normcdf (each z value is converted to p value)
            plv_diff_signif(p_value_difference > threshold) = 0;
            plv_signif_allPairs(iPair, :) = plv_diff_signif'; % save significant plv for this chan pair
            p_values_allPairs(iPair, :) = p_value_difference;
            
            %% the cluster correction on the permuted data
            for permi = 1:n_permutes
                
                % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
                fakecorrsz = (permuted_vals_diff(permi,:)-mean(permuted_vals_diff,1)) ./ std(permuted_vals_diff,[],1) ;
                %                 fakecorrsz(abs(fakecorrsz)<norminv(1-threshold))=0;
                fakecorrsz(2 * (1 - normcdf(abs(fakecorrsz))) > threshold) = 0;
                
                % get number of elements in largest supra-threshold cluster
                clustinfo = bwconncomp(fakecorrsz);
                max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
                % using cellfun here eliminates the need for a slower loop over cells
            end
            
            % apply cluster-level corrected threshold
            plv_diff_clust_corr = plv_diff_signif; % uncorrected pixel-level threshold
            
            % find islands and remove those smaller than cluster size threshold
            clustinfo = bwconncomp(plv_diff_clust_corr);
            clust_info = cellfun(@numel,clustinfo.PixelIdxList);
            clust_threshold = prctile(max_clust_info,100-threshold*100);
            
            % identify clusters to remove
            whichclusters2remove = find(clust_info<clust_threshold);
            
            % remove clusters
            for i=1:length(whichclusters2remove)
                plv_diff_clust_corr(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
            end
            
            plv_signif_allPairs_clustcorr(iPair, :) = plv_diff_clust_corr'; % save significant plv after cluster corr for this chan pair
        end
        
        % Save results
        strVariableFolder = [ basedir pacienti(p).folder '\' subfolder '\PLV_permut_stat\'];
        mkdir(strVariableFolder)
        %         save([strVariableFolder,'PLV_' ROI1 '-' ROI2 str2save condition '_trials_', datestr(now,'YYYY-mm'), '.mat'],...
        %             'ROI_chanpairs','plv_signif_allPairs', 'plv_signif_allPairs_clustcorr','p_values_allPairs', 'PLVDelay', 'PLVBaseline','dataDelay', 'dataBaseline', 'n_permutes', 'threshold', '-v7.3')
        save([strVariableFolder,'PLV_' ROI1 '-' ROI2 str2save condition '_trials_', datestr(now,'YYYY-mm'), '.mat'],...
            'ROI_chanpairs','plv_signif_allPairs', 'plv_signif_allPairs_clustcorr','p_values_allPairs', 'PLVCond1', 'PLVCond2','dataCond1', 'dataCond2', 'n_permutes', 'threshold', '-v7.3')
        
        
    end
end
toc
end
