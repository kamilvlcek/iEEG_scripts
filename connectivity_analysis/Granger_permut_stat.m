function Granger_permut_stat(ROI1, ROI2, condition, fsample, n_permutes, threshold)
% does Granger analysis between ROI1 and ROI2 for the last 1.5 sec of delay for ch pairs that have significant PLV last 1.5s delay > bs
% using permut stat (shuffling trials in one channel), compares 2 directions: if chn A -> cnh B is greater than chn B -> chn A
% first, set up pacienti_memact by selecting which patients to analyze

%%% these parametrs should be configured each time:
% ROI1 - e.g. 'VTC'
% ROI2 - e.g. 'IPL'
% condition - can be 'same' , 'diff' or 'all' (trials of both conditions)
%%% optional:
% fsample - new sample frequency (for Granger we need to downsample), default 40 Hz (then freq range would be [0 20])
% n_permutes - number of permutations, default = 200
% threshold - alpha level for p-value, default = 0.05

tic
if(~exist('condition','var')) || isempty(condition), condition = 'all'; end % default - all trials
if(~exist('fsample','var')) || isempty(fsample), fsample = 40; end % default resample to 40 Hz
if(~exist('n_permutes','var')) || isempty(n_permutes), n_permutes = 200; end
if(~exist('threshold','var')) || isempty(threshold), threshold = 0.05; end

ft_defaults % set up fieldtrip

data_preproc = 'MemAct_refBipo -1.5-5.9 bdel fieldtrip_2024-02.mat';

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
        [chan_labels, ROI_labels, ROI_chanpairs] = select_chan_ROI_pair([], ROI1, ROI2, patient_id, 1, patient_path); % select only significant chan pairs with PLV delay > bs
        
        if isempty(ROI_chanpairs) % if no significant ch pairs were found, skip this patient
            continue
        end
        
        itrials_good = get_correct_trials(condition, TrialInformationTable); % good trials for this condition or for both conditions
        cfg = [];
        cfg.channel = chan_labels;
        cfg.trials = itrials_good;
        cfg.latency = [4.4, 5.9+1/data.fsample];    % last 1.5 sec of delay, if fsample=512, 768 time points
        dataDelay = ft_selectdata(cfg,data);
        
        % save ROI labels to the data
        for i = 1:numel(dataDelay.channelInfo)
            dataDelay.channelInfo(i).ROI = ROI_labels{i};
        end
        
        %% granger analysis
        % first downsample
        cfg = [];
        cfg.resamplefs = fsample;
        dataDelayResampled = ft_resampledata(cfg, dataDelay);
        
        % real Granger for delay
        freqcfg             = [];
        freqcfg.method      = 'mtmfft';
        freqcfg.taper       = 'dpss';
        freqcfg.output      = 'fourier';
        freqcfg.tapsmofrq   = 2; % the amount of spectral smoothing through multi-tapering, Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
        freqcfg.pad         = 10;
        freqcfg.keeptrials = 'yes';
        freqDelay       = ft_freqanalysis(freqcfg, dataDelayResampled);
        
        channelcmbToRun = [chan_labels(ROI_chanpairs(:,1)) chan_labels(ROI_chanpairs(:,2))]; % labels of chan in pairs
        grangercfg = [];
        grangercfg.channelcmb = channelcmbToRun;
        grangercfg.method  = 'granger';
        grangercfg.granger.sfmethod = 'bivariate';
        grangercfg.conditional = 'no';
        grangerDelay = ft_connectivityanalysis(grangercfg, freqDelay); % grangerDelay.grangerspctrm  - (n_chpairs*2) x freq
        
        %% permutat stat
        % by shuffling trials in one channel
        
        % initialize to store the output for this patient
        Granger_signif_allPairs = zeros(size(ROI_chanpairs,1), numel(grangerDelay.freq)); % chn pairs x significant Granger diff
        p_values_allPairs = zeros(size(ROI_chanpairs,1), numel(grangerDelay.freq)); % chn pairs x p values for Granger diff
        Granger_signif_allPairs_clustcorr = zeros(size(ROI_chanpairs,1), numel(grangerDelay.freq)); % chn pairs x significant Granger diff after cluster correction
        
        nChannel_ToChange = 2; % to mix trials in the second channel (in one pair)
        iPairGranger = 1; % index of ch pair in granger spectrum (as each chan pair has two rows of values)
        
        for iPair = 1:size(ROI_chanpairs,1)
                       
            % Channel indexes
            iChannel_1 = ROI_chanpairs(iPair,1);
            iChannel_2 = ROI_chanpairs(iPair,2);
            
            fprintf('Calculating %s with  %s right now. Elapsed pairs for this patient are %d..... \n',dataDelay.label{iChannel_1, 1},dataDelay.label{iChannel_2, 1}, size(ROI_chanpairs,1)-iPair);
            
            % select data for one pair of chan
            cfg = [];
            cfg.channel = [iChannel_1, iChannel_2];
            dataChPair = ft_selectdata(cfg, dataDelayResampled);
            
            % initialize for Granger distribution of null hypothesis values
            permuted_vals_diff = zeros(n_permutes, numel(grangerDelay.freq)); % n_permutes x freq;
            
            % permute
            for permi=1:n_permutes
                
                % random permutation
                random_trials = randperm(length(dataChPair.trial));
                dataMixedTrials = dataChPair;
                
                % mixed trials in the second channel
                for nTrial = 1:length(random_trials)
                    dataMixedTrials.trial{nTrial}(nChannel_ToChange,:) = dataChPair.trial{random_trials(nTrial)}(nChannel_ToChange,:);
                end
                
                % Granger for mixed trials 
                freqDelay_random       = ft_freqanalysis(freqcfg, dataMixedTrials);                
                grangercfg.channelcmb = [];
                grangerDelay_random = ft_connectivityanalysis(grangercfg, freqDelay_random);
                
                % difference of 2 spectrums of 2 directions
                permuted_vals_diff(permi,:) = grangerDelay_random.grangerspctrm(1,:) - grangerDelay_random.grangerspctrm(2,:);
                
            end
            
            % % plot the distribution of H0 values for one freq (just for checking purposes)
            % figure(1), clf
            % histogram(permuted_vals_diff(:,10),20)
            
            % real Granger difference between 2 directions
            realGranger_diff = grangerDelay.grangerspctrm(iPairGranger,:)-grangerDelay.grangerspctrm(iPairGranger+1,:); % frequency x 1
            
            % transform to z values and set to zero values below thershold
            zmap = (realGranger_diff-mean(permuted_vals_diff,1))./std(permuted_vals_diff);
            Granger_diff_signif = realGranger_diff;
            p_value_difference = 2 * (1 - normcdf(abs(zmap))); % compute two-tailed p-values using normcdf (each z value is converted to p value)
            Granger_diff_signif(p_value_difference > threshold) = 0;
            Granger_signif_allPairs(iPair, :) = Granger_diff_signif; % save significant granger for this chan pair
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
            Granger_diff_clust_corr = Granger_diff_signif; % uncorrected pixel-level threshold
            
            % find islands and remove those smaller than cluster size threshold
            clustinfo = bwconncomp(Granger_diff_clust_corr);
            clust_info = cellfun(@numel,clustinfo.PixelIdxList);
            clust_threshold = prctile(max_clust_info,100-threshold*100);
            
            % identify clusters to remove
            whichclusters2remove = find(clust_info<clust_threshold);
            
            % remove clusters
            for i=1:length(whichclusters2remove)
                Granger_diff_clust_corr(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
            end
            
            Granger_signif_allPairs_clustcorr(iPair, :) = Granger_diff_clust_corr; % save significant granger after cluster corr for this chan pair           
        
            iPairGranger = iPairGranger+2; % next index of ch pair in granger spectrum 
        end
        
        % Save results
        strVariableFolder = [ basedir pacienti(p).folder '\' subfolder '\Granger_permut_stat\'];
        mkdir(strVariableFolder)
        save([strVariableFolder,'Granger_' ROI1 '-' ROI2 ' last 1.5s delay ' condition '_trials_', datestr(now,'YYYY-mm'), '.mat'],...
            'ROI_chanpairs','Granger_signif_allPairs', 'Granger_signif_allPairs_clustcorr','p_values_allPairs', 'grangerDelay', 'freqDelay','dataDelay', 'dataDelayResampled', 'n_permutes', 'threshold', '-v7.3')
        
    end
end
toc
end