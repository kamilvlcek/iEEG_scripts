function Granger_permut_stat(ROI1, ROI2, condition, period, significant, fsample, n_permutes, threshold, stat)
% does Granger analysis between ROI1 and ROI2 for ch pairs that have significant PLV or all
% possibility to use permutation stat (shuffling trials in one channel), compares 2 directions: if chn A -> cnh B is greater than chn B -> chn A
% first, set up pacienti_memact by selecting which patients to analyze

%%% these parameters should be configured each time:
% ROI1 - e.g. 'VTC'
% ROI2 - e.g. 'IPL'
% condition - can be 'same' , 'diff' or 'all' (trials of both conditions)
% period - which period of the task to use
%   1 - last 1.9 s of delay (delay 2) 
%   2 - first 1.9 s of delay (delay 1)
%   3 - encoding 1.9 s
%   4 - baseline 

%%% optional:
% significant - if 1, selects only significant channel pairs from PLV analysis; % default all ROI1-ROI2 chan pairs
% fsample - new sampling frequency for downsampling (default: 40 Hz (then freq range would be [0 20]))
% n_permutes - number of permutations (default: 200)
% threshold - alpha level for p-value (default: 0.05)
% stat - if 0, compute Granger without stats; if 1, run permutation stats (default: 0)
% --------------------
% by Sofiia Moraresku
% September 2024
% --------------------

tic
if(~exist('condition','var')) || isempty(condition), condition = 'all'; end % default - all trials
if(~exist('significant','var')) || isempty(significant), significant = 0; end % default all ROI1-ROI2 chan pairs
if(~exist('fsample','var')) || isempty(fsample), fsample = 40; end % default resample to 40 Hz
if(~exist('n_permutes','var')) || isempty(n_permutes), n_permutes = 200; end % can be increased for more stable results
if(~exist('threshold','var')) || isempty(threshold), threshold = 0.05; end
if(~exist('stat','var')) || isempty(stat), stat = 0; end % default - without statistics

ft_defaults % set up fieldtrip

data_preproc = 'MemAct_refBipo -2.0-7.9 adel fieldtrip_2024-04.mat';

setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact, some are marked not to analyze as they don't have ROIs we want

for p = 1:numel(pacienti)
    if pacienti(p).todo
        % load the prepared file in fieldtrip format
        patient_path = [basedir pacienti(p).folder '\' subfolder];
        load([patient_path '\' data_preproc]);
        
%         patient_id = regexp(pacienti(p).folder,' ','split');
%         patient_id = patient_id{1};
        patient_id = pacienti(p).folder; % in the table with all implanted channels, a full name of patient (as the name of the folder) is used
        
        % select channels, trials and period for the analysis
        [chan_labels, ROI_labels, ROI_chanpairs] = select_chan_ROI_pair([], ROI1, ROI2, patient_id, significant, patient_path); % select channels for this patient and ROIs
            
            if isempty(ROI_chanpairs)
                continue; % if for this patient no ch pairs in 2 ROIs were found, switch to the next
            end
            
            itrials_good = get_correct_trials(condition,0, TrialInformationTable); % good trials for this condition or for both conditions
            cfg = [];
            cfg.channel = chan_labels;
            cfg.trials = itrials_good;
            dataROI = ft_selectdata(cfg,data);
            
            %  period of the task to analyze
            pcfg = [];
            switch period
                case 1
                    pcfg.latency = [4, 5.9-1/data.fsample]; % last 1.9 sec of delay
                    str_period = 'last 1.9s delay';
                case 2
                    pcfg.latency = [2.1+1/data.fsample, 4];    % first 1.9 sec of delay
                    str_period = 'first 1.9s delay';
                case 3
                    pcfg.latency = [0.1+1/data.fsample, 2]; % 1.9s encoding
                    str_period = '1.9s encoding';
                case 4
                    pcfg.latency = [-1.9, -1/data.fsample]; % 1.9s baseline
                    str_period = '1.9s baseline';
            end
        
      dataPeriod = ft_selectdata(pcfg, dataROI);
      dataPeriod.period = str_period;
        % save ROI labels to the data
        for i = 1:numel(dataPeriod.channelInfo)
            dataPeriod.channelInfo(i).ROI = ROI_labels{i};
        end
        
        %% granger analysis
        % first downsample
        cfg = [];
        cfg.resamplefs = fsample;
        dataPeriodResampled = ft_resampledata(cfg, dataPeriod);
        
        % real Granger for delay
        freqcfg             = [];
        freqcfg.method      = 'mtmfft';
        freqcfg.taper       = 'dpss';
        freqcfg.output      = 'fourier';
        freqcfg.tapsmofrq   = 2; % the amount of spectral smoothing through multi-tapering, Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
        freqcfg.pad         = 20; % the same as in Dimakopoulos et al, 2022 elife
        freqcfg.keeptrials = 'yes';
        freqPeriod       = ft_freqanalysis(freqcfg, dataPeriodResampled);
        
        channelcmbToRun = [chan_labels(ROI_chanpairs(:,1)) chan_labels(ROI_chanpairs(:,2))]; % labels of chan in pairs
        grangercfg = [];
        grangercfg.channelcmb = channelcmbToRun;
        grangercfg.method  = 'granger';
        grangercfg.granger.sfmethod = 'bivariate';
        grangercfg.conditional = 'no';
        grangerPeriod = ft_connectivityanalysis(grangercfg, freqPeriod); % grangerPeriod.grangerspctrm  - (n_chpairs*2) x freq
        
        %% permutat stat
        % by shuffling trials in one channel
        if stat
        % initialize to store the output for this patient
        Granger_signif_allPairs = zeros(size(ROI_chanpairs,1), numel(grangerPeriod.freq)); % chn pairs x significant Granger diff
        p_values_allPairs = zeros(size(ROI_chanpairs,1), numel(grangerPeriod.freq)); % chn pairs x p values for Granger diff
        Granger_signif_allPairs_clustcorr = zeros(size(ROI_chanpairs,1), numel(grangerPeriod.freq)); % chn pairs x significant Granger diff after cluster correction
        
        nChannel_ToChange = 2; % to mix trials in the second channel (in one pair)
        iPairGranger = 1; % index of ch pair in granger spectrum (as each chan pair has two rows of values)
        
        for iPair = 1:size(ROI_chanpairs,1)
                       
            % Channel indexes
            iChannel_1 = ROI_chanpairs(iPair,1);
            iChannel_2 = ROI_chanpairs(iPair,2);
            
            fprintf('Calculating %s with  %s right now. Elapsed pairs for this patient are %d..... \n',dataPeriod.label{iChannel_1, 1},dataPeriod.label{iChannel_2, 1}, size(ROI_chanpairs,1)-iPair);
            
            % select data for one pair of chan
            cfg = [];
            cfg.channel = [iChannel_1, iChannel_2];
            dataChPair = ft_selectdata(cfg, dataPeriodResampled);
            
            % initialize for Granger distribution of null hypothesis values
            permuted_vals_diff = zeros(n_permutes, numel(grangerPeriod.freq)); % n_permutes x freq;
            
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
            realGranger_diff = grangerPeriod.grangerspctrm(iPairGranger,:)-grangerPeriod.grangerspctrm(iPairGranger+1,:); % frequency x 1
            
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
        end
        
        % Save results
        strVariableFolder = [ basedir pacienti(p).folder '\' subfolder '\Granger_permut_stat\'];
        mkdir(strVariableFolder)

        if stat
            save([strVariableFolder,'Granger_' ROI1 '-' ROI2 '_' str_period '_' condition '_trials_' num2str(n_permutes) 'permut_', datestr(now,'YYYY-mm'), '.mat'],...
                'ROI_chanpairs','Granger_signif_allPairs', 'Granger_signif_allPairs_clustcorr','p_values_allPairs', 'grangerPeriod','freqPeriod','dataPeriod', 'dataPeriodResampled',  'n_permutes', 'threshold', '-v7.3')
        else
            if significant == 0
                str_sign = 'all_pairs_';
            else
                str_sign = 'PLV_sign_pairs_';
            end
            save([strVariableFolder,'Granger_' ROI1 '-' ROI2 '_' str_period '_' condition '_trials_' str_sign, datestr(now,'YYYY-mm'), '.mat'],...
                'ROI_chanpairs', 'grangerPeriod','freqPeriod','dataPeriod', 'dataPeriodResampled', '-v7.3')
        end
        
    end
end
toc
end