function PLV_analysis_ROI(ROI1, ROI2, condition, freq)
%%% computes PLV between ROI1 and ROI2 in each patient for channels that showed alpha power increase to any condition during the delay
%%% here no statistic is computed
%%% plots and saves some figures to the patient's folder:
%%% - the heatmap of average PLV in alpha range (7-14 Hz) for all between ROI connections
%%% - PLV spectrum for delay (last 1.5 sec), encoding (1.5 sec) and  baseline (1.5 sec) in individual chan pairs that have average alpha PLV > 0.3

% first, set up pacienti_memact by selecting which patients to analyze
% ROI1 - e.g. 'LOC'
% ROI2 - e.g. 'IPL'
% condition - can be 'same' , 'diff' or 'all' (trials of both conditions)
% freq - freq range for which compute PLV

if(~exist('freq','var')) || isempty(freq), freq = [4:40]; end % default 4-40 Hz
ft_defaults % set up fieldtrip

data_preproc = 'MemAct_refBipo -1.5-5.9 bdel fieldtrip_2024-02.mat';
alphaIncrease = 'E:\work\PhD\MemoryActions\results\iEEG\delayed condition\January2024_9pat\Response2XLS_CM_Memact_CHilbert_8-13Hz_AnyResp_OTP_-0.5-5.9_refBipo_Ep2024-01_encod+del_CHMult_2024-01-25_13-01-46.xls';
alphaIncreaseChan = readtable(alphaIncrease, 'ReadRowNames',true); % table with all channels from all patients showing alpha increase during the delay

setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact, some are marked not to analyze as they don't have ROIs we want

for p = 1:numel(pacienti)
    if pacienti(p).todo
        % load the prepared file in fieldtrip format
        load([basedir pacienti(p).folder '\' subfolder '\' data_preproc]);
        
        % find channels for this patient and ROIs in the table if there are some
        patient_id = regexp(pacienti(p).folder,' ','split');
        patient_id = patient_id{1};
        ichanROI = find(strcmp(alphaIncreaseChan.pacient, patient_id) & ...
            (strcmp(alphaIncreaseChan.mybrainlabel, ROI1) | strcmp(alphaIncreaseChan.mybrainlabel, ROI2))); % indexes of channels in LOC and IPL with alpha increase in any condition
        chan_labels = strrep(alphaIncreaseChan.name(ichanROI), [patient_id ' '], ''); % find labels of these chan
        ROI_labels = alphaIncreaseChan.mybrainlabel(ichanROI); % store also ROI labels in the data
        
        if numel(unique(ROI_labels)) < 2
            continue; % if for this patient no channels in 2 ROIs were found, switch to the next
        end
        
        % select channels and trials for the analysis
        itrials_good = get_correct_trials(condition, TrialInformationTable); % good trials for this condition or for both conditions
        cfg = [];
        cfg.channel = chan_labels; % select channels in LOC and IPL
        cfg.trials = itrials_good; % for now, we don't reject individual epochs for each channel with spikes
        dataROI = ft_selectdata(cfg,data);
        
        % save ROI labels
        for i = 1:numel(dataROI.channelInfo)
            dataROI.channelInfo(i).ROI = ROI_labels{i};
        end
        
        % select separate periods
        cfg = [];
        cfg.latency = [2.4+1/dataROI.fsample, 3.9];    % last 1.5 sec of delay, if fsample=512, 768 time points 
        dataDelay = ft_selectdata(cfg,dataROI); % delay
        
        cfg = [];
        cfg.latency = [1/dataROI.fsample, 1.5];
        dataEncod = ft_selectdata(cfg,dataROI); % encoding 1.5 sec, 768 time points
        
        cfg = [];
        cfg.latency = [-1.5, -1/dataROI.fsample]; % baseline [-1.5 0], 768 time points
        dataBaseline = ft_selectdata(cfg,dataROI); 
        
        % Frequency and PLV analysis for delay period
        cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'dpss';
        cfg.output      = 'fourier';
        cfg.foi         = freq;
        cfg.tapsmofrq   = 2;
        cfg.pad         = 2;
        cfg.keeptrials = 'yes';
        FreqDelay   = ft_freqanalysis(cfg, dataDelay);
        
        cfg             = [];
        cfg.method      = 'plv';
        PLVDelay  = ft_connectivityanalysis(cfg, FreqDelay);
        
        % Frequency and PLV analysis for baseline period
        cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'dpss';
        cfg.output      = 'fourier';
        cfg.foi         = freq;
        cfg.tapsmofrq   = 2;
        cfg.pad         = 2;
        cfg.keeptrials = 'yes';
        FreqBaseline   = ft_freqanalysis(cfg, dataBaseline);
        
        cfg             = [];
        cfg.method      = 'plv';
        PLVBaseline  = ft_connectivityanalysis(cfg, FreqBaseline);
        
        % Frequency and PLV analysis for encoding period
        cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'dpss';
        cfg.output      = 'fourier';
        cfg.foi         = freq;
        cfg.tapsmofrq   = 2;
        cfg.pad         = 2;
        cfg.keeptrials = 'yes';
        FreqEncod   = ft_freqanalysis(cfg, dataEncod);
        
        cfg             = [];
        cfg.method      = 'plv';
        PLVEncod  = ft_connectivityanalysis(cfg, FreqEncod);
        
        % create the heatmap of PLV in alpha range for all between ROI connections
        cfg = [];
        cfg.frequency   = [7 14]; % select plv for alpha range 7-14 Hz and average
        cfg.avgoverfreq = 'yes';
        PLVDelayAlpha = ft_selectdata(cfg,PLVDelay);
        
        labels = {dataDelay.channelInfo.ROI}; % all ROI labels for channels in the data
        ROI1_chan_idx = find(contains(labels, ROI1)); % channel indexes for ROI1 - loc
        ROI2_chan_idx = find(contains(labels,ROI2)); % channel indexes for ROI2 -ipl
        
        figure(1); imagesc(PLVDelayAlpha.plvspctrm(ROI2_chan_idx, ROI1_chan_idx), [0 0.5])
        set(gca,'Xtick',[1 : numel(ROI1_chan_idx)],'XTickLabel',PLVDelayAlpha.label(ROI1_chan_idx))
        set(gca,'Ytick',[1 : numel(ROI2_chan_idx)],'YTickLabel',PLVDelayAlpha.label(ROI2_chan_idx))
        
        % Add row labels above the plot
        row_label_pos = 1:size(PLVDelayAlpha.plvspctrm(ROI2_chan_idx, ROI1_chan_idx), 1);
        text(repmat(size(data, 2) + 0.2, size(PLVDelayAlpha.plvspctrm(ROI2_chan_idx, ROI1_chan_idx), 1), 1), row_label_pos, {dataDelay.channelInfo(ROI2_chan_idx).ROI}, ...
            'HorizontalAlignment', 'right');
        
        % Add column labels below the plot
        col_label_pos = 1:size(PLVDelayAlpha.plvspctrm(ROI2_chan_idx, ROI1_chan_idx), 2);
        text(col_label_pos, repmat(0.8, size(PLVDelayAlpha.plvspctrm(ROI2_chan_idx, ROI1_chan_idx), 2), 1), {dataDelay.channelInfo(ROI1_chan_idx).ROI}, ...
            'HorizontalAlignment', 'center');
        
        % Adjust the axes to make room for labels
        axis on;
        axis tight;
        colorbar
        title('mean PLV in 7-14 Hz');
        
        % save figure
        heatmap_filename = [ basedir pacienti(p).folder '\' subfolder '\'  'heatmap_PLV_7-14Hz_' ROI1 '-' ROI2 '_' condition ' trials ' setup.suffix '_' datestr(now,'YYYY-mm') '.fig'];
        savefig(heatmap_filename);
        close all
        
        % find ch pairs with PLV > 0.3 in between ROI ch pairs
        [chn2, chn1] = find(PLVDelayAlpha.plvspctrm(ROI2_chan_idx, ROI1_chan_idx) > 0.3); % PLV > 0.3
        chnpairs_highPLV = cell(numel(chn1), 2);
        for i = 1:numel(chn1)
            chnpairs_highPLV(i,1) = PLVDelayAlpha.label(ROI2_chan_idx(chn2(i))); % original labels of chan
            chnpairs_highPLV(i,2) = PLVDelayAlpha.label(ROI1_chan_idx(chn1(i)));
        end
        
        % plot separate chan pairs
        % (LOC - IPL) ROI1-ROI2
        fig_filename = [ basedir pacienti(p).folder '\' subfolder '\'  'PLV_' ROI1 '-' ROI2 '_' condition ' trials ' setup.suffix '_' datestr(now,'YYYY-mm')];        
        
        for i = 1:numel(chn1)
            figure(i); plot(PLVDelay.freq, squeeze(PLVDelay.plvspctrm(ROI2_chan_idx(chn2(i)),ROI1_chan_idx(chn1(i)),:)),'b', 'LineWidth', 0.9)
            hold on
            plot(PLVBaseline.freq, squeeze(PLVBaseline.plvspctrm(ROI2_chan_idx(chn2(i)),ROI1_chan_idx(chn1(i)),:)),'r', 'LineWidth', 0.9)    
            hold on
            plot(PLVEncod.freq, squeeze(PLVEncod.plvspctrm(ROI2_chan_idx(chn2(i)),ROI1_chan_idx(chn1(i)),:)),'m', 'LineWidth', 0.9)
            xlim([freq(1) freq(end)])
            xlabel('Frequency, Hz')
            ylabel('PLV')
            legend('delay 1.5 sec', 'baseline 1.5 sec', 'encoding 1.5 sec');
            title(['PLV between ' PLVDelay.label{ROI1_chan_idx(chn1(i)), 1} ' ' ROI1 ' and ' PLVDelay.label{ROI2_chan_idx(chn2(i)), 1} ' ' ROI2]);
            figi_filename = [fig_filename '_chnpair' num2str(i) '.fig'];
            savefig(figi_filename);
            close all
        end
        
        % save data 
        dataBS_filename = [ basedir pacienti(p).folder '\' subfolder '\'  'PLV_' ROI1 '-' ROI2 '_' condition ' trials baseline 1.5s ' setup.suffix '_' datestr(now,'YYYY-mm') '.mat'];
        save(dataBS_filename, 'dataBaseline', 'FreqBaseline', 'PLVBaseline');
        dataEncod_filename = [ basedir pacienti(p).folder '\' subfolder '\'  'PLV_' ROI1 '-' ROI2 '_' condition ' trials encoding 1.5s ' setup.suffix '_' datestr(now,'YYYY-mm') '.mat'];
        save(dataEncod_filename, 'dataEncod', 'FreqEncod', 'PLVEncod');
        dataDelay_filename = [ basedir pacienti(p).folder '\' subfolder '\'  'PLV_' ROI1 '-' ROI2 '_' condition ' trials delay 1.5s ' setup.suffix '_' datestr(now,'YYYY-mm') '.mat'];
        save(dataDelay_filename, 'dataDelay', 'FreqDelay', 'PLVDelay');
    end
end

end