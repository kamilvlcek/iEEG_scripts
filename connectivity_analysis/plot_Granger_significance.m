function plot_Granger_significance(granger_data)
% visualization of channel pairs with significant Granger difference between 2 directions,
% plots all significant chan pairs and saves them to the patient's folder
% granger_data - filename of granger data which we want to plot

if(~exist('granger_data','var')) || isempty(granger_data), granger_data = 'Granger_VTC-IPL last 1.5s delay all_trials_2024-02.mat'; end
granger_folder = 'Granger_permut_stat';
setup = setup_memact(1); % setup for the delayed epochs
basedir = setup.basedir; % folder where the data of all patients stored
subfolder = setup.subfolder;
[pacienti] = pacienti_memact(); % struct with all patients in memact, some are marked not to analyze as they don't have ROIs we want

for p = 1:numel(pacienti)
    if pacienti(p).todo
        if isfile([basedir pacienti(p).folder '\' subfolder '\' granger_folder '\' granger_data])
            % load the PLV data with statistics
            load([basedir pacienti(p).folder '\' subfolder '\' granger_folder '\' granger_data]);
            
            % first find indexes of chan pairs with significant PLV difference
            significant_chanPairs = find(sum(Granger_signif_allPairs_clustcorr,2) ~= 0);
            
            iPairsGranger = 1:2:2*size(ROI_chanpairs, 1); % as we have 2 row of values for each ch pair
            
            % plot separate chan pairs with significant difference between 2 directions
            % (LOC - IPL) ROI1-ROI2
            fig_name_part = split(granger_data, '.mat');
            fig_filename = [basedir pacienti(p).folder '\' subfolder '\' granger_folder '\' fig_name_part{1}];
            
            % find the max granger value for setting ylim
            maxgranger1 = max(max(max(squeeze(grangerDelay.grangerspctrm(iPairsGranger(significant_chanPairs'), :)))));
            maxgranger2 = max(max(max(squeeze(grangerDelay.grangerspctrm(iPairsGranger(significant_chanPairs')+1, :)))));
            maxgranger = max(maxgranger1, maxgranger2);
            
            for i = 1:numel(significant_chanPairs)
                
                ipair = ROI_chanpairs(significant_chanPairs(i), :); % real indexes of channels in the data
                ipairGranger = iPairsGranger(significant_chanPairs(i)); % index of ch pair in granger spectrum
                significancePos = Granger_signif_allPairs_clustcorr(significant_chanPairs(i),:)>0;
                significanceNeg = Granger_signif_allPairs_clustcorr(significant_chanPairs(i),:)<0;
                
                figure(i); plot(grangerDelay.freq, grangerDelay.grangerspctrm(ipairGranger,:),'b', 'LineWidth', 1)
                hold on
                plot(grangerDelay.freq, grangerDelay.grangerspctrm(ipairGranger+1,:),'r', 'LineWidth', 1)
                hold on
                ylim([0 maxgranger*1.05])
                yLimits = ylim;
                y = yLimits(2) * 0.01;
                
                % plot significance
                ySignifPosValues = ones(1, numel(grangerDelay.freq)) * NaN; % Initialize with NaN
                ySignifPosValues(significancePos) = y; % Assign y values only at significant indices
                plot(grangerDelay.freq, ySignifPosValues, '.-b', 'LineWidth', 3, 'MarkerSize', 10);
                hold on
                
                ySignifNegValues = ones(1, numel(grangerDelay.freq)) * NaN; % Initialize with NaN
                ySignifNegValues(significanceNeg) = y; % Assign y values only at significant indices
                plot(grangerDelay.freq, ySignifNegValues, '.-r', 'LineWidth', 3, 'MarkerSize', 10);
                
                xlabel('Frequency, Hz')
                ylabel('Granger')
                chlabel1 = split(grangerDelay.labelcmb{ipairGranger,1}, '[');
                chlabel2 = split(grangerDelay.labelcmb{ipairGranger+1,1}, '[');
                
                legend(['from ' chlabel1{1} dataDelay.channelInfo(ipair(2)).ROI ' to ' chlabel2{1} dataDelay.channelInfo(ipair(1)).ROI], ['from ' chlabel2{1} dataDelay.channelInfo(ipair(1)).ROI ' to ' chlabel1{1} dataDelay.channelInfo(ipair(2)).ROI])
                title(['spectral Granger between ' chlabel1{1} ' ' dataDelay.channelInfo(ipair(2)).ROI ' and ' chlabel2{1} ' ' dataDelay.channelInfo(ipair(1)).ROI]);
                
                figi_filename = [fig_filename '_chnpair_' chlabel1{1} ' ' chlabel2{1} '_' num2str(i) '.fig'];
                savefig(figi_filename);
                close all
            end
            
        end
    end
end
end