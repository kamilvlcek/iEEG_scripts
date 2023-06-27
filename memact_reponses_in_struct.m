% script to produce CM for memact across selected structures
% computes statistics for all channels, with fdr correction 
% then removes channels with non-significant responses relative to baseline
% Sofiia 23.6.2023 

CB=CBrainPlot();

% select only structures in occipito-temporo-parietal region
selectedStruct = {'Cun','TOTZ', 'LgG','IOG','MOG','SOG','LPHT','PHG','FG', 'FFG', 'FuG', 'ITG', 'MTG', 'STG', 'SPL', 'SMG', 'AnG', 'IPL'};
exlcude = {'MFG', 'IFG', 'SFG'};
[PAC,ChNum] =  CB.StructFind({}, selectedStruct,'memact','b', exlcude); 


%% create CHilbertMulti and create extracts
CHilbert_filename = 'MemAct CHilbert 50-150Hz -0.5-2.0 refBipo Ep2023-06_CHilb.mat'; % from which files to extract

extracts_label = 'AnyResp_OTP'; % label for the extract files
CMlabel = 'Memact_AnyResp_OTP_50-150Hz'; % choose label for CHilbertMulti file

CM = CHilbertMulti();
filenames = CM.ExtractData(PAC,'memact',CHilbert_filename,extracts_label,1); % create extracts for these all channels, - overwrite the old ones 

%% import extract and compute statistics
CM.ImportExtract(filenames,CMlabel);  % import all these extracts to one CHilbertMulti object
setup = setup_memact(); % load setup for memact
CM.ResponseSearchMulti(0.1,setup.stat_kats,[],struct('chn',2, 'fdr',1)); % compute all contrasts, with fdr correction over all channels

%% store contrasts to markings
CM.SetStatActive(1); 
CM.SelChannelStat ({1, 2, 3}, [1 2 3], 0); % immed_same x bs, immed_diff x bs, immed_diff x immed_same
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();

%% save the file
CHilbertMulti_filename = 'CM Memact CHilbert 50-150Hz AllChan_OTP -0.5-2.0 refBipo Ep2023-06_CiEEG.mat'; % the name of the resulting CM file
CM.Save(['d:\EEG\motol\CHilbertMulti\' CHilbertMulti_filename]); % save the file

%% remove non significant channels
markings = [1 2]; % immed_same x bs, immed_diff x bs
iCh = any(CM.plotRCh.selCh(:,markings),2); % index of channel with first or second marking (f or g)
CM.RemoveChannels(find(~iCh)); % remove all channels without this first or second marking (channel which are not significant vs baseline)

%% save file 
CHilbertMulti_filename = 'CM Memact CHilbert 50-150Hz AnyResp_OTP_v2 -0.5-2.0 refBipo Ep2023-06_CiEEG.mat'; % the name of the resulting CM file
CM.Save(['d:\EEG\motol\CHilbertMulti\' CHilbertMulti_filename]); % save the file