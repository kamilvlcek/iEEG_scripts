% script to produce CM for memact across selected structures
% computes statistics for all channels, with fdr correction 
% then removes channels with non-significant responses relative to baseline

CB=CBrainPlot();

% select only structures in occipito-temporo-parietal region
% selectedStruct = {'Cun','TOTZ', 'LgG','IOG','MOG','SOG','LPHT','PHG','FG', 'FFG', 'FuG', 'ITG', 'MTG', 'STG', 'SPL', 'SMG', 'AnG', 'IPL'};
% including also TP, Hi and Ent and precun
selectedStruct = {'Cun','TOTZ', 'LgG','IOG','MOG','SOG','LPHT','PHG','FG', 'FFG', 'FuG', 'ITG', 'MTG', 'STG', 'TP','Hi','Ent','PCun', 'precun' 'SPL', 'SMG', 'AnG', 'IPL'};
exlcude = {'MFG', 'IFG', 'SFG'};
[PAC,ChNum] =  CB.StructFind({}, selectedStruct,'memact','b', exlcude); 

%% create CHilbertMulti and create extracts
% CHilbert_filename = 'MemAct CHilbert 50-150Hz -0.5-2.0 refBipo Ep2023-06_CHilb.mat'; % from which files to extract
% CHilbert_filename = 'MemAct CHilbert 2-150Hz -0.5-2.0 refBipo Ep2023-06_CHilb.mat';
% CHilbert_filename = 'MemAct CHilbert 8-13Hz -0.5-2.0 refBipo Ep2023-06_CHilb.mat';
% CHilbert_filename = 'MemAct CHilbert 15-30Hz -0.5-2.0 refBipo Ep2023-06_CHilb.mat';
CHilbert_filename = 'MemAct CMorlet 4-8HzM -0.8-2.0 refBipo Ep2023-06_CHilb.mat';

extracts_label = 'AnyResp_OTP_MTL'; % label for the extract files
% CMlabel = 'Memact_AnyResp_OTP_MTL_50-150Hz'; % choose label for CHilbertMulti file
% CMlabel = 'Memact_AnyResp_OTP_MTL_8-13Hz'; % choose label for CHilbertMulti file
% CMlabel = 'Memact_AnyResp_OTP_MTL_15-30Hz'; 
CMlabel = 'Memact_AnyResp_OTP_MTL_4-8HzM'; 

CM = CHilbertMulti();
filenames = CM.ExtractData(PAC,'memact',CHilbert_filename,extracts_label,0); % create extracts for these all channels, - overwrite the old ones 

%% import extract and compute statistics
CM.ImportExtract(filenames,CMlabel);  % import all these extracts to one CHilbertMulti object
setup = setup_memact(); % load setup for memact
CM.ResponseSearchMulti(0.1,setup.stat_kats,[],struct('chn',2, 'fdr',1)); % compute all contrasts, with fdr correction over all channels

%% store contrasts to markings
% alpha band
CM.SetStatActive(1);
CM.SelChannelStat ({1, 2, 3}, [1 2 3], 0);
% 7-15Hz: fgh = {{[immed_same]},{[immed_diff]},{[immed_diffXimmed_same]}}, (178  240   14)
% 8-13 Hz: fgh = {{[immed_same]},{[immed_diff]},{[immed_diffXimmed_same]}}, (155  233   10)
% 8-13 Hz with 7 patients fgh = {{[immed_same]},{[immed_diff]},{[immed_diffXimmed_same]}}, (235  350   28)
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();

%% gamma band
CM.SetStatActive(1); 
CM.SelChannelStat ({1, 2, 3}, [1 2 3], 0); % immed_same x bs, immed_diff x bs, immed_diff x immed_same; 
% fgh = (126  139   10), with 'chn', 2
% (115  113   20), with 'chn', 1
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();
% 7 patients: fgh = {{[immed_same]},{[immed_diff]},{[immed_diffXimmed_same]}}, (195  202   12)

%% store contrasts to markings
% beta band 28.06.2023
CM.SetStatActive(1); 
CM.SelChannelStat ({1, 2, 3}, [1 2 3], 0); % immed_same x bs, immed_diff x bs, immed_diff x immed_same;
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();
% 7 patients: fgh = {{[immed_same]},{[immed_diff]},{[immed_diffXimmed_same]}}, (160  285   23)

%% theta band 28.06.2023
CM.SetStatActive(1); 
CM.SelChannelStat ({1, 2, 3}, [1 2 3], 0); % immed_same x bs, immed_diff x bs, immed_diff x immed_same;
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();
% 7 patients: fgh = {{[immed_same]},{[immed_diff]},{[immed_diffXimmed_same]}}, (332  334   13)

%% save the file
% CHilbertMulti_filename = 'CM Memact CHilbert 50-150Hz AllChan_OTP -0.5-2.0 refBipo Ep2023-06_CiEEG.mat'; % the name of the resulting CM file
% CHilbertMulti_filename = 'CM Memact CHilbert 50-150Hz AllChan_OTP_stat_chn-1 -0.5-2.0 refBipo Ep2023-06_CiEEG.mat';
% CHilbertMulti_filename = 'CM Memact CHilbert 2-150Hz AllChan_OTP -0.5-2.0 refBipo Ep2023-06_CiEEG.mat';
CHilbertMulti_filename = 'CM Memact CHilbert 8-13Hz AllChan_OTP -0.5-2.0 refBipo Ep2023-06_CiEEG.mat';
CM.Save(['d:\EEG\motol\CHilbertMulti\MemAct\' CHilbertMulti_filename]); % save the file

%% remove non-significant channels
markings = [1 2]; % immed_same x bs, immed_diff x bs
iCh = any(CM.plotRCh.selCh(:,markings),2); % index of channel with first or second marking (f or g)
CM.RemoveChannels(find(~iCh)); % remove all channels without this first or second marking (channel which are not significant vs baseline)

%% save file 
% CHilbertMulti_filename = 'CM Memact CHilbert 50-150Hz AnyResp_OTP -0.5-2.0 refBipo Ep2023-06_CiEEG.mat'; % the name of the resulting CM file
% CHilbertMulti_filename = 'CM Memact CHilbert 8-13Hz AnyResp_OTP_MTL -0.5-2.0 refBipo Ep2023-06_CiEEG.mat';
% CHilbertMulti_filename = 'CM Memact CHilbert 50-150Hz AnyResp_OTP_MTL -0.5-2.0 refBipo Ep2023-06_CiEEG.mat';
% CHilbertMulti_filename = 'CM Memact CHilbert 15-30Hz AnyResp_OTP_MTL -0.5-2.0 refBipo Ep2023-06_CiEEG.mat';
CHilbertMulti_filename = 'CM Memact CMorlet 4-8HzM AnyResp_OTP_MTL -0.8-2.0 refBipo Ep2023-06_CiEEG.mat';
% CM.Save(['d:\EEG\motol\CHilbertMulti\MemAct\' CHilbertMulti_filename]); % save the file
CM.Save(['d:\eeg\motol\MorletMulti\MemAct\' CHilbertMulti_filename]);

%% recompute statistics only across active channels
setup = setup_memact(); % load setup for memact
CM.ResponseSearchMulti(0.1,setup.stat_kats,[],struct('chn',2, 'fdr',1));

%% gamma band - only active channels
CM.SetStatActive(1);
CM.SelChannelStat ({1, 2, 3}, [1 2 3], 0); % immed_same x bs, immed_diff x bs, immed_diff x immed_same; 
% 6 patients, fgh = (140  142   15)
% 7 patients: (212  221   23)
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();
CM.SelChannelStat ({{2 '&~1'}}, [4], 1, 1); % only positive change
% j = {{[immed_diffANimmed_same]}}, (24) % only immed_diff vs bs
CM.SelChannelStat ({{1 '&~2'}}, [5], 1, 1);
% k = {{[immed_sameANimmed_diff]}}, (9) % only immed_same vs bs
CM.SelChannelStat ({{1 '&2' '&~3'}}, [6], 1, 1);
% l = {{[immed_sameAimmed_diff]}}, (89) % general chan

%% alpha band - only active channels
CM.SetStatActive(1);
CM.SelChannelStat ({1, 2, 3}, [1 2 3], 0); % immed_same x bs, immed_diff x bs, immed_diff x immed_same; 
% 8-13 Hz 7 patients: 244  357   36
CM.SelChannelStat ({{2 '&~1'}}, [4], 1, 1); % 8-13 Hz
% j = {{[immed_diffANimmed_same]}}, (77)
CM.SelChannelStat ({{2 '&~1'}}, [5], 1, -1);
% k = {{[immed_diffANimmed_same]}}, (41)
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();

CM.SetSelChActive(2);
CM.SelChannelStat ({1, 2}, [1 2], 0);
% fg = {{[immed_same]},{[immed_diff]}}, (244  357)
CM.SelChannelStat ({3}, [3], 1, 1); % 
% h = {{[immed_diffXimmed_same]}}, (8) immed_diff > immed_same
CM.SelChannelStat ({3}, [4], 1, -1);
% j = {{[immed_diffXimmed_same]}}, (28) immed_diff < immed_same

%% beta band - only active channels
CM.SetStatActive(1);
CM.SelChannelStat ({1, 2, 3}, [1 2 3], 0); % immed_same x bs, immed_diff x bs, immed_diff x immed_same; 
% 7 patients: (182  292   33)
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();

%% theta band - only active channels
CM.SetStatActive(1);
CM.SelChannelStat ({1, 2}, [1 2], 0); % immed_same x bs, immed_diff x bs, immed_diff x immed_same; 
% 7 patients: (342  348)
CM.SelChannelStat ({3}, [3], 1, 1); % 
%h = {{[immed_diffXimmed_same]}}, (11) immed_diff > immed_same
CM.SelChannelStat ({3}, [4], 1, -1);
%j = {{[immed_diffXimmed_same]}}, (9) immed_diff < immed_same
CM.SetSelChActive(2);
CM.SetSelChActive(1); % set first set of markings
CM.PlotResponseCh();

%% import rejected channels from previous CM
CM.RejectChannelsImport('d:\eeg\motol\CHIlbertMulti\MemAct\CM Memact CHilbert 50-150Hz AnyResp_OTP -0.5-2.0 refBipo Ep2023-06_CiEEG.mat');
