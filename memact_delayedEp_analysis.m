% script to produce CM for memact test in delayed epochs across selected structures
% first searches for files before delay, than searches for files after delay and then joins them and normalize epochs
% computes statistics for all channels, with fdr correction 
% then removes channels with non-significant responses relative to baseline and recomputes statistic with FDR again between cats
% Sofiia 26.07.2023 

CB=CBrainPlot();

% select structures only in parietal and temporal cortex
selectedStruct = {'TOTZ', 'LgG', 'LPHT','PHG','FG', 'FFG', 'FuG', 'ITG', 'MTG', 'STG', 'TP','Hi','Ent','PCun', 'precun', 'SPL', 'SMG', 'AnG', 'IPL'};
exlcude = {'MFG', 'IFG', 'SFG'};
[PAC,ChNum] =  CB.StructFind({}, selectedStruct,'memact','b', exlcude); 
% found 663 channels in 7 pacients

%% create CHilbertMulti before delay and create extracts
CHilbert_filenamebd = 'MemAct CHilbert 50-150Hz -0.5-4.0 refBipo Ep2023-07bdel_CHilb.mat';
extracts_label = 'AllChan_PT'; % label for the extract files
CMlabelbd = 'Memact_AllChan_PT_50-150Hz bdel'; 
CMbd = CHilbertMulti();
filenamesbd = CMbd.ExtractData(PAC,'memact',CHilbert_filenamebd,extracts_label,0); % create extracts for these all channels, - doesn't overwrite the old ones 

%% import extracts for CM before delay
CMbd.ImportExtract(filenamesbd,CMlabelbd);  % import all these extracts to one CHilbertMulti object

%% create CHilbertMulti after delay and create extracts
CHilbert_filenamead = 'MemAct CHilbert 50-150Hz -1.9-2.0 refBipo Ep2023-07adel_CHilb.mat';
extracts_label = 'AllChan_PT'; % label for the extract files
CMlabelad = 'Memact_AllChan_PT_50-150Hz bdel'; 
CMad = CHilbertMulti();
filenamesad = CMad.ExtractData(PAC,'memact',CHilbert_filenamead,extracts_label,0); % create extracts for these all channels, - doesn't overwrite the old ones 

%% import extracts for CM after delay
CMad.ImportExtract(filenamesad,CMlabelad);  % import all these extracts to one CHilbertMulti object

%% join two objects and normalize epochs, compute statistics
CMbd.AppendData(CMad);
CMbd.NormalizeEpochs([-0.5 -0.1]); % normalize by mean baseline -0.5 -0.1 before encoding phase

setup = setup_memact(1); % load setup for memact with delayed epochs
CMbd.ResponseSearchMulti(0.1,setup.stat_kats,[],struct('chn',2, 'fdr',1)); % compute all contrasts, with fdr correction over all channels
CMbd.PlotResponseCh();

%% save the file with all chan in parieto-temporal cortex
CHilbertMulti_filename = 'CM Memact CHilbert 50-150Hz AllChan_PT -0.5-7.9 refBipo Ep2023-07delayed_CiEEG.mat';
CMbd.Save(['d:\EEG\motol\CHilbertMulti\MemAct\' CHilbertMulti_filename]); % save the file

%% store contrasts to markings
CMbd.SetStatActive(1);
CMbd.SelChannelStat ({1, 2, 3}, [1 2 3], 0);
CMbd.SetSelChActive(2);
CMbd.SetSelChActive(1); % set first set of markings
CMbd.PlotResponseCh();

%% remove non-significant channels
markings = [1 2]; % del_same x bs, del_diff x bs
iCh = any(CMbd.plotRCh.selCh(:,markings),2); % index of channel with first or second marking (f or g)
CMbd.RemoveChannels(find(~iCh)); % remove all channels without this first or second marking (channel which are not significant vs baseline)
% 288 channels removed

%% recompute statistics only across active channels
CMbd.ResponseSearchMulti(0.1,setup.stat_kats,[],struct('chn',2, 'fdr',1));

%% contrasts
CMbd.SetStatActive(1);
CMbd.SelChannelStat ({1, 2, 3}, [1 2 3], 0); % del_same x bs, del_diff x bs, del_same x del_diff; 
% 7 patients: (311 287 0)
CMbd.SetSelChActive(2);
CMbd.SetSelChActive(1); % set first set of markings
CMbd.PlotResponseCh();

%% save the file with active chan 
CHilbertMulti_filename = 'CM Memact CHilbert 50-150Hz AnyResp_PT -0.5-7.9 refBipo Ep2023-07delayed_CiEEG.mat';
CMbd.Save(['d:\EEG\motol\CHilbertMulti\MemAct\' CHilbertMulti_filename]); % save the file
