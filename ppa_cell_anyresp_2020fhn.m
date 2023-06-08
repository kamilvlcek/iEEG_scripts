%11.7.2020 All channels 2707
%script used to produce CM PPA AnyResp 50-150Hz refBipo -0.2-0.8 Ep2018-08 2020-07-13_CHMult.mat 
%used in the final version of vlcek et al 2020 fhn.pdf
%computes statistics for all channels, with fdr correction 
%then removes channels with non significant response relative to baseline for Scene or Object

CB=CBrainPlot();
pactodo = {'VT53','VT55','VT56','VT57','VT58','VT59'};
[PAC,ChNum] =  CB.StructFind([],[],'ppa','b',[],pactodo); %find all bipolar channels across all patients for the PPA test 
%you can edit the PAC file at this stage - e.g. remove some patients
%PAC(1:2707)=[];

%% create CHilbertMulti and create extracts
%CHilbert_filename = 'PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11_CHilb.mat'; %from which files to extract
CHilbert_filename = 'PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2023-02_CiEEG.mat'; %from which files to extract

extracts_label = 'AnyResp'; %label for the extract files
CMlabel = 'PPA_AnyResp_50-150Hz'; %choose label for CHilbertMulti file

CM = CHilbertMulti();
filenames = CM.ExtractData(PAC,'ppa',CHilbert_filename,extracts_label,1); % create extracts for these all channels, - overwrite the old ones. 
% check - all extracts created?

%% import extract and compute statistics
CM.ImportExtract(filenames,CMlabel);  %import all these extracts to one CHilbertMulti object
setup = setup_ppa(); %load setup for PPA
CM.ResponseSearchMulti(0.1,setup.stat_kats,[],struct('chn',2, 'fdr',1)); %compute all contrasts, with fdr correction over all channels.

%% store contrasts to markings
%default first set of markings
CM.SetStatActive(1); 
%set contrast 1 with categories: [2  3  1], opakovani: {}, 2=Face x 3=Object x 1=Scene
CM.SelChannelStat({1, 2, 3, 6 },[3 2 1 4],0,1); %kategory,marks,add,signum - {Face, Object,Scene,SceneXObject},[h g f j ]
% fghj = {[Scene],[Object],[Face],[SceneXObject]}, (270  359  302   70)
CM.SetStatActive(2);
%set contrast 2 with categories: [1  2  3], opakovani: {}, 1=Scene x 2=Face x 3=Object
CM.SelChannelStat ({ 5}, [5], 1,1);  %kategory,marks,add,signu - {ObjectXScene} , [ k]
% k = {[ObjectXScene]}, (105)

CM.SetSelChActive(2); %set second set of markings
CM.SetStatActive(1);
%set contrast 1 with categories: [2  3  1], opakovani: {}, 2=Face x 3=Object x 1=Scene
CM.SelChannelStat ({1, 2, 3, {2 '&~3'} , {3 '&~2'} , {'&',2,3}}, [3 2 1 4 5 6], 0,1);
%fghjkl = {[Scene],[Object],[Face],[ObjectANScene],[SceneANObject],[ObjectAScene]}, (270F  359O  302S  146   89  213)
CM.SetSelChActive(1); %set first set of markings

%% save the file
%CHilbertMulti_filename = 'CM PPA AnyResp 50-150Hz refBipo -0.2-0.8 Ep2018-08 2020-07-13'; %the name of the resulting CM file
CHilbertMulti_filename = 'CM PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2023-02_CiEEG.mat'; %the name of the resulting CM file
%CHilbertMulti_previous_filename = 'CM PPA AnyResp 50-150Hz refBipo -0.2-0.8 Ep2018-08 2020-01-03_CiEEG.mat'; %CM file from which to import
CM.Save(['d:\EEG\motol\CHilbertMulti\' CHilbertMulti_filename]); %save the file

%% plot channel
CM.PlotResponseCh; % check if everything is OK

%% remove non significant channels
CHilbertMulti_previous_filename = []; %CM file from which to import
markings = [1 2]; %markins to keep (first set of markings - [Scene],[Object],[Face],[SceneXObject],[ObjectXScene])
iCh = any(CM.plotRCh.selCh(:,markings),2); %index of channel with first or second marking (f or g= Scene or Object)

CM.RemoveChannels(find(~iCh)); %#ok<FNDSB> %remove all channels without this first or second marking

if ~isempty(CHilbertMulti_previous_filename)
    %import list of rejected channels from previous file                         
    CM.RejectChannelsImport(['d:\EEG\motol\CHilbertMulti\' CHilbertMulti_previous_filename]);

    %set third set of markings - to be able to select just new files, not present in CHilbertMulti_previous_filename
    CM.SetSelChActive(3);                                   
    %stores comparisons of channels with a previous file f=same, g=only in the new file
    CM.CompareChannels('d:\EEG\motol\CHilbertMulti\CM PPA AnyResp 50-150Hz refBipo -0.2-0.8 Ep2018-08 2020-01-03_CiEEG.mat');
    CM.CH.FilterChannels({},{},'g'); %show only the new channels
end

%% save file 
CM.Save(); %save the file again