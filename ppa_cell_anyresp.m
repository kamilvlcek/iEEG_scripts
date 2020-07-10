%find responding channels 
CB = CBrainPlot;
filename = 'PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11_CHilb.mat';
CB.IntervalyResp('ppa',[0 0.8],filename,1,1);

%% fdr correction
CB.FdrVALS(1,'fdr1');
% number of removed/left channels in 7 categories (across all patients):
% Face:	 78/207,
% Object:	 74/285,
% Scene:	 62/251,
% ObjectXFace:	 61/117,
% SceneXFace:	 43/98,
% SceneXObject:	 36/55,
% AllEl:	 0/2707,

%% create extracts and import them
% only objects and scenes
PAC = CB.MergePAC(CB.PAC{1,2},CB.PAC{1,3});  %2=Face x 3=Object x 1=Scene
CM = CHilbertMulti();
filenames = CM.ExtractData(PAC,'ppa','PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11_CHilb.mat','AnyResp',1); %overwrite
% extrahovano 26 souboru z 26

CMlabel = ['PPA_AnyResp07_50-150Hz']; %choose label for CHilbertMulti file
CM.ImportExtract(filenames,CMlabel);
%% compute statistics
setup = setup_ppa();
CM.ResponseSearchMulti(0.1,setup.stat_kats,[],struct('chn',2, 'fdr',1)); 
CM.PlotResponseCh;
%% store statistics it to markings
CM.SetStatActive(1); 
%nastavena statistika 1 s kats: [2  3  1], opakovani: {}, 2=Face x 3=Object x 1=Scene
CM.SelChannelStat({1, 2, 3, 6 },[3 2 1 4],0,1); %kategorie,marks,add,signum - 
% fghj = {[Scene],[Object],[Face],[SceneXObject]}, (205  310  295   88)
CM.SetStatActive(2);
%nastavena statistika 2 s kats: [1  2  3], opakovani: {}, 1=Scene x 2=Face x 3=Object
CM.SelChannelStat ({ 5}, [5], 1,1);  %kategorie,marks,add,signu
% k = {[ObjectXScene]}, (142)
CM.SetSelChActive(2);
CM.SetStatActive(1);
%nastavena statistika 1 s kats: [2  3  1], opakovani: {}, 2=Face x 3=Object x 1=Scene
CM.SelChannelStat ({1, 2, 3, {2 '&~3'} , {3 '&~2'} , {'&',2,3}}, [3 2 1 4 5 6], 0,1);
%fghjkl = {[Scene],[Object],[Face],[ObjectANScene],[SceneANObject],[ObjectAScene]}, (205face  310 Object  295Scene   64   49  246)
CM.SetSelChActive(1);

%% save file
CM.Save('d:\eeg\motol\CHilbertMulti\PPA\CM PPA AnyResp 50-150Hz refBipo -0.2-0.8 Ep2018-08 2020-07-07');