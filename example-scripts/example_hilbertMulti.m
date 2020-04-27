%% Loading
patient_file = 'd:\eeg\motol\CHilbertMulti\PPA\CM PPA AnyResFE 50-150Hz refBipo -0.2-0.8 Ep2018-08 1-65 2020-03-31_CHMult.mat';

% origFile = open(patient_file);
% file has no information about referencing
CM = CHilbertMultiL(patient_file);
CM.PlotResponseCh;

%% 
% hilbertM.plotresponsefrequency(15);
CM.PlotResponseFreqMean();% 1:4,1:3

%% Wilcox second order
% Baseline
wp = CM.wilcoxaveragebaseline('categories', 1);
%'channels', 1:5, 'baseline', [-0.2 0], 'frequencies', 1:5, 'categories', {'Scene'}
plotpintime(wp(:, :, 2), [0.0 0.8], CM.Hfmean);