%% Loading
patient_file = 'd:\EEG\motol\CHilbertMulti\PPA\CM PPA AnyResp 50-150Hz refBipo -0.2-0.8 Ep2018-08 2020-07-13_CHMult.mat'; %6.8.2020

% origFile = open(patient_file);
% file has no information about referencing
CM = CHilbertMultiL(patient_file);
CM.PlotResponseCh;

%% 
% hilbertM.plotresponsefrequency(15);
CM.PlotResponseFreqMean([1 3], [0 .5]); %6.8.2020