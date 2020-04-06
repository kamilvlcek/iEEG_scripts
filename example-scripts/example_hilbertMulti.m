%% Loading
patient_file = '..\example-data\CM PPA AnyResFE 50-150Hz refBipo -0.2-0.8 Ep2018-08 1-65 2020-03-31_CHilb.mat';

origFile = open(patient_file);
% file has no information about referencing
hilbertM = CHilbertMultiL(patient_file);

%% 
hilbertM.plotresponsefrequency(15);