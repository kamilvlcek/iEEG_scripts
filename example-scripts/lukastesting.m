patient_file = '..\example-data\p073_PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11 FEX_CHilb.mat';

open(patient_file)
% Notes
% file has no information about referencing
% Stats are run immediately after initialization - not recommended
hilbert = CHilbert(patient_file);

% Existing plot response for a single channel and multiple categories
hilbert.PlotResponseFreq(20, [0:3]);

% Newly added plot response which shoudl "average" given channels
hilbert.PlotResponseFreq([1:4], [0:3]);

% Wilcox against a baseline 

wp = hilbert.wilcoxbaseline([-0.2 0], [0.01 0.8]);
