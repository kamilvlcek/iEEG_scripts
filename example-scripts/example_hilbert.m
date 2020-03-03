patient_file = '..\example-data\p073_PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11 FEX_CHilb.mat';

open(patient_file)
% file has no information about referencing
hilbert = CHilbertL(patient_file);

% Existing plot response for a single channel and multiple categories
hilbert.PlotResponseFreq(20, 0:3);
% Doesn't work with the cells because I don't understand how the categories are being passed
hilbert.plotresponsefrequency(20, 0:3);
hilbert.plotresponsefrequency(1:5, 0:3);

% Newly added plot response which shoudl "average" given channels
hilbert.PlotResponseFreqMean(1:4, 0:3);
hilbert.plotresponsefrequency(1:4, 0:3);

%% Wilcox 
%against a baseline for all frequencies
wBaseline = hilbert.wilcoxbaseline([-0.2 0], [0.01 0.8], [], {'Ovoce'});
% Wilcox agains baseline returs a 4D matrix with
% time x channel x frequency x category p values.
% the following plots the p values for time x channel for the
% first frequency (52 Hz) and first category ('Ovoce')
plotpintime(wBaseline(:, :, 1, 1), [0.01 0.8]);
% This plots the time x Frequencies for the 18the channel for category ('Ovoce')
plotpintime(squeeze(wBaseline(:, 18, :, 1)), [0.01 0.8], hilbert.Hfmean);

% Comparing for Scene
wBaseline = hilbert.wilcoxbaseline([-0.2 0], [0.01 0.8], [], {'Object'});
plotpintime(squeeze(wBaseline(:, 1, :, 4)), [0.01 0.8], hilbert.Hfmean);

% Wilcox category comparison
wCategory = hilbert.wilcoxcategories([{'Ovoce'} {'Scene'}], [0.0 0.8]);
% PLots P value for difference between Ovoce and Scene for 18th channel for
% all frequencies
plotpintime(squeeze(wCategory(:, 18, :)), [0.0 0.8], hilbert.Hfmean);

wCategory = hilbert.wilcoxcategories([{'Ovoce'} {'Scene'}], [0.0 0.5]);
plotpintime(squeeze(wCategory(:, 18, :)), [0.0 0.8])

%% Testing
[~, ~, ~, iEpochs] = hilbert.CategoryData(1);
per = permute(hilbert.HFreqEpochs, [1 3 4 2]);
sel = per(:,:,iEpochs);
res = reshape(sel, [64 20 51 20]);
% Checking how it was reshaped
all(sel(:,1,1) == res(:, 1,1,1))
all(sel(:,1,2) == res(:, 1,2,1))
all(sel(:,1,52) == res(:, 1,1,2))
% K, this was reshaped as I expected - e.g. sel is a vector of 51 events
% for each channel (20) and reshaping it keeps the first two dimensions and
% then splits it into 51 x 20