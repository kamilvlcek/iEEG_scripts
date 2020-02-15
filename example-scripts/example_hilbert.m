patient_file = '..\example-data\p073_PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11 FEX_CHilb.mat';

open(patient_file)
% file has no information about referencing
hilbert = CHilbertL(patient_file);

% Existing plot response for a single channel and multiple categories
hilbert.PlotResponseFreq(20, [0:3]);
% Doesn't work with the cells because I don't understand how the categories are being passed
hilbert.plotresponsefrequency(20, [0:3]);
hilbert.plotresponsefrequency([1:5], [0:3]);

% Newly added plot response which shoudl "average" given channels
hilbert.PlotResponseFreqMean([1:4], [0:3]);
hilbert.plotresponsefrequency([1:4], [0:3]);

% Wilcox against a baseline 
wBaseline = hilbert.wilcoxbaseline([-0.2 0], [0.01 0.8], [], {'Ovoce'});
freq1 = wBaseline(:, :, 1, hilbert.PsyData.CategoryNum({'Ovoce'}) + 1);
plotpintime(freq1, [0.01 0.8]);

% Wilcox category comparison
wCategory = hilbert.wilcoxcategories([{'Object'} {'Scene'}], [0.01 0.8]);
plotpintime(wCategory, [0.01 0.8])

wCategory = hilbert.wilcoxcategories([{'Ovoce'} {'Scene'}], [0.0 0.2]);
plotpintime(wCategory, [0.0 0.2])

wCategory = hilbert.wilcoxcategories([{'Ovoce'} {'Scene'}], [0.0 0.5]);
plotpintime(wCategory, [0.0 0.5])

wCategory = hilbert.wilcoxcategories([{'Ovoce'} {'Scene'}], [], 5);
plotpintime(wCategory, hilbert.epochtime(1:2))