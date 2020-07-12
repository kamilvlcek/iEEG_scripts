%% Loading
patient_file = '..\example-data\p073_PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11 FEX_CHilb.mat';
hilbert = CHilbertL(patient_file);

%% Second order wilcox for category response vs baseline comparisons
channels = 1:5;
frequencies = 1:50;
categorySelect = 1;
wp = hilbert.wilcoxaveragebaseline('channels', 1:5, 'baseline', [-0.2 0], 'categories', categorySelect - 1);
categoryData = wp(:, :, categorySelect);

title = ['Significant difference from baseline in ', hilbert.PsyData.CategoryName(categorySelect - 1)];
figure('Name', title);
plt = imagesc([0.0 0.8], hilbert.Hfmean, (1 - categoryData)', [0.95 1]);
set(plt, 'AlphaData', ~isnan(categoryData'));
axis xy;
xlabel('time');
colorbar;

hilbert.plotresponsefrequency(1:5, 0);