patient_file = '..\example-data\p073_PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11 FEX_CHilb.mat';
open(patient_file)
% Notes
% file has no information about referencing
% 
hilbert = CHilbert(filename); % filename kterehokoli v tech dvou
% Stats are run immediately after initialization - not recommended
% 

hilbert.PlotResponseCh(20)