testname = 'ppa'; %menrot, ppa, aedist
reference = 'refBipo'; %'refBipo', refHead, refEle
IntervalyRespSignum = 1;
intervals = [(0:.1:.7)' , (0.1:.1:.8)' ]; 
filename = 'PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2018-08_CHilb.mat';
kontrast = 5; %5; %ktera ze statistik se ma zobrazit v mozku
%% 1. HLEDANI SIGNIFIKANCI
CB = CBrainPlot; %vytvorim tridu
[ pacienti, setup  ] = pacienti_setup_load( testname );
%stat = setup.stat_kats{kontrast};  
CB.IntervalyResp(testname,intervals,filename,kontrast,IntervalyRespSignum); %ziskam signif rozdily pro kategorie a mezi kategoriemi pro vsechny pacienty   

%% 2. Zobrazeni do mozku
CB.PlotBrain3DConfig(struct('Names',1,'NoNames',0,'signum',1,'overwrite',1));
% TODO - presunout PlotBrain3DConfig do PlotBrain3D
% TODO - PlotBrain3D - udelat vypis vsech kanalu a v kterych kategoriich a intervalech jsou uvedene 
CB.PlotBrain3D(); %udela vsechny kategorie
