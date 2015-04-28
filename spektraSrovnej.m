fdr = 1;
ch = 23; % 23=PPA, p73
freq = 50:10:150;
T = 0:0.1:size(EEG.data,2)/EEG.srate; %cas zacatku a konce epochy
dataset1 = 1;
fprintf('*** %s\n',ALLEEG(dataset1).setname);
HHScene = spektra(ALLEEG(dataset1),ch,freq);

dataset2 = 2;
fprintf('*** %s\n',ALLEEG(dataset2).setname);
HHNonScene = spektra(ALLEEG(dataset2),ch,freq);

W = WilcoxM(HHScene,HHNonScene,fdr);
W = klouzaveokno(W,8,'mean');
figure('Name','W map hilbert Scene vs NonScene');
imagesc(T,freq, 1-W,[ iff(fdr,0.95,0.99) 1]);%mapa, od p>0.05 bude modra barva 
colorbar;
axis xy;
pp = min(W,[],2); %signifikance vsech frekvenci
fprintf('frekvenci se signif rozdilem: %i\n',sum(pp<0.05));

