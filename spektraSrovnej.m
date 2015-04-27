fdr = 1;
ch = 23;
freq = 50:10:150;
T = 0:0.1:size(EEG.data,2)/EEG.srate; %cas zacatku a konce epochy
fprintf('*** %s\n',ALLEEG(1).setname);
HHScene = spektra(ALLEEG(1),ch,freq);

fprintf('*** %s\n',ALLEEG(2).setname);
HHNonScene = spektra(ALLEEG(2),ch,freq);

W = WilcoxM(HHScene,HHNonScene,fdr);
W = klouzaveokno(W,25);
figure('Name','W map hilbert Scene vs NonScene');
imagesc(T,freq, 1-W,[ iff(fdr,0.95,0.99) 1]);%mapa, od p>0.05 bude modra barva 
colorbar;
axis xy;

