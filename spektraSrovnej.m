fdr = 1;
W = WilcoxM(allHHScene,allHHNonScene,fdr);
figure('Name','W map hilbert Scene vs NonScene');
imagesc(T,F, 1-W,[ iff(fdr,0.95,0.99) 1]);%mapa, od p>0.05 bude modra barva 
colorbar;
axis xy;