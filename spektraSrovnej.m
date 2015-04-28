fdr = 1;
channels = 1:102; % 23=PPA, p73

%ch = 23; 
freq = 50:10:150;
dataset1 = 1;
dataset2 = 2;
pocetvzorku = size(ALLEEG(dataset1).data,2);
T = 0: (1/ALLEEG(dataset1).srate) : (pocetvzorku-1)/ALLEEG(dataset1).srate; %cas zacatku a konce epochy


HHprumery = zeros( [numel(channels) pocetvzorku size(ALLEEG(dataset1).data,3)]); % prumerna hilbertova obalka: kanaly cas epochy
HHppp = zeros(numel(channels),1); %pocet frekvenci se signif rozdilem scene vs nonscene
for ch = channels
    fprintf(' ---- CHANNEL %i ----- \n',ch);
    fprintf('*** %s\n',ALLEEG(dataset1).setname);
    HHScene = spektra(ALLEEG(dataset1),ch,freq,false,false);

    fprintf('*** %s\n',ALLEEG(dataset2).setname);
    HHNonScene = spektra(ALLEEG(dataset2),ch,freq,false,false);

    W = WilcoxM(HHScene,HHNonScene,fdr);
    W = klouzaveokno(W,8,'mean');    
    
    pp = min(W,[],2); %signifikance vsech frekvenci
    ppp = sum(pp<0.05);
    HHppp(ch)=ppp;
    HHprumery(ch,:,:)= mean(HHScene,1);
    fprintf('frekvence se signif rozdilem: %i\n',ppp);
    
    if ppp > 0 %obrazek kreslim, jen pokud by tam bylo neco videt
        figure('Name',['W map hilbert Scene vs NonScene, channel ' num2str(ch)]);
        subplot(2,1,1);
        imagesc(T,freq, 1-W,[ iff(fdr,0.95,0.99) 1]);%mapa, od p>0.05 bude modra barva 
        colorbar;
        axis xy;
        subplot(2,1,2);
        prumernyHilbert = squeeze(mean(mean(HHScene,1),3));
        plot(T,prumernyHilbert,'Color','red');
        hold on;
        plot(T, squeeze(mean(mean(HHNonScene,1),3)));
        eventlatency = ALLEEG(dataset1).event(1,1).latency/ALLEEG(dataset1).srate;
        line([eventlatency eventlatency], [0 2.5],'LineWidth',2,'Color','black');
        title(['W map hilbert Scene vs NonScene, channel ' num2str(ch)]);
    end
end