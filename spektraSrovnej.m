fdr = 1;
channels = [90] ; % 23=PPA, p73 
els = [ 11 21 31 43 53 64 75 92 101]; %hranice elektrod pro p73 - konce
%ch = 23; 
freq = 50:10:150;
dataset1 = 1;
dataset2 = 5;
pocetvzorku = size(ALLEEG(dataset1).data,2);
T = 0: (1/ALLEEG(dataset1).srate) : (pocetvzorku-1)/ALLEEG(dataset1).srate; %cas zacatku a konce epochy

HHprumery = zeros( [numel(channels) pocetvzorku size(ALLEEG(dataset1).data,3)]); % prumerna hilbertova obalka: kanaly cas epochy
HHppp = zeros(numel(channels),1); %pocet frekvenci se signif rozdilem scene vs nonscene
for ch = channels
    fprintf(' ---- CHANNEL %i ----- \n',ch);
    fprintf('*** %s\n',ALLEEG(dataset1).setname);
    CondA = bipolarRef(ALLEEG(dataset1).data,els,0);  
    HHScene = spektra(ALLEEG(dataset1),CondA,ch,freq,true,true);

    fprintf('*** %s\n',ALLEEG(dataset2).setname);
    CondB = bipolarRef(ALLEEG(dataset2).data,els,0);  
    HHNonScene = spektra(ALLEEG(dataset2),CondB,ch,freq,true,true);

    W = WilcoxM(HHScene,HHNonScene,fdr);
    W = klouzaveokno(W,8,'mean');    
    
    pp = min(W,[],2); %signifikance vsech frekvenci
    ppp = sum(pp<0.05);
    HHppp(ch)=ppp;
    HHprumery(ch,:,:)= mean(HHScene,1);      
    fprintf('frekvence se signif rozdilem: %i\n',ppp);
    
    if ppp > 0 %obrazek kreslim, jen pokud by tam bylo neco videt
        figure('Name',['W map hilbert Scene vs NonScene, channel ' num2str(ch)]);
        %1. plot vsech frekvenci a jejich signif rozdilu
        subplot(2,1,1);
        imagesc(T,freq, 1-W,[ iff(fdr,0.95,0.99) 1]);%mapa, od p>0.05 bude modra barva 
        colorbar;
        axis xy;
        
        %2. plot prumeru vsech frekvenci
        subplot(2,1,2);
        meanHilbertA = squeeze(mean(mean(HHScene,1),3));
        plot(T,meanHilbertA,'Color','red','LineWidth',2);
        hold on;
        meanHilbertB = squeeze(mean(mean(HHNonScene,1),3));
        plot(T, meanHilbertB,'Color','blue','LineWidth',2);
        eventlatency = ALLEEG(dataset1).event(1,1).latency/ALLEEG(dataset1).srate;
        line([eventlatency eventlatency], [0 2.5],'LineWidth',2,'Color','black');
        title(['W map hilbert Scene vs NonScene, channel ' num2str(ch)]);
        %signifikantni rozdil prumeru
        Wm = WilcoxM(mean(HHScene,1),mean(HHNonScene,1),fdr);
        plot(T,Wm,'Color',[.5 .5 .5]);
        iWm = Wm<=0.05;
        plot(T(iWm),Wm(iWm),'m.');
        colorbar; %to je tam jen kvuli stejne sirce hodniho a dolniho grafu
    end
end