%24.4.2015 - zkousim porovnat vysledky funkce spectrogram ze spektrem z EEGlabu
% viz d:\EEG\motol\pacienti\p68 Daenemark\p68 spectra.cdr
wsize_sec = 0.06;
ch = 23;
freq = 10:10:150;

wsize = floor(wsize_sec * EEG.srate);
woverlap = 2; %pro jine hodnoty nechapu pocet hodnot, ktere vraci spectrogram
allPP = zeros( [ size(EEG.data,3) size(freq,2)  floor(size(EEG.data,2)/wsize*woverlap-1) ] ); %pocet hodnot spektrogramu zkousim odhadnout 
allHH = zeros( [ size(EEG.data,3) size(freq,2)-1  size(EEG.data,2)] ); %epochy frekvence cas

%vypocet prumerne sily ve frekvencnich pasmech
ffPP = zeros( [size(EEG.data,3) size(freq,2)  ]); %prumerny spectrogram pro kazdou epochu a frekvenci
ffHH = zeros( [size(EEG.data,3) size(freq,2)-1]); %prumerny hilbert za kazdou epochu a frekvenci
for epocha = 1: size(EEG.data,3);
    %spectrogram
    [~,~,~,P]=spectrogram(double(EEG.data(ch,:,epocha)),wsize,wsize/woverlap,freq,EEG.srate,'yaxis');
    PP = 10*log10(abs(P));
    %hilbert
    for f = 1:numel(freq)
        ffPP(epocha,f)=mean(PP(f,:));
        if f < numel(freq) %hilbert je vzdy od do, takze pro posledni neni definovany
            ffHH(epocha,f)=mean(hilbertJirka(double(EEG.data(ch,:,epocha)),freq(f), freq(f+1),EEG.srate));
        end
    end
end
nPP = mean(ffPP,1); %prumer spectrogram pro kazdou frekvenci pres vsechny epochy
nHH = mean(ffHH,1); %prumer hilbert pro kazdou frekvenci pres vsechny epochy


fprintf('ep ');
for epocha = 1: size(EEG.data,3);
    %spectrogram
    [S,F,T,P]=spectrogram(double(EEG.data(ch,:,epocha)),wsize,wsize/woverlap,freq,EEG.srate,'yaxis');
    PP = 10*log10(abs(P));
    for f = 1:numel(freq)
        PP(f,:) = PP(f,:) ./ nPP(f); 
    end
    allPP(epocha,:,:)=PP; %frekvence=radky, cas=sloupce
    
    %hilbert
    for f = 1:numel(freq)-1
        H = hilbertJirka(double(EEG.data(ch,:,epocha)),freq(f), freq(f+1),EEG.srate);
        H = H ./ nHH(f);
        allHH(epocha,f,:)=H;
    end
    
    fprintf('%i ',epocha);
end
meanPP = squeeze(mean(allPP,1));
figure('Name','prumerny spektrogram');
imagesc(T,F,meanPP);
axis xy;
hold on;
line([eventlatency eventlatency], [freq(1) freq(end)],'LineWidth',2,'Color','black');
colorbar;
caxis([-0,2.5]);

meanHH = squeeze(mean(allHH,1));
figure('Name','prumerny hilbert');
imagesc(T,F,meanHH);
axis xy;
hold on;
eventlatency = EEG.event(1,1).latency/EEG.srate;
line([eventlatency eventlatency], [freq(1) freq(end)],'LineWidth',2,'Color','black');
colorbar;
caxis([-0,2.5]);

disp('hotovo');
