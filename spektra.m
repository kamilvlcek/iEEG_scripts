function [allHH,T,F]=spektra(EEG, ch,freq, Pyes)
%24.4.2015 - zkousim porovnat vysledky funkce spectrogram ze spektrem z EEGlabu
% viz d:\EEG\motol\pacienti\p68 Daenemark\p68 spectra.cdr
wsize_sec = 0.06;
%ch = 23;
%freq = 10:10:150;
if nargin < 4 
    Pyes = 0;
end
eventlatency = EEG.event(1,1).latency/EEG.srate;
allHH = zeros( [  size(freq,2)-1  size(EEG.data,2) size(EEG.data,3)] ); % frekvence cas epochy

%vypocet prumerne sily ve frekvencnich pasmech
fprintf('prumery za frekvence:');
ffHH = zeros( [ size(freq,2)-1 size(EEG.data,3)]); %prumerny hilbert za kazdou frekvenci a epochu

for epocha = 1: size(EEG.data,3);
    %hilbert
    for f = 1:numel(freq)
        if f < numel(freq) %hilbert je vzdy od do, takze pro posledni neni definovany
            ffHH(f,epocha)=mean(hilbertJirka(double(EEG.data(ch,:,epocha)),freq(f), freq(f+1),EEG.srate));
        end
    end
end

nHH = mean(ffHH,2); %prumer hilbert pro kazdou frekvenci pres vsechny epochy
fprintf('done\n');

% prubehy frekvenci v case
fprintf('hilbert epocha: ');
for epocha = 1: size(EEG.data,3)
    %hilbert
    for f = 1:numel(freq)-1
        H = hilbertJirka(double(EEG.data(ch,:,epocha)),freq(f), freq(f+1),EEG.srate);
        H = H ./ nHH(f);
        allHH(f,:,epocha)=H;
    end
    
    fprintf('%i ',epocha);
end
fprintf('\n'); %konec radku

%grafy
T = 0:0.1:size(EEG.data,2)/EEG.srate; %cas zacatku a konce epochy
meanHH = squeeze(mean(allHH,3)); %prumer pres vsechny epochy
figure('Name','prumerny hilbert');
imagesc(T,freq,meanHH);
axis xy;
hold on;
line([eventlatency eventlatency], [freq(1) freq(end)],'LineWidth',2,'Color','black');
colorbar;
caxis([-0,2.5]);

%statistika
C = mean(allHH(:,1:EEG.event(1,1).latency-1,:),2); %prumer pres cas 
fdr = 1;
W = WilcoxA(allHH,C);
figure('Name','W map hilbert');
imagesc(T,freq, 1-W,[ iff(fdr,0.95,0.99) 1]);%mapa, od p>0.05 bude modra barva 
colorbar;
axis xy;

if Pyes
    wsize = floor(wsize_sec * EEG.srate);
    woverlap = 2; %pro jine hodnoty nechapu pocet hodnot, ktere vraci spectrogram
    allPP = zeros( [  size(freq,2)  floor(size(EEG.data,2)/wsize*woverlap-1) size(EEG.data,3) ] ); %pocet hodnot spektrogramu zkousim odhadnout 
    ffPP = zeros( [ size(freq,2)  size(EEG.data,3) ]); %prumerny spectrogram pro kazdou frekvenci a epochu
    for epocha = 1: size(EEG.data,3);
        %spectrogram
        [~,~,~,P]=spectrogram(double(EEG.data(ch,:,epocha)),wsize,wsize/woverlap,freq,EEG.srate,'yaxis');
        PP = 10*log10(abs(P));
        for f = 1:numel(freq)
            ffPP(f,epocha)=mean(PP(f,:));
        end
    end
    nPP = mean(ffPP,2); %prumer spectrogram pro kazdou frekvenci pres vsechny epochy
    
    fprintf('spectrogram epocha: ');
    for epocha = 1: size(EEG.data,3)
         %spectrogram
        [~,F,T,P]=spectrogram(double(EEG.data(ch,:,epocha)),wsize,wsize/woverlap,freq,EEG.srate,'yaxis');
        PP = 10*log10(abs(P));
        for f = 1:numel(freq)
            PP(f,:) = PP(f,:) ./ nPP(f); 
        end
        allPP(:,:,epocha)=PP; %frekvence=radky, cas=sloupce
    end
    fprintf('\n'); %konec radku

    meanPP = squeeze(mean(allPP,3)); %prumer pres vsechny epochy
    figure('Name','prumerny spektrogram');
    imagesc(T,F,meanPP);
    axis xy;
    hold on;
    line([eventlatency eventlatency], [freq(1) freq(end)],'LineWidth',2,'Color','black');
    colorbar;
    caxis([-0,2.5]);
end



