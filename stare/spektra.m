function [allHH,T,F]=spektra(EEG, EEGdata, ch,freq, wilcoxontest, obrazky)
%zobrazi spekrogramy pomoci Hilberovy obalky, nemeni referenci, predpoklada bipolarni
%napøíklad 
%[allHH,T,F]=spektra(EEG, EEG.data, 23,10:10:150,1,1);
%Vstupy:
%EEG - data z EEGlabu - record
%EEGdata - 3D matrix -  kanály cas epochy
%ch - èíslo kanálu
%freq - seznam frekvenèních pásem 
%wilcoxontest - bool
%obrazky - bool
%
%Výstupy:
%allHH - 3D matrix - frekvence cas epochy
%T - èas - výstypu z funkce spektrogram, pokus se provádí
%F - frekvence - výstypu z funkce spektrogram, pokus se provádí
%
%24.4.2015 - zkousim porovnat vysledky funkce spectrogram ze spektrem z EEGlabu


Pyes = 0; %jestli provadet analyzu funkci spectrogram
%wilcoxontest = 0;
%ch = 23;
%freq = 10:10:150;
eventlatency = EEG.event(1,1).latency/EEG.srate;
allHH = zeros( [  size(freq,2)-1  size(EEGdata,2) size(EEGdata,3)] ); % frekvence cas epochy
T = 0:0.1:size(EEGdata,2)/EEG.srate; %cas zacatku a konce epochy
F = freq; %vystupni parametr
    
%vypocet prumerne sily ve frekvencnich pasmech
fprintf('prumery za frekvence:');
ffHH = zeros( [ size(freq,2)-1 size(EEGdata,3)]); %prumerny hilbert za kazdou frekvenci a epochu

for epocha = 1: size(EEGdata,3);
    %hilbert
    for f = 1:numel(freq)
        if f < numel(freq) %hilbert je vzdy od do, takze pro posledni neni definovany
            ffHH(f,epocha)=mean(hilbertJirka(double(EEGdata(ch,:,epocha)),freq(f), freq(f+1),EEG.srate));
        end
    end
end

nHH = mean(ffHH,2); %prumer hilbert pro kazdou frekvenci pres vsechny epochy
fprintf('done\n');

% prubehy frekvenci v case
fprintf('hilbert epocha: ');
for epocha = 1: size(EEGdata,3)
    %hilbert
    for f = 1:numel(freq)-1
        H = hilbertJirka(double(EEGdata(ch,:,epocha)),freq(f), freq(f+1),EEG.srate);
        H = H ./ nHH(f);
        allHH(f,:,epocha)=H;
    end
    
    if mod(epocha,10)==0, fprintf(' %i',epocha); end
end
fprintf('\n'); %konec radku
  
%grafy
if obrazky
    meanHH = squeeze(mean(allHH,3)); %prumer pres vsechny epochy, frekvence x cas
    figure('Name','prumerny hilbert');
    imagesc(T,freq,meanHH);
    axis xy;
    hold on;
    line([eventlatency eventlatency], [freq(1) freq(end)],'LineWidth',2,'Color','black');
    colorbar;
    caxis([-0,2.5]);
    ylabel('freq [Hz]');
    xlabel('time [s]');
    title([EEG.setname ': ' num2str(ch)]);
    
    meanHH = squeeze(mean(allHH,1))'; %prumer pres vsechny frekvence, otocene: epochy x cas
    figure('Name','prumerna frekvence');
    imagesc(T,1:size(EEGdata,3),meanHH);
    axis xy;
    hold on;
    line([eventlatency eventlatency], [1 size(EEGdata,3)],'LineWidth',2,'Color','black');
    colorbar;
    caxis([-0,2.5]);
    ylabel('epochs');
    xlabel('time [s]');
    title([EEG.setname ': ' num2str(ch)]);
end

%statistika
if wilcoxontest
    C = mean(allHH(:,1:EEG.event(1,1).latency-1,:),2); %prumer pres cas 
    fdr = 1;
    W = WilcoxA(allHH,C);
    figure('Name','W map hilbert');
    imagesc(T,freq, 1-W,[ iff(fdr,0.95,0.99) 1]);%mapa, od p>0.05 bude modra barva 
    colorbar;
    axis xy;
    ylabel('freq [Hz]');
    xlabel('time [s]');
    title([EEG.setname ': ' num2str(ch)]);
end

if Pyes
    wsize_sec = 0.06; %#ok<UNRCH> %sirka okna spektrogramu v sec 
    wsize = floor(wsize_sec * EEG.srate);
    normalizovat = 0; %jestli se maji frekvence spectrogramu normalizovat na prumer
    woverlap = 2; %pro jine hodnoty nechapu pocet hodnot, ktere vraci spectrogram
    allPP = zeros( [  size(freq,2)  floor(size(EEGdata,2)/wsize*woverlap-1) size(EEGdata,3) ] ); %pocet hodnot spektrogramu zkousim odhadnout 
    
    if normalizovat
        ffPP = zeros( [ size(freq,2)  size(EEGdata,3) ]);  %prumerny spectrogram pro kazdou frekvenci a epochu 
        for epocha = 1: size(EEGdata,3);
            %spectrogram
            [~,~,~,P]=spectrogram(double(EEGdata(ch,:,epocha)),wsize,wsize/woverlap,freq,EEG.srate,'yaxis');
            PP = 10*log10(abs(P));
            for f = 1:numel(freq)
                ffPP(f,epocha)=mean(PP(f,:));
            end
        end
        nPP = mean(ffPP,2); %prumer spectrogram pro kazdou frekvenci pres vsechny epochy
    end
    
    fprintf('spectrogram epocha: ');
    for epocha = 1: size(EEGdata,3)
         %spectrogram
        [~,F,T,P]=spectrogram(double(EEGdata(ch,:,epocha)),wsize,wsize/woverlap,freq,EEG.srate,'yaxis');
        PP = 10*log10(abs(P));
        if normalizovat 
            for f = 1:numel(freq) 
                PP(f,:) = PP(f,:) ./ nPP(f); 
            end
        end
        allPP(:,:,epocha)=PP; %frekvence=radky, cas=sloupce
    end
    fprintf('\n'); %konec radku

    if obrazky
        meanPP = squeeze(mean(allPP,3)); %prumer pres vsechny epochy
        figure('Name','prumerny spektrogram');
        imagesc(T,F,meanPP);
        axis xy;
        hold on;
        line([eventlatency eventlatency], [freq(1) freq(end)],'LineWidth',2,'Color','black');
        colorbar;
        if normalizovat, caxis([-0,2.5]); end; 
    end
end