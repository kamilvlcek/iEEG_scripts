%24.4.2015 - zkousim porovnat vysledky funkce spectrogram ze spektrem z EEGlabu
% viz d:\EEG\motol\pacienti\p68 Daenemark\p68 spectra.cdr
wsize_sec = 0.06;
ch = 43;
freq = 10:10:150;

wsize = floor(wsize_sec * EEG.srate);
woverlap = 2; %pro jine hodnoty nechapu pocet hodnot, ktere vraci spectrogram
allPP = zeros( [ size(EEG.data,3) size(freq,2)  floor(size(EEG.data,2)/wsize*woverlap-1) ] ); %pocet hodnot spektrogramu zkousim odhadnout 
allHH = zeros( [ size(EEG.data,3) size(freq,2)-1  size(EEG.data,2)] ); %epochy frekvence cas

ffHH = zeros( [size(EEG.data,3) size(freq,2)-1]); %primerny hilbert za kazdou epochu a frekvence
for epocha = 1: size(EEG.data,3);
    for f = 1:numel(freq)-1
        ffHH(epocha,f)=mean(hilbertJirka(double(EEG.data(ch,:,epocha)),freq(f), freq(f+1),EEG.srate));
    end
end
nHH = mean(ffHH,1); %prumer hilbert pro kazdou frekvenci
fprintf('ep ');
for epocha = 1: size(EEG.data,3);
    %spectrogram
    [S,F,T,P]=spectrogram(double(EEG.data(ch,:,epocha)),wsize,wsize/woverlap,freq,EEG.srate,'yaxis');
    PP = 10*log10(abs(P));
    allPP(epocha,:,:)=PP;
    
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

meanHH = squeeze(mean(allHH,1));
figure('Name','prumerny hilbert');
imagesc(T,F,meanHH);
axis xy;

disp('hotovo');
