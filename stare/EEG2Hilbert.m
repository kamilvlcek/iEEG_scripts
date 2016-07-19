function [ EEG ] = EEG2Hilbert( EEG, channels, freq )
%EEG2HILBERT prevede vsechny kanaly na prumer hilbertovych obalek
%   freq je seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi
%   EEG je dataset z EEGlabu

% http://www.scholarpedia.org/article/Hilbert_transform_for_brain_waves
if size(EEG.data,3) > 1
   %EEG.data = squeeze(EEG.data(:,:,1)); %pokud uz rozdelene na epochy, beru jen prvni epochu
   kresli = 1; %pokud mam rozdelena data na epochy, chci kreslit spektrogram prvni epochy
else
   kresli = 0;
end

for ch = channels %jednotlive elektrody
    fprintf('channel %i: Hz ',ch);
    HH = zeros(numel(freq)-1,size(EEG.data,2));
    HT = zeros([size(EEG.data,2) numel(freq) size(EEG.data,3)]); %o jeden rozmer vic nez EEG data - pro epochovana data: time, freq, epoch,
    
    for fno = 1:numel(freq)-1 %seznam frekvenci
        loF = freq(fno);
        hiF = freq(fno+1)-1;

        for epoch = 1:size(EEG.data,3) %epochy u rozdelenych dat, u kontinualnich pouze jednou
            hh  = hilbertJirka(double(EEG.data(ch,:,epoch)),loF,hiF,EEG.srate);
            HT(:,fno,epoch) = hh;
            
            HH(fno,:) = HH(fno,:) + (hh./mean(hh)).*100; %podil prumeru 100 = prumerna hodnota
        end
        HH(fno,:) = HH(fno,:) ./ size(EEG.data,3); %pokud vic epoch - secital jsem obalky a ted pocitam prumer
        fprintf('%i ',freq(fno));
        %disp(freq(fno));
        %plot(HH(fno,1:5000));
        %hold on;
    end
        
    wsize_sec = 0.125; %velikost okna pro vypocet frekvence ve vterinasch
    wsize = wsize_sec * EEG.srate;  %velikost okna v poctu hodnot
    [~,~,~,P] = spectrogram(double(EEG.data(ch,:,1)),wsize,wsize/2,freq,EEG.srate,'yaxis');
        %P ma rozmery freq x pocethodnot/128
    spectrogram(double(EEG.data(ch,:,1)),wsize,wsize/2,5:10:145,EEG.srate,'yaxis'); 
    htplot(HT,1);
    
    M = mean(HH,1);
    if kresli
        figure('Name',['Channel ' num2str(ch)]);
        imagesc(HH); % spektrogram prvni epochy
        title(['Channel ' num2str(ch)]);
    else
        EEG.data(ch,:) = M; %ukladat obalky chci jen u dat nerozdelenych na epochy
    end
    %plot(M(1:5000),'r');
       
    disp('OK');
end

end

function [] = htplot(HT, epoch)

htf = figure('Name','Hilbert freq'); %#ok<NASGU>
for freq = 1:size(HT,3) %frekvence
    plot(HT(:,freq,epoch));
    hold on;
end

end



