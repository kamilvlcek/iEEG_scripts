function [ EEG ] = EEG2Hilbert( EEG, channels, freq )
%EEG2HILBERT prevede vsechny kanaly na prumer hilbertovych obalek
%   freq je seznam freqvenci pro ktere se ma delat prumer - lo, .., ..., .., hi

% http://www.scholarpedia.org/article/Hilbert_transform_for_brain_waves
for ch = channels
    fprintf('channel %i: Hz ',ch);
    HH = zeros(numel(freq)-1,size(EEG.data,2));
    %figure;
    for fno = 1:numel(freq)-1
        loF = freq(fno);
        hiF = freq(fno+1)-1;
        
        hh  = hilbertJirka(double(EEG.data(ch,:)),loF,hiF,EEG.srate);
        HH(fno,:) = (hh./mean(hh)); %podil prumeru 1 = prumerna hodnota
        fprintf('%i ',freq(fno));
        %disp(freq(fno));
        %plot(HH(fno,1:5000));
        %hold on;
    end   
    
    M = mean(HH,1);
    %plot(M(1:5000),'r');
    
    EEG.data(ch,:) = M;
    disp('OK');
end

end

