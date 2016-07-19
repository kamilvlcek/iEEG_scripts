function [ freqPow ] = hilbertJirka( rawData, loF, hiF, srate )
%HILBERTJIRKA hilbertJirka( rawData, loF, hiF, srate )
%   vrati Power vybraneho frekvencniho pasma

% Hilberova obalka podle Jirky Hammera - mail 1.8.2014
    %loF = pasmo; hiF = pasmo + 10; srate = 1000;
    %rawData = yy;

    %1) filrovani v gamma pasmu: loF = 60, hiF = 100, srate = sampling rate
    Wn = [loF, hiF]/(srate/2); % normalized bandpass frequencies
    n = 4; % butterworth order
    [b,a] = butter(n, Wn); % returns polynoms of Butterw. filter
    filtData = filtfilt(b, a, rawData);

     %2) analyticky signal pomoci Hilbertovy transformace:
    tmp = hilbert(filtData);

     %3) gamma amplitude, power
    %gammaAmp = abs(tmp);
    freqPow = abs(tmp).^2;
end

