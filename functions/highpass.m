function [ x_filt ] = highpass( loF,srate,x_raw)
%HIGHPASS high pass filter od Jirky Hammera
%   Pouzivam celkem standardni high-pass filter pomoci IIR Butterworth filtru, v matlabu takto 
%   (x_raw jsou vstupni data, vektor nebo matice, v pripade matice to filtruje po sloupcich):


%loF = 0.15; % lower cutoff, in [Hz]
Wn = loF/(srate/2); % normalized bandpass frequencies
n = 3; % butterworth order
[b,a] = butter(n, Wn, 'high'); % returns polynoms of Butterw. filter
x_filt = filtfilt(b, a, x_raw); 

%figure; hold on; plot(x_raw, 'r'); plot(x_filt, 'b'); 

end

