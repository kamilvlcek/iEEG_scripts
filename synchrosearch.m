function [ ] = synchrosearch( d,fs, mults, channel,start,konec )
%SYNCHROSEARCH skript na hledani synchronizacniho pulsu v datech
%   vykresli jeden kanal a ceka v pauze na mezernik - ukoncim ctrl-c
%   start a konec jsou ve vterinach podobne jako u dataview
%   channel je kanal, kde se zacina hledat

fh = figure('Name','Synchrosearch');
if ~exist('start','var') 
    start = 0;
end
if ~exist('konec','var') 
    konec = size(d,1)/fs;    %konec souboru ve vterinach
end
if ~exist('mults','var') || isempty(mults)
    %mults = ones(1,size(d,2));
else
    d = bsxfun(@times,double(d), mults);
end
delka = konec - start; %delka souboru ve vterinach
disp(['start ' num2str(start) ' delka ' num2str(delka) ' konec ' num2str(konec)]);
pause on;

if ~exist('channel','var') %muzu zadat kanal se synchronizaci - 15.4.2015
    channel = 1; %synchronizace byva 2 kanaly pred koncem - pred EKG
end

%kod od Jirky Hammer na high pass filter - potom jsou videt lip synchropulsy
loF = 0.15; % lower cutoff, in [Hz]
Wn = loF/(fs/2); % normalized bandpass frequencies
n = 3; % butterworth order
[b,a] = butter(n, Wn, 'high'); % returns polynoms of Butterw. filter

for ch = channel:size(d,2)
    x_raw = d(start*fs+1:(start+delka)*fs,ch);
    plot( start:1/fs:start+delka-1/fs, x_raw,'b' );   %puvodni eeg data
    axis([start start+delka -3500 3500])
    
    hold on;
    x_filt = filtfilt(b, a, x_raw);
    plot( start:1/fs:start+delka-1/fs, x_filt,'r' );  %zfiltrovana eeg data
    xlabel('seconds');
    title([' channel ' num2str(ch)]);
    hold off;
    pause;

end

pause off;
end

