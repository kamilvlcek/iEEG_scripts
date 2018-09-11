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
    mults = ones(1,size(d,2));
end
delka = konec - start; %delka souboru ve vterinach
disp(['start ' num2str(start) ' delka ' num2str(delka) ' konec ' num2str(konec)]);
pause on;

if ~exist('channel','var') %muzu zadat kanal se synchronizaci - 15.4.2015
    channel = 1; %synchronizace byva 2 kanaly pred koncem - pred EKG
end


for ch = channel:size(d,2)
    plot( start:1/fs:start+delka-1/fs,  d(start*fs+1:(start+delka)*fs,ch) .* mults(1,ch));   
    axis([start start+delka -3000 3000])
    xlabel('seconds');
    title([' channel ' num2str(ch)]);
    pause;

end

pause off;
end

