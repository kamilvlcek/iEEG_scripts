function [ ] = synchrosearch( d,fs,start,konec, channel )
%SYNCHROSEARCH skript na hledani synchronizacniho pulsu v datech
%   vykresli jeden kanal a ceka v pauze na mezernik - ukoncim ctrl-c
%   start a konec jsou ve vterinach podobne jako u dataview
%   channel je kanal, kde se zacina hledat

figure('Name','Synchrosearch');
delka = konec - start;
pause on;

if ~exist('channel','var') %muzu zadat kanal se synchronizaci - 15.4.2015
    channel = 1; %synchronizace byva 2 kanaly pred koncem - pred EKG
end 

for ch = channel:size(d,2)
    plot( start:1/fs:start+delka-1/fs,  d(start*fs:(start+delka)*fs-1,ch));   
    axis([start start+delka -3000 3000])
    text(start,1000,[' channel ' num2str(ch)]);
    pause;

end

pause off;