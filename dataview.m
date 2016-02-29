function [ tsZAC ] = dataview(d,H,tabs,start,konec, channel)
%DATAVIEW zobrazi synchronizacni puls s polu s anotacema
%   pracuje s datama EEG z Motola
%   delka i start jsou v sekundach
%   nepovinny channel je cislo kanalu se synchronizaci

if ~exist('channel','var') %muzu zadat kanal se synchronizaci - 15.4.2015
    channel = size(d,2)-2; %synchronizace byva 2 kanaly pred koncem - pred EKG
end  
if ~exist('konec','var')
    delka = H.records-1;  %23.6.2015 - kdyz neudam delku, zobrazuje cely zaznam
else
    delka = konec - start;
end
if start==0, start = 1; end;

figure('Name','Synchronizace');

%procenta = delka/size(d,1)*100;
%for x = start:delka:H.records; %x je index zacatku
%    if x+delka-1 > H.records
%       break; %kdyz uz delkou presahuju konec zaznamu, ukoncim cyklus
%    end
x = start;
    plot( x:1/H.samplerate(1):x+delka-1/H.samplerate(1),  d(x*H.samplerate(1):(x+delka)*H.samplerate(1)-1,channel));
    axis([x x+delka -3000 3000])
    sekund_zac = x; % cas v sekundach,cisla s desetinnymi teckami se do grafu na osu x nevejdou
    sekund_konec = (x+delka);
    %sekund_delka = delka/H.samplerate(1);
    %set(gca,'XTickLabel',cellstr(int2str( ( sekund_zac : (sekund_delka/10) : sekund_konec)'))');
    xlabel(['sekundy z ' num2str(H.records) ' s celkove']); 
    zobrazenych = 0;
    for a = 1:numel(H.annotation.starttime)
        if H.annotation.starttime(a)>sekund_zac && H.annotation.starttime(a)<sekund_konec
            anotace_index = H.annotation.starttime(a);
            line([anotace_index anotace_index],[-3000 3000],'Color','red');
            text(anotace_index,2900-zobrazenych*80,H.annotation.event{a},'Color','red');
            zobrazenych = zobrazenych+1;
        end
    end
    tsZAC = sec2tabs(x,H.samplerate(1),tabs); %timestamp zacatku
    disp( ['zacatek: ' num2str(x) 's, timestamp: ' datestr(tsZAC,'dd-mmm-yyyy HH:MM:SS.FFF') ', ' num2str(tsZAC,'%11.15f')]);
   %keyboard; %zastavi a muzu se divat na promenne, pokracuju pomoci return
%end

end

function [timestamp]=sec2tabs(sec, fs,tabs)
    timestamp = tabs(sec*fs);
end


