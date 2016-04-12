function [ iSTART,iEND ] = dataview(d,H,tabs,fs, start,konec,  channel)
%DATAVIEW zobrazi synchronizacni puls s polu s anotacema
%   pracuje s datama EEG z Motola
%   delka i start jsou v sekundach
%   nepovinny channel je cislo kanalu se synchronizaci

if ~exist('channel','var') %muzu zadat kanal se synchronizaci - 15.4.2015
    channel = size(d,2)-2; %synchronizace byva 2 kanaly pred koncem - pred EKG
end  
if ~exist('konec','var')
    konec = H.records;
    delka = H.records-start;  %23.6.2015 - kdyz neudam delku, zobrazuje cely zaznam
else
    delka = konec - start;
end
if start==0, start = 1; end;
%fs = H.samplerate(1); %sampling rate
figure('Name','Synchronizace');

%procenta = delka/size(d,1)*100;
%for x = start:delka:H.records; %x je index zacatku
%    if x+delka-1 > H.records
%       break; %kdyz uz delkou presahuju konec zaznamu, ukoncim cyklus
%    end
x = start;
    plot( x:1/fs:x+delka-1/fs,  d(x*fs:(x+delka)*fs-1,channel));
    axis([x x+delka -3000 3000])
    sekund_zac = x; % cas v sekundach,cisla s desetinnymi teckami se do grafu na osu x nevejdou
    sekund_konec = (x+delka);
    %sekund_delka = delka/fs;
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
    
    iSTART = x*fs; %index zacatku v poli d a tabs
    tsSTART = tabs(iSTART); %timestamp zacatku
    disp( ['zacatek: ' num2str(x) 's, timestamp: ' datestr(tsSTART,'dd-mmm-yyyy HH:MM:SS.FFF') ', iSTART: ' num2str(iSTART)]);
    iEND = konec*fs; %index konce v poli d a tabs
    tsEND = tabs(iEND); %timestamp konce
    disp( ['konec: ' num2str(konec) 's, timestamp: ' datestr(tsEND,'dd-mmm-yyyy HH:MM:SS.FFF') ', iEND: ' num2str(iEND)]);
   %keyboard; %zastavi a muzu se divat na promenne, pokracuju pomoci return
%end

end



