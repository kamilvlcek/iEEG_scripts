function [ iSTART,iEND,pauzyvdatech ] = dataview(d,tabs,fs,mults, start,konec,evts, annotations,  channel,Events)
%DATAVIEW zobrazi synchronizacni puls s polu s anotacema
%   pracuje s datama EEG z Motola
%   delka i start jsou v sekundach
%   nepovinny channel je cislo kanalu se synchronizaci
%  (d,tabs,fs,mults,start,konec, evts, annotations,  channel,evts,Events)
if size(d,2)==1
    channel = 1; %pokud je d jednorozmerne, neresim cislo kanalu
elseif ~exist('channel','var') || isempty(channel) %muzu zadat kanal se synchronizaci - 15.4.2015
    channel = size(d,2)-2; %synchronizace byva 2 kanaly pred koncem - pred EKG
end  
if ~exist('start','var') || start == 0
    start = 0;
    iStart = 1; %index v poli tabs
else
    iStart = round(start*fs); %index v poli tabs
end
zaznamvterin = (tabs(end)-tabs(1))*24*3600; %delka zaznamu ve vterinach
if ~exist('konec','var') || konec == 0
    konec = floor(zaznamvterin); 
    delka = konec -start;  %23.6.2015 - kdyz neudam delku, zobrazuje cely zaznam
    iKonec = size(tabs,1);%index v poli tabs
else
    delka = konec - start;
    iKonec = min(size(tabs,1),round(konec*fs));%index v poli tabs
end
if ~exist('mults','var') || numel(mults)<size(d,2)
    mults = ones(1,size(d,2));
end

%fs = H.samplerate(1); %sampling rate

pauzyvdatech = kontrolacasu(tabs(iStart:iKonec),fs);
if numel(pauzyvdatech) > 0
    disp('chyby datech mezi zacatkem a koncem');
    m=input('Do you want to continue, y/n [n]:','s');
    if m~='y'  
        return; 
    else
        iP = pauzyvdatech(:,1)<iKonec & pauzyvdatech(:,1)>=iStart;
        delka = delka - round(sum((tabs(pauzyvdatech(iP,1)+1)-tabs(pauzyvdatech(iP,1)))*24*3600));        
    end    
end

figure('Name','Synchronizace');
%procenta = delka/size(d,1)*100;
%for x = start:delka:H.records; %x je index zacatku
%    if x+delka-1 > H.records
%       break; %kdyz uz delkou presahuju konec zaznamu, ukoncim cyklus
%    end
x = start;
data = d(x*fs+1:(x+delka)*fs,channel) .* mults(1,channel);
ymax = 3000;
if sum(data>3000)>5e4
    ymax = 5e5;
    disp('range increased to 5e5');
end
yrange = [-ymax ymax];
    xtime = x:1/fs:x+delka-1/fs; %cas ve vterinach - osa x zobrazenych dat
    plot( xtime, data ); %,'-o'
    axis([x x+delka yrange])    
    sekund_zac = x; % cas v sekundach,cisla s desetinnymi teckami se do grafu na osu x nevejdou
    sekund_konec = (x+delka);
    %sekund_delka = delka/fs;
    %set(gca,'XTickLabel',cellstr(int2str( ( sekund_zac : (sekund_delka/10) : sekund_konec)'))');
    xlabel(['sekundy z ' num2str(zaznamvterin) ' s celkove']); 
    ylabel(['channel ' num2str(channel)]); 
    zobrazenych = 0;
    if exist('annotations','var') && isstruct(annotations)
        for a = 1:numel(annotations.starttime)
            if annotations.starttime(a)>sekund_zac && annotations.starttime(a)<sekund_konec
                anotace_index = annotations.starttime(a);
                line([anotace_index anotace_index],yrange,'Color','red');
                text(anotace_index,2900-zobrazenych*80,annotations.event{a},'Color','red');
                zobrazenych = zobrazenych+1;
            end
        end
    end
    if exist('evts','var') && isstruct(evts)
        [~, idx] = unique({evts.dateStr}',  'stable'); %'rows', jen unikatni dateStr ve strukture. 
        evts = evts(idx);
        for a = 1:numel(evts)
            secs = find(tabs >= datenum(evts(a).dateStr),1) / fs; %v kolika vterichan od zacatku tabs
            if ~isempty(secs) && secs>sekund_zac && secs<sekund_konec && isfield(evts(a),'annotation') && ~isempty(evts(a).annotation)                
                line([secs secs],yrange,'Color','red');
                text(secs,2900-zobrazenych*80,evts(a).annotation,'Color','red');
                text(secs,2000-zobrazenych*80, evts(a).dateStr,'Color','magenta');
                zobrazenych = zobrazenych+1;
            end
        end
    end
    if exist('Events','var') && isstruct(Events)
        for k = 1:numel(Events.c.timestamps)
            secs = find(tabs >= Events.c.timestamps(k),1) / fs;
            line([secs secs],yrange,'Color','black');
            text(secs,100,['c' num2str(k)],'Color','black');
        end
        for k = 1:numel(Events.g.timestamps)
            secs = find(tabs >= Events.g.timestamps(k),1) / fs;
            line([secs secs],yrange,'Color','red');
            text(secs,200,['g' num2str(k)],'Color','red');
        end
        for k = 1:numel(Events.e.timestamps)
            secs = find(tabs >= Events.e.timestamps(k),1) / fs;
            line([secs secs],yrange,'Color','magenta');
            text(secs,300,'e','Color','magenta');
        end
    end
    iSTART = x*fs+1; %index zacatku v poli d a tabs
    tsSTART = tabs(iSTART); %timestamp zacatku
    disp( ['zacatek: ' num2str(x) 's, timestamp: ' datestr(tsSTART,'dd-mmm-yyyy HH:MM:SS.FFF') ', iSTART: ' num2str(iSTART)]);
    iEND = konec*fs; %index konce v poli d a tabs
    if(numel(pauzyvdatech))>0
        iEND = iEND - floor(sum( (tabs(pauzyvdatech(iP)+1)-tabs(pauzyvdatech(iP)))*24*3600*fs) );
    end
    tsEND = tabs(iEND); %timestamp konce
    disp( ['konec: ' num2str(konec) 's, timestamp: ' datestr(tsEND,'dd-mmm-yyyy HH:MM:SS.FFF') ', iEND: ' num2str(iEND)]);
    %keyboard; %zastavi a muzu se divat na promenne, pokracuju pomoci return
    %casove zna�ky po 10 minut�ch
    tabsX = tabs(iSTART:iEND); %tabs odpovidajici ose X=zobrazenym datum 
    for ix = 1:600*fs:numel(tabsX)        
%         text(xtime(ix),2500,datestr(tabs(ix),'dd-mmm-yyyy HH:MM:SS.FFF'),'Color','black');
        text(xtime(ix),2500,datestr(tabsX(ix),'dd-mmm-yyyy HH:MM:SS.FFF'),'Color','black');        
    end
    text(xtime(end),1000,datestr(tabsX(end),'dd-mmm-yyyy HH:MM:SS.FFF'),'Color','black');
    if numel(pauzyvdatech) > 0
       for ip = 1:size(pauzyvdatech,1)
           line([pauzyvdatech(ip,2) pauzyvdatech(ip,2)],yrange,'Color','red');
           text(pauzyvdatech(ip,2),2000-ip*200,datestr(tabs(pauzyvdatech(ip,1)+1),'dd-mmm-yyyy HH:MM:SS.FFF'),'Color','red');
           text(pauzyvdatech(ip,2),2090-ip*200,datestr(tabs(pauzyvdatech(ip,1)),'dd-mmm-yyyy HH:MM:SS.FFF'),'Color','black');
           text(pauzyvdatech(ip,2),1000-ip*200,num2str(pauzyvdatech(ip,3),'%.2f sec'),'Color','red');
       end       
    end
%end

end

function [dlouhe] = kontrolacasu(tabs,fs)
tabsdiff = diff(tabs)*24*3600*1000;       %rozdily timestampu v ms
dlouhe = find(tabsdiff > (1/fs*1000*2) | tabsdiff<0); %najit rozdily v casech vetsi nez dvojnasobek vzorkovaci frekvence
dlouhe(:,2) = dlouhe(:,1)/fs;
dlouhe(:,3) = (tabs(dlouhe(:,1)+1)-tabs(dlouhe(:,1)))*24*3600; %delka rozdilu v s
if numel(dlouhe)  > 0
    disp(['pauzy v datech: ' num2str(numel(dlouhe))]);
    disp(['prvni pauza: i=' num2str(dlouhe(1)) '=' num2str(dlouhe(1)/fs) 's, timediff=' num2str(tabsdiff(dlouhe(1))) ' ms']);    
end
end



