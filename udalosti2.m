function [ UU ] = udalosti2( d, fs, tabs, nahoru, mults,kresli, interval, diffsec,threshold,exclude)
%UDALOSTI vrati indexy udalosti v poli d a jejich timestampy
%pouziva se pro import dat z Motola do EEGlabu
%udalosti2( d, fs, tabs, nahoru, mults,kresli, interval, diffsec,threshold)
%  fs je pocet udalosti v jedne s, sampling rate
%  nahoru je 1 pokud je udalost do pozitivnich hodnot
%  prepoklada hodnoty ve sloupci d(:,1)
%  kresli=1 - vykresli synchro kanal se zachycenymi udalostmi
%  interval - muzu nepovinne zadat cast souboru k detekci udalosti - obsahuje timestampy
%  diffsec - nejmensi interval mezi nasledujicimi udalostim v sekundach, default 0.125
%  threshold - vyska trigerovaciho pulsu, default je 2000

if ~exist('inverval','var')  || isempty(interval)
    interval = [ tabs(1) tabs(end)];
end
if ~exist('kresli','var') || isempty(kresli)
    kresli = 1; %defaultni hodnota
end
if ~exist('mults','var') || numel(mults)<size(d,2)
    mults = ones(1,size(d,2));
end
if ~exist('diffsec','var') || isempty(diffsec)
    diffsec = 0.125; %nejmensi interval mezi nasledujicimi udalostim
end
if ~exist('threshold','var')  || isempty(threshold)
    threshold = 2000; %hodnota kterou vsechny udalosti prekracuji
end
if ~exist('exclude','var') 
    exclude = []; %hodnota kterou vsechny udalosti prekracuji
end
if size(d,2) > 1    
    LPT = size(d,2)-2; %synchronizace byva 2 kanaly pred koncem - pred EKG
    d = double(d(:,LPT)) .* mults(LPT); %predpokladam d a mults se stejnym poctem kanalu 
    LPT = 1; %zredukoval jsem d na jen jednorozmerny vektor ze synchronizacniho kanalu
else
    LPT = 1; %cislo kanalu synchronize - d ted obsahuje jen jeden sloupec
    d = double(d) * mults(LPT); 
end

if nahoru 
    U = find(d(:,LPT)>threshold);
else
    U = find(d(:,LPT)<-threshold); %zatim vsechny hodnoty, ktere prekracuji limit
end

assert( size(U,1) > 1, ['nenalezeny zadne synchronizacni udalosti v synchro kanalu ']);  %#ok<NBRAK>
% do druheho sloupce dam rozdil soucasne proti predchozi radce
U(2:end,2)=U(2:end,1)-U(1:end-1); 
U(1,2) = 1000; %prvni hodnotu nechci smazat, nema definovany casovy rozdil vuci predchozi

%smazu udalosti blize nez 0.125 sec od predchozi - nasel jsem spravnou udalost 0.23 od predchozi
iU = U(:,2)/fs < diffsec; % minimalni casovy interval
U(iU,:)=[]; 
%smazu udalosti ktere jsou sum - do tretiho a ctvrteho sloupce dam uroven sumu pred a po udalosti
for j = 1:size(U,1)
    %pole 3-5 jsou tu jen kvuli obrazku kresli=2
    U(j,3)=std(d(  U(j,1)-200:U(j,1)-100  ,LPT)); %stdev 100 bodu pred zacatkem udalosti
    U(j,4)=std(d(  U(j,1)+100:U(j,1)+200  ,LPT)); %stdev 100 bodu 100 po zacatku udalosti
    U(j,5)=mean(d(  U(j,1)+20:U(j,1)+80  ,LPT)); %prumer 60 bodu hned po zacatku udalosti - je tam spicka dolu?
    U(j,6)=min(d(  U(j,1):U(j,1)+100  ,LPT)); %minimum 100 bodu hned po zacatku udalosti - je tam spicka dolu?
    U(j,7)=tabs(U(j,1)); %caszacatku + U(j,1)/fs/24/3600; %cas udalosti h:m:s jako float
end
%iU =  U(:,3)>200  | U(:,4)> 200; %| U(:,6) < -1500; %U(:,3)<5 || U(:,4)<5 
%U(iU,:)=[]; %**********

%vymazu udalosti mimo zvolene casovy interval
iU = U(:,7)<interval(1) | U(:,7)>interval(2);
U(iU,:)=[]; %**********

if ~isempty(exclude)
    U(exclude,:) = [];
end
if kresli == 1 %defaultni obrazek
    figure('Name','Udalosti na synchronizacnim pulsu'); 
    subplot(2,1,1); %prvni obrazek se sychronizacnimi pulsy a pres nej udalostmi
    secs = [1:size(d,1)]'/fs; %#ok<NBRAK> %sekundy zaznamy
    plot(secs,d(:,LPT)); %osa x jsou sekundy zaznamu - 26.4.2017 drive index v poli d a tabs
    hold on;
    for j = 1:size(U,1)
        line([U(j,1) U(j,1)]/fs,[0 iff(nahoru==1,threshold,-threshold)],'LineWidth',2,'Color',iff(nahoru==1,'red','green'));
        %text( U(j,1)/fs, iff(nahoru==1,threshold,-threshold)/mod(j,5) ,num2str(j),'Color',iff(nahoru==1,'red','green'));
        % ten text to hrozne zdrzuje, tak si ho necham v komentari
    end
    title('synchronizacni puls a na nem udalosti');
    xlabel('sec');
    
    subplot(2,1,2); %druhy obrazek s intervaly mezi pulsy
    plot(U(:,1)/fs,U(:,2)/fs,'.-'); % 4.7.2017 - na ose x jsou vteriny zaznamu
    ylim([-0.1 10]); 
    %xlim([-10 size(U,1)+10]);
    title('intervaly mezi udalostmi v sekundach');
    ylabel('sec');
    %jeste vypisu do spodniho graf pocty pulsu v blocich
    blokyzac = [1 ; find(U(:,2)./512>5)]; %zacatky bloku, oddelenych delsi casovou mezerou (5 s)
    blokydelka = [blokyzac(2:end) ; size(U,1)+1]-blokyzac(1:end); %delka bloku = konce -  zacatky, jeden konec pridavam jako velikost
    blokyzac(:,2) = U(blokyzac(:,1),1)./fs; %cas zacatku bloku ve vterinach
    for b = 1:size(blokyzac,1)
        text(blokyzac(b,2),5,num2str(blokydelka(b)));
        text(blokyzac(b,2),6,num2str(blokyzac(b,1)));
    end
end;


%obrazek sumu vsech udalosti
if kresli==2
    %asi hlavne ladici obrazek pro detekci udalosti, zobrazuje std, std, mean
    figure('Name','sum pred, sum po, prumer hned po');
    plot(U(:,3:5),'-O');
end

UU = [U(:,1), U(:,7)]; %index v poli d a timestamps
end

