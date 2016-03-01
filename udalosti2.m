function [ UU ] = udalosti2( d, fs, tabs, nahoru, kresli, interval )
%UDALOSTI vrati indexy udalosti v poli d a jejich timestampy
%pouziva se pro import dat z Motola do EEGlabu
%  fs je pocet udalosti v jedne s, sampling rate
%  nahoru je 1 pokud je udalost do pozitivnich hodnot
%  prepoklada hodnoty ve sloupci d(:,1)
%  kresli=1 - vykresli synchro kanal se zachycenymi udalostmi
%  interval - muzu nepovinne zadat cast souboru k detekci udalosti - obsahuje timestampy

if ~exist('inverval','var') 
    interval = [ tabs(1) tabs(end)];
end
if ~exist('kresli','var') 
    kresli = 1; %defaultni hodnota
end
LPT = 1; %cislo kanalu synchronize - d ted obsahuje jen jeden sloupec
Th = 2000; %me kterou vsechny udalosti prekracuji
if nahoru 
    U = find(d(:,LPT)>Th);
else
    U = find(d(:,LPT)<-Th); %zatim vsechny hodnoty, ktere prekracuji limit
end

% do druheho sloupce dam rozdil soucasne proti predchozi radce
U(2:end,2)=U(2:end,1)-U(1:end-1); 
U(1,2) = 1000; %prvni hodnotu nechci smazat, nema definovany casovy rozdil vuci predchozi

%smazu udalosti blize nez 0.25 sec od predchozi
iU = U(:,2)< fs/4; %vice nez 0.25 sec
U(iU,:)=[]; 
%smazu udalosti ktere jsou sum - do tretiho a ctvrteho sloupce dam uroven sumu pred a po udalosti
for j = 1:size(U,1)
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

if kresli == 1 %defaultni obrazek
    figure('Name','Udalosti na synchronizacnim pulsu'); 
    plot(d(:,LPT)); %osa x je index v poli d a tabs
    for j = 1:size(U,1)
        line([U(j,1) U(j,1)],[0 iff(nahoru==1,Th,-Th)],'LineWidth',2,'Color',iff(nahoru==1,'red','green'));
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

