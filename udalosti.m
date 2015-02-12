function [ UU ] = udalosti( d, H, nahoru, kresli, interval )
%UDALOSTI vrati casy udalosti v poli d
%  cas je pocet udalosti v jedne s, sampling rate
%  nahoru je 1 pokud je udalost do pozitivnich hodnot
%  prepoklada hodnoty ve sloupci 1 
%global d;
LPT = 1; %cislo kanalu synchronize - d ted obsahuje jen jeden sloupec
Th = 2000; %me kterou vsechny udalosti prekracuji
cas = H.samplerate(1,1);
caszacatku = datenum(strrep(H.starttime, '.', ':')); 
if nahoru 
    U = find(d(:,LPT)>Th);
else
    U = find(d(:,LPT)<-Th); %zatim vsechny hodnoty, ktere prekracuji limit
end

% do druheho sloupce dam rozdil soucasne proti predchozi radce
U(2:end,2)=U(2:end,1)-U(1:end-1); 

%smazu udalosti blize nez 0.25 sec od predchozi
iU = U(:,2)< cas/4; %vice nez 0.25 sec
U(iU,:)=[]; 
%smazu udalosti ktere jsou sum - do tretiho a ctvrteho sloupce dam uroven sumu pred a po udalosti
for j = 1:size(U,1)
    U(j,3)=std(d(  U(j,1)-200:U(j,1)-100  ,LPT)); %stdev 100 bodu pred zacatkem udalosti
    U(j,4)=std(d(  U(j,1)+100:U(j,1)+200  ,LPT)); %stdev 100 bodu 100 po zacatku udalosti
    U(j,5)=mean(d(  U(j,1)+20:U(j,1)+80  ,LPT)); %prumer 60 bodu hned po zacatku udalosti - je tam spicka dolu?
    U(j,6)=min(d(  U(j,1):U(j,1)+100  ,LPT)); %minimum 100 bodu hned po zacatku udalosti - je tam spicka dolu?
    U(j,7)=caszacatku + U(j,1)/cas/24/3600; %cas udalosti h:m:s jako float
end
iU =  U(:,3)>200  | U(:,4)> 200; %| U(:,6) < -1500; %U(:,3)<5 || U(:,4)<5 
    %kupodivu to porad mase spatne udalosti Daenemark AEDIst - 7.3.2014
%U(iU,:)=[]; %**********

%vymazu udalosti mimo zvolene casovy interval
iU = U(:,7)<interval(1) | U(:,7)>interval(2);
U(iU,:)=[]; %**********


%obrazek sumu vsech udalosti
if kresli==2
    figure('Name','sum pred, sum po, prumer hned po');
    plot(U(:,3:5),'-O');
end

%pause on;
if kresli == 1
    figure(); 
    %hold on; %**********
end;

for j = 1:size(U,1)-1;
    u = U(j,1);
    u1= U(j+1,1);
    if kresli == 1
        plot(d(u-500:u1+500,LPT));
        casden = datestr(caszacatku+u/cas/24/3600,'HH:MM:SS.FFF'); %num2str(u/cas)
        casdiff = num2str(U(j,2)/cas);
        disp([ num2str(j) ' - ' casden ' s, ' casdiff 's diff - sum:' num2str(U(j,3)) ' ' num2str(U(j,4)) ]); 
        pause; %**********
    end
end

UU = [U(:,1)/cas, U(:,7)];
end

