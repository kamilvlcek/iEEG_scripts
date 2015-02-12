function [ eegdata, UU, intervalTime,radkyZacatky ] = importIEEG ( d,H, LPT,t, interval) 
%IMPORT  vytvori matlab pole, ktere se da nacist do EEG Lab
%   LPT je cislo kanalu se synchronizacnimi pulsy
%   interval Zajmu jsou timestampy, mezi kterymi me zaznam zajima


radkyZacatky = zeros(size(interval,1),1); %zacatky intervalu [2,sec] = delky predchoziho intervalu
eegdata = zeros(size(d,2),0);
UU = zeros(0,4);
intervalTime = zeros(0,1);
kresli_RT_ITI = 1;
for radka = 1:size(interval,1) %cyklus pro zadane casove intervaly v datech
    %treti parametr = 1=nahoru, 0=dolu
    kresli_udalosti = 0; %1=pauzovat u vsech udalosti, 2 = sum udalosti
    U1 = udalosti(d(:,LPT),H,0,kresli_udalosti, interval(radka,:) ); %casy podnetu - kupodivu u p68 jsou dolu!
    U2 = udalosti(d(:,LPT),H,1,kresli_udalosti,interval(radka,:) ); %casy reakci
    %prvni sloupec je cas od zacatku v sekundach (cislo radku/ vzorkovaci frekvenci)
    %druhy sloupec je denni cas jako timestamp 
    if size(U1,1)~=size(U2,1)
        disp(['radka ' num2str(radka) ': ruzne pocty zacatku a koncu udalosti']);
        keyboard
    end
    U = [U1(:,1) U2(:,1) U2(:,1)-U1(:,1) zeros(size(U1,1),1)]; %k tomu jste reakcni casy
    %sloupce U: cas podnetu od zacatku v sekundach, cas reakce od zacatku v sekundach, reakcni cas,
    %cas ITI pred timto trialem
    U(2:end,4)=U(2:end,1)-U(1:end-1,2); 
    if kresli_RT_ITI
        figure('Name','Udalosti: RT a ITI');
        plot(U(:,3),'-or'); %vykreslim casy udalosti - reakcni cas
        disp([ 'udalosti: ' num2str(size(U,1)) ] );
        hold on;
        plot(U(:,4),'-+'); %vykreslim casy udalosti - cas ITI pred timto trialem
    end
    %zatim jsou casy v sekundach od zacatku 

    %kontrola posledni udalosti s reakcnim casem <0
    %figure;
    %u1 = U(78,1) * H.samplerate(1,1); 
    %u2 = U(78,2) * H.samplerate(1,1);
    %plot(d(u2-500:u1+1500,LPT));

    %mazu udalost s reakcnim casem zapornym
    %iU = U(:,3)<0;
    %U(iU,:)=[];

    %t je cas v sekundach od zacatku zaznamu
    zaznamTSz = datenum(strrep(H.starttime, '.', ':')); %vytvori timestamp zacatku zaznamu
    t(:,2)=zaznamTSz + (t(:,1)-t(1,1))/24/3600;  %pridam druhy sloupec do t - timestamp od zacatku zaznamu
                    %t(1,1) nemusi byt 0 - proto odecitam - 26.5.2014 - p69
    iD = t(:,2)>=interval(radka,1) & t(:,2)<=interval(radka,2);
    eegdata = [ eegdata d(iD,:)' ]; %#ok<AGROW> %pridam EEG dalsi data za predchozi
    intervalTime0 = t(iD,1)-t(1,1); %cas od zacatku zaznamu [HODNE,sec] - vyber pouze interval

    intervalSec = (interval(radka,:)-zaznamTSz)*24*3600; %prevedu rozdil timestampu na sekundy [1,sec]
    U(:,1)=U(:,1)-intervalSec(1)+radkyZacatky(radka,1);
    U(:,2)=U(:,2)-intervalSec(1)+radkyZacatky(radka,1); %chci sekundy od zacatku intervalu zajmu
    
    %pridam udalosti do celkoveho souhrnu udalosti
    UU = [ UU; U]; %#ok<AGROW>
    %cas vsecho hodnot od zacatku PRVNIHO intervalu [sec] - pridam dalsi radky
    intervalTime = [intervalTime; round2( (intervalTime0(:,1)-intervalSec(1)+radkyZacatky(radka,1) ),0.0001)];   %#ok<AGROW>
    
    if size(interval,1) > radka
        radkyZacatky(radka+1,1)=intervalSec(2)-intervalSec(1);
    end
    %jeste potrebuju dat 
end
figure('Name','udalosti v zaznamu');
plot(intervalTime,eegdata(size(d,2)-2,:))
hold on;
for j = 1:size(UU,1) %cervene stimuli
    line([UU(j,1) UU(j,1)],[-2000 2000],'LineWidth',2,'Color','red'); 
end
for j = 1:size(UU,1) %zlute reakce
    line([UU(j,2) UU(j,2)],[-2000 2000],'LineWidth',1,'Color','yellow'); 
end
for j = 1:size(radkyZacatky,1)
    line([radkyZacatky(j,1) radkyZacatky(j,1)],[-3000 3000],'LineWidth',4,'Color','green'); 
end
end

