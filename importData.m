function [eegdata, U, T , R] = importData(H,d,t,sec_oddo,prvninahoru, LPT)
%IMPORTDATA - importuje data do EEGlab formatu
%sec_oddo jsou sekundy od a sekundy do napriklad [3200 3800]
%prvninahoru - jestli prvni synchronizacni puls=podnet jde nahoru=1, nebo dolu=0
%LPT je cislo kanalu se synchronizaci, default  size(d,2)-2, po


%clear all; %smaze promenne
%load('..\pacienti\Nebuchlova\VT2_4.mat'); %d, H,
%load('..\pacienti\Daenemark p68\VT3_1.mat')
%disp('importovano... ');
if ~exist('LPT','var'), LPT=size(d,2)-2; end
pause on;

%vykreslim synchronizaci pulsy - kontrola, jestli mam spravne kanal LPT
kresli = 0;
if kresli 
    figure;
    U = find(d(:,LPT)>2000); %synchronizacni pulsy prekracuji 2000
    cas = H.samplerate(1,1); %sampling rate
    for j = 1:size(U,1);
        u = U(j,1); 
        plot(d(u-500:u+2500,LPT));
        disp([ num2str(j) ' - ' num2str(u/cas) ' s']); %pause;
        pause;
    end
end

%prvni a druhe provedeni AEDist
%interval = [ datenum('9:37:54') datenum('9:46:50'); datenum('9:47:25') datenum('9:53:00') ]; %p68
%p69 JK
%interval = [ datenum('11:11:38.2') datenum('11:17:18'); datenum('11:18:34') datenum('11:24:11') ]; %p69
%interval2 = [1230 1600; 1658 2000]/24/3600 +datenum(strrep(H.starttime, '.', ':')); %vyjadreni v sekundach zaznamu
%datenum vraci cislo, kde 1 znamena 1 den

%p71 Str
%interval = [488 1180; 1320 2610]/24/3600 +datenum(strrep(H.starttime, '.', ':')); %vyjadreni v sekundach zaznamu

%p73 VT6 1.11.2014
%interval = [590 1250; 1370 2010]/24/3600 +datenum(strrep(H.starttime, '.', ':')); %vyjadreni v sekundach zaznamu

%p73 VT6 PPA lokalizer 19.11.2014
%interval = [3750 3750+700]/24/3600 +datenum(strrep(H.starttime, '.', ':')); %vyjadreni v sekundach zaznamu

%p79 PPA lokalizer 18.6.2015
interval = sec_oddo/24/3600 +datenum(strrep(H.starttime, '.', ':')); %vyjadreni v sekundach zaznamu


% prectu si z U1 pomoci datestr(U1(180,2),'HH:MM:SS.FFF')

[ eegdata, U, T , R] = importIEEG ( d,H, LPT,t, interval,prvninahoru );
%save('AEDist.mat','eegdata','U', 'T','R'); %ulozim tyto promene abych je mohl nacist po spusteni EEG labu
