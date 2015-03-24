% skript na zmenseni vzorkovaci frekvence dat z Motola
clear all;
%dir = 'D:\eeg\motol\pacienti\bezdek p49\';
dir = 'D:\eeg\motol\pacienti\Sruma\';
%soubor = 'VT7_arena'; %maly soubor
%soubor = 'VT7_arena1'; % 4GB soubor!!
soubor = 'Sruma INV wifi 8kHz EEG 3 test Arena_2'; %3GB soubor !!

podil = 8; %kolikrat snizit frekvenci
frek = '1000'; %vysledna frekvence

%prvni pulka souboru - s celym najednou se spatne pracuje
load([ dir soubor '.mat']); 
dpul = floor(size(d,1)/(podil*2))*podil; % aby delitelne 16ti
els = size(d,2);
d(dpul+1:end,:)=[]; %smazu druhou pulku
dc1 = zeros(ceil(size(d,1)/podil),els);
for j = 1:els %musim decimovat kazdou elektrodu zvlast
    dc1(:,j) = decimate(d(:,j),podil); %na 500 Hz z 8000 Hz
end
clear d;

%ostatni promenne krome d nerozdeluju
if exist('t', 'var')
        t = downsample(t,podil); %èas v sekundách
end
if exist('tabs', 'var')
    tabs = downsample(tabs,podil); %èas formatu    28-Jan-2014 11:35:45.000
end
if exist('timerel', 'var')
    timerel = downsample(timerel,podil); % co to je ???
end

%druha pulka souboru
load([ dir soubor '.mat'],'d'); %maly soubor, nactu jen d
d(1:dpul,:)=[]; %smazu prvni pulku souboru

dc2 = zeros(ceil(size(d,1)/podil),els);
for j = 1:els %musim decimovat kazdou elektrodu zvlast
    dc2(:,j) = decimate(d(:,j),podil); %na 500 Hz z 8000 Hz
end

clear d;

%spojim obe pulky souboru
d = vertcat(dc1,dc2); 
clear('dc1','dc2','j','els','podil','dpul'); %vymazu promenne, ktere uz nepotrebuju

%data ulozim do noveho souboru
save([dir soubor '_' frek 'hz']); %,'H','d','t','tabs','timerel' - 12.2.2015 - nekdy chybi tabs, takze musim ulozit vsechno co je

