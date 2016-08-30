function [fs, delka] = zmensidata(filename,podil)
% ZMENSIDATA funkce na zmenseni vzorkovaci frekvence dat z Motola
% prevedeno na funkce 30.8.2016
% vraci vyslednou frekvenci a delku dat

%prvni pulka souboru - s celym najednou se spatne pracuje
load(filename,'d','tabs','fs'); 

assert( rem(fs,podil) ==0, 'vysledna frekvence musi byt cele cislo'); %#ok<NODEF>
frek = num2str(fs/podil); %vysledna frekvence
fs = fs / podil;

disp('first half of d ... ');
dpul = floor(size(d,1)/(podil*2)) * podil; %#ok<NODEF> % aby delitelne podilem - pulka zaokrouhlena dolu
els = size(d,2);
d(dpul+1:end,:)=[]; %smazu druhou pulku

dc1 = zeros(ceil(size(d,1)/podil),els); % d uz obsahuje jen svoji prvni pulku a delka je delitelna podilem
for j = 1:els %musim decimovat kazdou elektrodu zvlast
    dc1(:,j) = decimate(d(:,j),podil); %na 500 Hz z 8000 Hz
end
clear d;

%ostatni promenne krome d nerozdeluju
if exist('tabs', 'var')
    tabs = downsample(tabs,podil); %#ok<NODEF> %èas formatu    28-Jan-2014 11:35:45.000
end

%druha pulka souboru
disp('second half of d ...');
load(filename,'d'); %maly soubor, nactu jen d
d(1:dpul,:)=[]; %smazu prvni pulku souboru

dc2 = zeros(ceil(size(d,1)/podil),els);
for j = 1:els %musim decimovat kazdou elektrodu zvlast
    dc2(:,j) = decimate(d(:,j),podil); %na 500 Hz z 8000 Hz
end

clear d;

%spojim obe pulky souboru
d = vertcat(dc1,dc2);  %#ok<NASGU>
delka = numel(tabs);


%data ulozim do noveho souboru
[pathstr,fname,ext] = fileparts(filename);  %ext bude .mat
newfilename = [pathstr '\' fname '_' frek 'hz' ext];
save(newfilename, 'd','tabs','fs','-v7.3'); 
disp(['saved as ' newfilename]);
end
