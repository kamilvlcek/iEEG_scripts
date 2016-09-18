function [alloego] = alloego_data(Events,smazat)
%ALLOEGO_DATA vytvori a vraci strukturu s udalostmi experimentu AlloEgo v BVA. 
% vychází ze struktury Events of Lukase 

%tabulka odpovedi pacienta
%musim profiltrovat stlacene klavesy a setoff goal aby to byla sekvence c off f g
%zatim to udelam jen pole smazat
if exist('smazat','var') && isstruct(smazat) %pole obsahuje indexy ke smazani , smazat = struct;
    if isfield(smazat,'f')
        Events.f.timestamps(smazat.f) = []; 
        disp(['smazano f BVA_frames:' mat2str(Events.f.BVA_frames(smazat.f))]); %napr. p83 smazat.f = [6 12 15];
    end
    if isfield(smazat,'g')
        Events.g.timestamps(smazat.g) = [];
        disp(['smazano g BVA_frames:' mat2str(Events.g.BVA_frames(smazat.g))]); %napr. p83 smazat.g = 11;
    end
end
    


assert(numel(Events.set_off_goal.timestamps) == numel(Events.f.timestamps),'ruzne delky setoffgoal a f');
assert(numel(Events.set_off_goal.timestamps) == numel(Events.c.timestamps),'ruzne delky setoffgoal a c');

data = zeros(numel(Events.set_off_goal.timestamps),9);
data(:,3)=ones(numel(Events.set_off_goal.timestamps),1); %vsechno spravne
data(:,4)=(Events.f.timestamps-Events.c.timestamps)*3600*24; %reakcni cas - rozdil mezi cas f a c
data(:,8)=Events.set_off_goal.timestamps; %casy zacatku pohybu
data(:,9)=Events.f.timestamps; %casy stlaceni f - zadne konce pohybu nejsou

assert(sum(data(:,8) < data(:,9)) == size(data,1),'stimuli nejsou vzdy driv nez odpovedi');
%nazvy sloupcu tabulky
sloupce = {};
sloupce.soubor=1;
sloupce.klavesa=2;
sloupce.spravne=3;
sloupce.rt = 4;
sloupce.opakovani=5;
sloupce.zpetnavazba=6;
sloupce.kategorie=7;
sloupce.ts_podnet=8;
sloupce.ts_odpoved=9;

%retezcove hodnoty kodu klavesa a faktoru testu, musi souhlasit s kody v AEdistData.php
klavesa = cell(1,2);
klavesa(1,:)={'None' 0};

podminka = cell(1,2);
podminka(1,:)={'vsechny' 0};

alloego = struct('data',data,'sloupce',sloupce);
alloego.strings.klavesa = klavesa;
alloego.strings.podminka = podminka; %kategorie, aby nazev byl stejny jako u PPA

%timestampy zacatku a konce dat z testu
alloego.interval = [Events.e.timestamps(1) max(Events.e.timestamps(end),Events.g.timestamps(end))];
disp( [ 'inverval: ' num2str(datestr(alloego.interval(1) ,'dd-mmm-yyyy HH:MM:SS.FFF')) ' - ' num2str(datestr(alloego.interval(2) ,'dd-mmm-yyyy HH:MM:SS.FFF')) ]);
end

