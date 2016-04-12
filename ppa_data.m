function [ppa] = ppa_data(pacientid,U1,U2,tabs,eegfile)
%PPA_DATA vytvori a vraci strukturu s udalostmi experimentu PPA
% pacientid - id pacienta, napriklad p85, podle pojmenovani vystupni tabulky AEdistData.php.
% U1 a U2 vystup z udalosti2() - casy synchropulsu k podnedu a odpovedi
% eegfile je nazev puvodniho souboru s eeg daty

dir = 'd:\prace\homolka\epileptici EEG\vysledky\PPalokalizer\';
%tabulka odpovedi pacienta
data = load([dir pacientid '_ppa.txt']);
if size(U1,1) ~= size(data,1) || size(U2,1) ~= size(data,1)
    disp('ruzne delky dat a  udalosti, nelze zpracovat!');
    return;
end
data(:,9)=U1(:,2);
data(:,10)=U2(:,2);

%nazvy sloupcu tabulky
sloupce = {};
sloupce.soubor=1;
sloupce.klavesa=2;
sloupce.spravne=3;
sloupce.rt = 4;
sloupce.opakovani_obrazku=5;
sloupce.cislo_obrazku=6;
sloupce.pauza=7;
sloupce.kategorie=8;
sloupce.ts_podnet=9;
sloupce.ts_odpoved=10;

%retezcove hodnoty kodu klavesa a faktoru testu, musi souhlasit s kody v AEdistData.php
klavesa = cell(2,2);
klavesa(1,:)={'None' -1};
klavesa(2,:)={'space' 1};

podminka = cell(4,2);
podminka(1,:)={'Ovoce' 0};
podminka(2,:)={'Scene' 1};
podminka(3,:)={'Face' 2};
podminka(4,:)={'Object' 3};

ppa = struct('data',data,'sloupce',sloupce);
ppa.strings.klavesa = klavesa;
ppa.strings.podminka = podminka;

%timestampy zacatku a konce dat z testu
ppa.interval = [tabs(1) tabs(end)];
%vypis pro kontrolu intevalu
disp(['PPA data od ' datestr(tabs(1),'dd-mmm-yyyy HH:MM:SS.FFF') ' do ' datestr(tabs(end),'dd-mmm-yyyy HH:MM:SS.FFF')]);

ppa.eegfile = eegfile;
ppa.pacientid = pacientid;
end

