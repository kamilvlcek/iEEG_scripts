function [aedist] = aedist_data(pacientid,U1,U2,tabs)
%AEDIST_DATA vytvori a vraci strukturu s udalostmi experimentu AEdist
% pacientid - id pacienta, napriklad p85, podle pojmenovani vystupni tabulky AEdistData.php.
% U1 a U2 vystup z udalosti2() - casy synchropulsu k podnedu a odpovedi

dir = 'd:\prace\homolka\epileptici EEG\vysledky\AEDist\';
%tabulka odpovedi pacienta
data = load([dir pacientid '_aedist.txt']);
if size(U1,1) ~= size(data,1) || size(U2,1) ~= size(data,1)
    disp('ruzne delky dat a  udalosti, nelze zpracovat!');
    return;
end
data(:,8)=U1(:,2);
data(:,9)=U2(:,2);

%nazvy sloupcu tabulky
sloupce = {};
sloupce.soubor=1;
sloupce.klavesa=2;
sloupce.spravne=3;
sloupce.rt = 4;
sloupce.opakovani=5;
sloupce.zpetnavazba=6;
sloupce.podminka=7;
sloupce.ts_podnet=8;
sloupce.ts_odpoved=9;

%retezcove hodnoty kodu klavesa a faktoru testu, musi souhlasit s kody v AEdistData.php
klavesa = cell(3,2);
klavesa(1,:)={'None' -1};
klavesa(2,:)={'left' 0};
klavesa(3,:)={'right' 1};

podminka = cell(3,2);
podminka(1,:)={'cervena' 0};
podminka(2,:)={'vy' 1};
podminka(3,:)={'znacka' 2};

aedist = struct('data',data,'sloupce',sloupce);
aedist.strings.klavesa = klavesa;
aedist.strings.podminka = podminka;

%timestampy zacatku a konce dat z testu
aedist.interval = [tabs(1) tabs(end)];
end

