function [menrot] = menrot_data(pacientid,U1,U2,tabs)
%AEDIST_DATA vytvori a vraci strukturu s udalostmi experimentu AEdist
% pacientid - id pacienta, napriklad p85, podle pojmenovani vystupni tabulky AEdistData.php.
% U1 a U2 vystup z udalosti2() - casy synchropulsu k podnedu a odpovedi

dir = 'd:\prace\homolka\epileptici EEG\vysledky\menrot\';
%tabulka odpovedi pacienta
data = load([dir pacientid '_menrot.txt']);

%musim ze dvou sloupcu [0 1] udelat jeden sloupec [0 1 2 3]
data(:,9) = data(:,7)*2+data(:,8); %kategorie vy/znacka a 2D/3D
data(:,7:8) = [];

if size(U1,1) ~= size(data,1) || size(U2,1) ~= size(data,1)
    disp(['ruzne delky dat (' num2str(size(data,1)) ') a  udalosti (' num2str(size(U1,1)) ',' num2str(size(U2,1)) '), nelze zpracovat!']);
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
sloupce.kategorie=7;
sloupce.ts_podnet=8;
sloupce.ts_odpoved=9;

%retezcove hodnoty kodu klavesa a faktoru testu, musi souhlasit s kody v AEdistData.php
klavesa = cell(3,2);
klavesa(1,:)={'None' -1};
klavesa(2,:)={'left' 0};
klavesa(3,:)={'right' 1};

podminka = cell(3,2);
podminka(1,:)={'vy-2D' 0};
podminka(2,:)={'vy-3D' 1};
podminka(3,:)={'znacka-2D' 2};
podminka(4,:)={'znacka-3D' 3};

menrot = struct('data',data,'sloupce',sloupce);
menrot.strings.klavesa = klavesa;
menrot.strings.podminka = podminka; %kategorie, aby nazev byl stejny jako u PPA

%timestampy zacatku a konce dat z testu
menrot.interval = [tabs(1) tabs(end)];
end

