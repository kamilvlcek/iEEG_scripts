function [fs, delka] = zmensidata(filename,podil,fsforce)
% ZMENSIDATA funkce na zmenseni vzorkovaci frekvence dat z Motola
% prevedeno na funkce 30.8.2016
% vraci vyslednou frekvenci a delku dat
% filename - pokud chci zpracovat vic souboru, filename je adresar (s \ na konci ) - pak zpracuje vsechny soubory 
% nebo filename = adresar + zacatek souboru (cele bez \ na konci) = maska - pak zpracuju jen vyber souboru

if ~exist('fsforce','var'), fsforce = [];end %vynucna vystupni fs. Na co?

[~,name,ext] = fileparts(filename);
if  isempty(ext)
    %asi se jedna o adresar - zpracuju postupne vsechny soubory
    adresar = filename;  
    if isempty(name) %filename koncil \ a jedna se tedy o adresar          
        files = dir(fullfile(adresar, '*.mat'));
    else %filename ma na konci zacatek = masku jmena
        files = dir([adresar '*.mat']);
    end
    for f = 1:numel(files)
        filename = [ files(f).folder filesep files(f).name];
        disp(['filename ' num2str(f) '/' num2str(numel(files)) ': ' filename]);
        [fs,delka]=zmensidata(filename,podil,fsforce);
    end
    
else
    %zpracovavam jeden soubor
    %prvni pulka souboru - s celym najednou se spatne pracuje
    load(filename);  
    if size(d,1)<2
        disp(['zaznam prilis kratky: ' num2str(size(d,1)) ' vzorku' ]);
        delka = numel(tabs);
        return;
    end
    if exist('mults', 'var')
        load(filename,'mults'); %pokud existuji, roznasobim to mults - decimate pracuje jen s double
        d = bsxfun(@times,double(d), mults); %#ok<NODEF> %rovnou to roznasobim mults
    end
    if exist('fsforce','var') && ~isempty(fsforce)
        fs = fsforce;
    end

    assert( rem(fs,podil) ==0, 'vysledna frekvence musi byt cele cislo'); 
    frek = num2str(fs/podil); %vysledna frekvence
    fs = fs / podil;

    disp('first half of d ... ');
    dpul = floor(size(d,1)/(podil*2)) * podil;  % aby delitelne podilem - pulka zaokrouhlena dolu
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
    if exist('mults', 'var') 
        d = bsxfun(@times,double(d), mults); %rovnou to roznasobim mults - decimate pracuje jen s double          
    elseif ~isfloat(d)
        disp('d neni float type a neexistuje mults');
        delka = numel(tabs);
        return;
    end
    
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
    save(newfilename, '-regexp', '^(?!(mults|dc1|dc2|delka|dpul|els|ext|frek|fsforce|j|fname|name|newfilename|pathstr|podil)$).','-v7.3');
    disp(['saved as ' newfilename]);

end
end
