function [newnames, fs, delka] = zmensidata(filenames,podil,yratio,fsforce)
% ZMENSIDATA funkce na zmenseni vzorkovaci frekvence dat z Motola
% prevedeno na funkce 30.8.2016
% vraci vyslednou frekvenci a delku dat
% filename - pokud chci zpracovat vic souboru, filename je adresar (s \ na konci ) - pak zpracuje vsechny soubory 
% nebo filename = adresar + zacatek souboru (cele bez \ na konci) = maska - pak zpracuju jen vyber souboru
newnames = {};
if ~exist('fsforce','var'), fsforce = []; end %vynucna vystupni fs. Na co?
if ~iscell(filenames), filenames = {filenames}; end %pokud jen jeden soubor, nemusi byt cell array
for ff = 1:numel(filenames) %muzu mit celou serii adresaru na zpracovani se stejnymi parametry - 2019/07
    filename = filenames{ff};    
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
            disp(['filename ' num2str(f) '/' num2str(numel(files)) '(' num2str(ff) '/' num2str(numel(filenames)) '): ' filename]);
            [newfilename,fs,delka]=zmensidata({filename},podil,yratio,fsforce);
            newnames = [newnames; newfilename]; %#ok<AGROW>
        end

    else
        %zpracovavam jeden soubor
        fprintf('processing file ... %d/%d: %s\n', ff, numel(filenames), filename);
        load(filename);   %#ok<LOAD>
        if size(d,1)<2
            disp(['data are too short: ' num2str(size(d,1)) ' samples' ]);
            delka = numel(tabs);
            %newnames = [newnames; filename]; %#ok<AGROW> %jednovzorkovy zaznam vubec nechci dal zpracovavat
            continue;
        end
        if exist('mults', 'var')
            load(filename,'mults'); %pokud existuji, roznasobim to mults - decimate pracuje jen s double
            d = bsxfun(@times,double(d), mults); %rovnou to roznasobim mults
        end
        if exist('yratio', 'var') %cim roznasobit data, aby byly stejne velike jako v systemu nickone - 0.1mV?
            d = d*yratio; %novy system Quantum ma data 100x vetsi - uV?
        end
        if exist('fsforce','var') && ~isempty(fsforce)
            fs = fsforce;
        end

        assert( rem(fs,podil) ==0, 'the resulting frequency should be an integral number'); 
        frek = num2str(fs/podil); %vysledna frekvence
        fs = fs / podil;

        els = size(d,2);       

        dc1 = zeros(ceil(size(d,1)/podil),els); % d uz obsahuje jen svoji prvni pulku a delka je delitelna podilem
        for j = 1:els %musim decimovat kazdou elektrodu zvlast
            dc1(:,j) = decimate(d(:,j),podil); %na 500 Hz z 8000 Hz
        end
        
        %ostatni promenne krome d nerozdeluju
        if exist('tabs', 'var')
            tabs = downsample(tabs,podil);  %èas formatu    28-Jan-2014 11:35:45.000
        end

        delka = numel(tabs);
        d = dc1;
        clear dc1;

        %data ulozim do noveho souboru
        [pathstr,fname,ext] = fileparts(filename);  %ext bude .mat
        newfilename = [pathstr '\' fname '_' frek 'hz' ext];    
        save(newfilename, ... 
            '-regexp', '^(?!(mults|dc1|dc2|delka|dpul|els|ext|frek|fsforce|j|fname|name|newfilename|pathstr|podil|ff|filenames|newnames)$).','-v7.3');
        disp(['saved as ' newfilename]);
        newnames = [newnames; newfilename]; %#ok<AGROW>
    end
end
end
