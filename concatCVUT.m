function [files] = concatCVUT(adresar,spojit,testrozdil)
%spoji dva rozdelene soubory
% CVUT data
% nefunguje to na notebooku protoze to potrebuje prilis mnoho pameti - 18.6.2015
% funguje na Amygdale - 23.6.2015

%spojit = {
%        'f:\eeg\motol\pacienti\p97 Novak VT13\VT13_2016-02-11_09-20_001.mat',...
%        'f:\eeg\motol\pacienti\p97 Novak VT13\VT13_2016-02-11_10-20_002.mat'
%        
%        };
if ~exist('testrozdil','var'), testrozdil = 1; end
if ~exist('spojit','var') || isempty(spojit)
    if numel(adresar) == 1
        %mam udelat vsechny soubory v adresari    
        files = dir(fullfile(adresar, '*.mat'));
        spojit = cell(numel(files),1);
        for f = 1:numel(files)
            spojit{f} = [ files(f).name];        
        end
    else
        spojit = adresar; %v adresar jsou nejspis soubory k spojeni vcetne plne cesty
        adresar = '';
    end
end
spojeno = 0;    %kolik souboru jsem spojil
j0 = 1; %index prvniho souboru, ktere spojuji
files = cell(0,4); %tam budu uklada jmena vytvorenych souboru
for j = 1:numel(spojit)
    disp([ num2str(j) '/' num2str(numel(spojit)) ': ' adresar spojit{j}]);
    load([adresar spojit{j}]);
    if ~exist('mults','var'), mults = []; end
    if ~exist('evts','var'), evts = []; end
    if ~exist('header','var'), header = []; end
    
    disp(['delka ' num2str(size(tabs,1)) ' vzorku']); %#ok<NODEF> 
    if j> j0 %pokud se nejedna o prvni soubor
        rozdil_sec = (tabs(1)-tabs0(end))*24*3600; %rozdil mezi koncem jednoho a zacatkem druheho
        disp(['rozdil ' num2str(rozdil_sec) ' sekund']); 
        if testrozdil && (rozdil_sec >= (tabs(3)-tabs(1))*24*3600 || rozdil_sec < 0) %rozdil je dvojnasobek normalniho rozdilu, nejspis se soubory nemaji spojit                        
            [filename,delka] = ulozdata(j0,spojeno,adresar,spojit,d0,tabs0,fs0,header0,mults0,evts0); 
            files = [files; {filename, spojeno, delka, rozdil_sec}]; %#ok<AGROW>
            j0 = j;      %aktualni soubor chci zpracovat jako prvni - pocitani zacinam od zacatku  
        end
    end
    if j0 == j %pokud jde o prvni soubor
        tabs0 = tabs;
        clear tabs;
        if exist('header','var') 
            header0 = header; 
            clear header;
        else
            header0 = [];
        end
        if exist('mults','var') %mults obsahuje jednu hodnotu pro kazdy kanal
            mults0 = mults; 
            clear mults;
        else
            mults0 = [];
        end
        if exist('fs','var') 
            disp(['fs=' num2str(fs)]); 
            fs0 = fs; 
            clear fs;
        else
            fs0 = [];
        end  
        if exist('evts','var')
            evts0 = evts; 
            clear evts;
        else
            evts0 = [];
        end
        d0 = d;  %#ok<NODEF>
        clear d;
        spojeno = 1; 
    else %nasledujici soubory
        
        tabs0 = [tabs0; tabs]; %#ok<AGROW>        
        clear tabs;
        d0 = [d0; d]; %#ok<NODEF,AGROW>
        clear d;
        if exist('evts','var')
            evts0 = [evts0, evts]; %#ok<AGROW>
        end
        if exist('header','var') 
            if isfield(header,'length') 
                header0.length = header0.length + header.length;
            end
            if isfield(header,'records') 
                header0.records = header0.records + header.records;
            end
        end 
        spojeno = spojeno +1; 
       
    end
    
end
d = d0;
clear d0;
tabs = tabs0;
clear tabs0;
if exist('header0','var') 
    header = header0; 
    clear header0;
end
if exist('mults0','var') 
    mults = mults0; 
    clear mults0;
end
if exist('fs0','var') 
    fs = fs0; 
    clear fs0;
end
if exist('evts0','var') 
    evts = remdupstruct(evts0);
    clear evts0;
end

[filename,delka] = ulozdata(j0,spojeno,adresar,spojit,d,tabs,fs,header,mults,evts);
files = [files; {filename, spojeno, delka, []}];

end 

%local function
function [filename,delka] = ulozdata(j0,spojeno,adresar,spojit,d,tabs,fs,header,mults,evts)     %#ok<INUSL,INUSD>
    delka = size(tabs,1);
    disp(['vysledna delka ' num2str(delka) ' vzorku']);
    if isempty(mults), clear mults; end
    if isempty(header), clear header; end
    if isempty(evts), clear evts; end
    %if spojeno > 1 - ulozim i jen jeden soubor, aby to bylo prehledne
        dot = strfind(spojit{j0},'.');
        filename = [adresar  spojit{j0}(1:dot(1)-1) '_' num2str(spojeno) '_concat.mat'];
        disp(['ukladam ze ' num2str(spojeno) ' souboru: ' filename]);
        save(filename, '-regexp', '^(?!(spojit|j|dot|OBJ,f|delka|filename|adresar)$).','-v7.3');
    %else 
    %    disp('neukladam - pouze jeden soubor');
    %end    
end



