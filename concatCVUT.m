function [delka] = concatCVUT(adresar,spojit)
%spoji dva rozdelene soubory
% CVUT data
% nefunguje to na notebooku protoze to potrebuje prilis mnoho pameti - 18.6.2015
% funguje na Amygdale - 23.6.2015

%spojit = {
%        'f:\eeg\motol\pacienti\p97 Novak VT13\VT13_2016-02-11_09-20_001.mat',...
%        'f:\eeg\motol\pacienti\p97 Novak VT13\VT13_2016-02-11_10-20_002.mat'
%        
%        };
if ~exist('spojit','var')  
    %mam udelat vsechny soubory v adresari    
    files = dir(fullfile(adresar, '*.mat'));
    spojit = cell(numel(files),1);
    for f = 1:numel(files)
        spojit{f} = [ files(f).name];        
    end
end
    
for j = 1:numel(spojit)
    disp([ num2str(j) '/' num2str(numel(spojit)) ': ' adresar spojit{j}]);
    load([adresar spojit{j}]);
    disp(['delka ' num2str(size(tabs,1))]); %#ok<NODEF>
    if j == 1 
        tabs0 = tabs;
        clear tabs;
        if exist('header','var') 
            header0 = header; %#ok<NODEF>
            clear header;
        end
        if exist('mults','var') %mults obsahuje jednu hodnotu pro kazdy kanal
            mults0 = mults; %#ok<NODEF>
            clear mults;
        end
        if exist('fs','var') 
            fs0 = fs; %#ok<NODEF>
            clear fs;
        end  
        if exist('evts','var')
            evts0 = evts; %#ok<NODEF>
            clear evts;
        end
        d0 = d; %#ok<NODEF>
        clear d;
       
    else
        disp(['rozdil ' num2str((tabs(1)-tabs0(end))*24*3600) ' sekund']);
        if (tabs(1)-tabs0(end))*24*3600 >= 1 %rozdil jedne vteriny je velmi zvlastni, nejspis se soubory nemaji spojit            
            m=input('Do you want to continue, Y/N [N]:','s');
            if m~='y',   break; end
        end
        tabs0 = [tabs0; tabs]; %#ok<AGROW>        
        clear tabs;
        d0 = [d0; d]; %#ok<NODEF,AGROW>
        clear d;
        if exist('evts','var')
            evts0 = [evts0, evts]; %#ok<AGROW,NODEF>
        end
        if exist('header','var') 
            header0.length = header0.length + header.length; %#ok<NODEF>
            header0.records = header0.records + header.records;
        end         
    end
    
end
d = d0;
clear d0;
tabs = tabs0;
clear tabs0;
if exist('header','var') 
    header = header0; %#ok<NASGU>
    clear header0;
end
if exist('mults','var') 
    mults = mults0; %#ok<NASGU>
    clear mults0;
end
if exist('fs','var') 
    fs = fs0; %#ok<NASGU>
    clear fs0;
end
if exist('evts','var')
    evts = remdupstruct(evts0); %#ok<NASGU>
    clear evts0;
end
delka = num2str(size(tabs,1));
disp(['vysledna delka ' delka]);
dot = strfind(spojit{1},'.');
disp(['ukladam ' adresar '\' spojit{1}(1:dot(1)-1) '_concat.mat']);

save([ adresar '\' spojit{1}(1:dot(1)-1) '_concat.mat'], '-regexp', '^(?!(spojit|j|dot|OBJ,f|delka|filename|adresar)$).','-v7.3');

end


