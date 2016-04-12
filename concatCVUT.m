%spoji dva rozdelene soubory
% CVUT data
% nefunguje to na notebooku protoze to potrebuje prilis mnoho pameti - 18.6.2015
% funguje na Amygdale - 23.6.2015
clear all;
spojit = {
        'f:\eeg\motol\pacienti\p97 Novak VT13\VT13_2016-02-11_09-20_001.mat',...
        'f:\eeg\motol\pacienti\p97 Novak VT13\VT13_2016-02-11_10-20_002.mat'
        
        };

for j = 1:numel(spojit)
    disp(spojit{j});
    load(spojit{j});
    disp(['delka ' num2str(size(tabs,1))]);
    if j == 1 
        tabs0 = tabs;
        clear tabs;
        if exist('header','var') 
            header0 = header;
            clear header;
        end
        if exist('mults','var') %mults obsahuje jednu hodnotu pro kazdy kanal
            mults0 = mults;
            clear mults;
        end
        if exist('fs','var') 
            fs0 = fs;
            clear fs;
        end        
        d0 = d;
        clear d;
       
    else
        tabs0 = [tabs0; tabs]; %#ok<AGROW>
        clear tabs;
        d0 = [d0; d]; %#ok<AGROW>
        clear d;
        if exist('header','var') 
            header0.length = header0.length + header.length;
            header0.records = header0.records + header.records;
        end         
    end
    
end
d = d0;
clear d0;
tabs = tabs0;
clear tabs0;
if exist('header','var') 
    header = header0;
    clear header0;
end
if exist('mults','var') 
    mults = mults0;
    clear mults0;
end
if exist('fs','var') 
    fs = fs0;
    clear fs0;
end
disp(['vysledna delka ' num2str(size(tabs,1))]);
dot = strfind(spojit{1},'.');
disp(['ukladam ' spojit{1}(1:dot(1)-1) '_concat.mat']);
save([ spojit{1}(1:dot(1)-1) '_concat.mat'], '-regexp', '^(?!(spojit|j|dot|OBJ)$).','-v7.3');
clear spojit j dot;

