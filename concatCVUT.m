%spoji dva rozdelene soubory
% CVUT data
% nefunguje to na notebooku protoze to potrebuje prilis mnoho pameti - 18.6.2015
clear all;
spojit = {
        'w:\eeg\motol\pacienti\P79 Pluhackova VT8\VT8_2015-04-09_09-46_001.mat',...
        'w:\eeg\motol\pacienti\P79 Pluhackova VT8\VT8_2015-04-09_10-46_002.mat' ...
        };

for j = 1:numel(spojit)
    disp(spojit{j});
    load(spojit{j});
    if j == 1 
        d0 = d;
        clear d;
        tabs0 = tabs;
        clear tabs;
        header0 = header;
        clear header;
    else
        d0 = [d0; d]; %#ok<AGROW>
        clear d;
        tabs0 = [tabs0; tabs]; %#ok<AGROW>
        clear tabs;
    end
    
end
d = d0;
clear d0;
tabs = tabs0;
clear tabs0;
header = header0;
clear header0;
dot = strfind(spojit{1},'.');

save( [ spojit{1}(1:dot(1)-1) '_concat'],header,d,tabs);


