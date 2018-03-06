function [str] = cell2str(cellvar)
%funkce pro prevedeni jednorozmerneho cell array do string. Vola rekurzi dalsi urovne
%kamil 6.3.2018
if iscell(cellvar)
    str = '{';
    for k = 1:numel(cellvar)
        if iscell(cellvar{k})
            str = [str cell2str(cellvar{k})]; %#ok<AGROW>
        else
            str = [str '[' num2str(cellvar{k}) ']']; %#ok<AGROW>
        end 
        if k < numel(cellvar)
            str = [str ',']; %#ok<AGROW>
        end
    end
    str = [str '}'];
else
    str = ['[' num2str(cellvar) ']'];
end
end

