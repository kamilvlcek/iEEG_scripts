function [str] = cell2str(cellvar,filename)
%funkce pro prevedeni jednorozmerneho cell array do string. Vola rekurzi dalsi urovne
%filename : 1 - chci vysledke kompatibilni s filename
%kamil 6.3.2018
if ~exist('filename','var'), filename = 0; end
if filename
    lbrace = '(';
    rbrace = ')';
    lbracket = '';
    rbracket = '';
    comma = '-';
else
    lbrace = '{';
    rbrace = '}';
    lbracket = '[';
    rbracket = ']';
    comma =  ',';
end

if iscell(cellvar)
    str = lbrace;
    for k = 1:numel(cellvar)
        if iscell(cellvar{k})
            str = [str cell2str(cellvar{k})]; %#ok<AGROW>
        else
            str = [str lbracket num2str(cellvar{k}) rbracket]; %#ok<AGROW>
        end 
        if k < numel(cellvar)
            str = [str comma]; %#ok<AGROW>
        end
    end
    str = [str rbrace];
else
    str = [lbracket num2str(cellvar) rbracket];
end
end

