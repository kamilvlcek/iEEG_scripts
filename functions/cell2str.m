function [str] = cell2str(cellvar,filename)
%funkce pro prevedeni jednorozmerneho cell array do string. Vola rekurzi dalsi urovne
%filename : 1 - chci vysledke kompatibilni s filename
%kamil 6.3.2018
if ~exist('filename','var'), filename = 0; end
if  filename %returns eg. (tt;6-0,7-0;6-0,7-1;6-1,7-0;6-1,7-1)
    lbrace = '(';
    rbrace = ')';
    lbracket = '';
    rbracket = '';
    comma = ';';
else % returnse.g. {tt|6-0,7-0|6-0,7-1|6-1,7-0|6-1,7-1}
    lbrace = '{';
    rbrace = '}';
    lbracket = ''; %[ - when using mat2str, no bracket necesary, it uses itself
    rbracket = ''; %]
    comma =  '|';
end

if iscell(cellvar)
    str = lbrace;
    for k = 1:numel(cellvar)
        if iscell(cellvar{k})
            str = [str cell2str(cellvar{k})]; %#ok<AGROW>
        elseif ischar(cellvar{k})
            str = [str cellvar{k}]; %#ok<AGROW>            
        else
            str = [str lbracket replace(mat2str(cellvar{k}'), {' ', '[', ']',';'}, {'-', '', '',','}) rbracket]; %#ok<AGROW> % 3.7.2023 - new form for cellvar{k} being a matrix 
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

