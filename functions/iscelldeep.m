function [ output] = iscelldeep( val )
%ISCELLDEEP vrati true pokud je val cell, a nejaky z jejich prvku je taky cell
%   kvuli stat_kats = {[0 1 2 3],{[0 1],[2 3]},{[0 2],[1 3]}};

if iscell(val)
    output = false;
    for j = 1:numel(val)        
        if iscell(val{j})
            output = true;
            break;
        end        
    end
else
    output= false;
end


end

