function [ output] = iscelldeep( val )
%ISCELLDEEP vrati true pokud je val cell, a nejaky z jejich prvku je taky cell
%   kvuli stat_kats = {[0 1 2 3],{[0 1],[2 3]},{[0 2],[1 3]}};

if iscell(val)
    for j = 1:numel(val)
        output = false;
        if iscell(val{j})
            output = true;
        end        
    end
else
    output= false;
end


end

