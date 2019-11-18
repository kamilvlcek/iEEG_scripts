function [ output] = iscelldeep( val )
%ISCELLDEEP vrati true pokud je val cell, a nejaky z jejich prvku je taky cell a zadny z prvku neni jedno cislo
%   kvuli stat_kats = {[0 1 2 3],{[0 1],[2 3]},{[0 2],[1 3]}};

if iscell(val)
    output = true; %false pouze pokud nektery z prvku neni dalsi pole. 
    for j = 1:numel(val)        
        if iscell(val{j})
            output = true;
            break;
        end        
        if numel(val{j}) < 2
          output = false; %pokud kterykoliv z prvku neni dalsi pole
          break;
        end
    end
else
    output= false;
end


end

