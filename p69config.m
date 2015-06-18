function [ p69 ] = p69config(  )
%P73CONFIG vrati konfiguraci pro p73
%   vraci struct
p69.els = [12 20 27 37 47 56 64 75 83 90 99 109 119 125]; %p69 - hranice elektrod konce
p69.channels = 1:125;
p69.name = 'p69';
p69.cas = [-100 1000]; %AEDist
end

