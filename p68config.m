function [ p68 ] = p68config(  )
%P68CONFIG vrati konfiguraci pro p68
%   vraci struct
p68.els = [ 10 18 28  41 54 64 72 82 91 100 114 125]; %p68 - hranice elektrod konce
p68.channels = 1:125;
p68.name = 'p68';
p68.cas = [-200 1000]; %AEDist
end

