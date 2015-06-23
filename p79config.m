function [ p79 ] = p79config(  )
%P68CONFIG vrati konfiguraci pro p79
%   vraci struct
p79.els = [ 11 21 29 39 49 55 63 64 75 88 99 107 119 126]; %p79 - hranice elektrod =konce
p79.channels = [1:63 65:126];
p79.name = 'p79';
p79.AEDist.cas = [-200 1000]; %??? kontrola AEDist
p79.lpt = 64; %synchronizacni puls
end

