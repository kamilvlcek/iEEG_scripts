function [ p71 ] = p71config(  )
%P71CONFIG vrati konfiguraci pro p73
%   vraci struct
p71.els = [ 16 32 48  64 77 89 105 115 125]; %hranice elektrod pro p71 - konce
p71.channels = 1:125;
p71.name = 'p71';
end

