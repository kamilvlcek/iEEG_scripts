function [ casstr ] = cas(t, index )
%CAS Summary of this function goes here
%   Detailed explanation goes here

casstr = datestr(t(index)/3600/24,'HH:MM:SS.FFF');

end

