function [  ] = anotace( H, zac)
%ANOTACE Summary of this function goes here
%   Detailed explanation goes here

disp(['celkova delka v sekundach ' num2str(H.records)]);
for a = zac: size(H.annotation.starttime)
   disp([ num2str(a) ' ' strtrim(H.annotation.event{a}) ' ... ' num2str( H.annotation.starttime(a)) ...
       ' s ... ' num2str( H.annotation.starttime(a)*H.samplerate(1))]);
end
end

