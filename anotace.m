function [  ] = anotace( H, tabs, zac)
%ANOTACE vypise vsechny anotace v zaznamu spolu s jejich casem
%   H, zacatek 1 nebo vice

disp(['celkova delka v sekundach ' num2str(H.records) ',' ...
    'od ' datestr(tabs(1),'dd-mmm-yyyy HH:MM:SS.FFF') ', do '  datestr(tabs(end),'dd-mmm-yyyy HH:MM:SS.FFF') ]);
for a = zac: size(H.annotation.starttime)
   disp([ num2str(a) ' ' strtrim(H.annotation.event{a}) ' ... ' num2str( H.annotation.starttime(a)) ...
       ' s ... ' num2str( H.annotation.starttime(a)*H.samplerate(1))]);
end
end

