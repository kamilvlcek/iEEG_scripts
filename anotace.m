function [  ] = anotace( H, tabs, zac)
%ANOTACE vypise vsechny anotace v zaznamu spolu s jejich casem
%   H, zacatek 1 nebo vice

if nargin < 3, zac=1; end %kdyz neuvedu zacatek, jede se od prvni anotace
disp(['celkova delka v sekundach ' num2str(H.records) ',' ...
    'od ' datestr(tabs(1),'dd-mmm-yyyy HH:MM:SS.FFF') ', do '  datestr(tabs(end),'dd-mmm-yyyy HH:MM:SS.FFF') ]);


for a = zac: size(H.annotation.starttime)
    
   disp([ num2str(a) ' ' strtrim(H.annotation.event{a}) ' ... ' num2str( H.annotation.starttime(a)) ...
       ' s ... ' num2str( H.annotation.starttime(a)*H.samplerate(1))]);
   cas = H.annotation.starttime(a);
   %cas0 = iff(a == zac,H.annotation.starttime(1),H.annotation.starttime(a-1));
   if a==zac, figure('Name','Anotace'); end
   plot([cas cas],[0 1]); %vykreslim cas anotaci do grafu
   hold on;
end


end

