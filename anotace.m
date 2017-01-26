function [  ] = anotace( H, tabs, evts, zac)
%ANOTACE vypise vsechny anotace v zaznamu spolu s jejich casem
%   H, zacatek 1 nebo vice

if ~exist('zac','var'), zac=1; end %kdyz neuvedu zacatek, jede se od prvni anotace

if exist('H','var') && ~isempty(H)
    disp(['celkova delka v sekundach ' num2str(H.records) ',' ...
        'od ' datestr(tabs(1),'dd-mmm-yyyy HH:MM:SS.FFF') ', do '  datestr(tabs(end),'dd-mmm-yyyy HH:MM:SS.FFF') ]);
    for a = zac: size(H.annotation.starttime) %tady jsou casy vsech anotaci

       disp([ num2str(a) ' ' strtrim(H.annotation.event{a}) ' ... ' num2str( H.annotation.starttime(a)) ...
           ' s ... ' num2str( H.annotation.starttime(a)*H.samplerate(1))]);
       cas = H.annotation.starttime(a); % v sekundach od zacatku zaznamu
       %cas0 = iff(a == zac,H.annotation.starttime(1),H.annotation.starttime(a-1));
       if a==zac, figure('Name','Anotace'); end
       plot([cas cas],[0 1]); %vykreslim cas anotaci do grafu
       hold on;
    end
elseif exist('evts','var') && ~isempty(evts)
    disp(['od ' datestr(tabs(1),'dd-mmm-yyyy HH:MM:SS.FFF') ', do '  datestr(tabs(end),'dd-mmm-yyyy HH:MM:SS.FFF') ]);
    for a = zac: size(evts,2)
        eventtime = datenum(evts(a).dateStr);
        if eventtime >= tabs(1) && eventtime <= tabs(end) % protoze petr jezdik dava do vsech souboru vsechny eventy
            if ~isempty(evts(a).annotation)
                disp([num2str(a) ', ' evts(a).dateStr ', ' evts(a).annotation   ]);
            end
        end
    end
end

end

