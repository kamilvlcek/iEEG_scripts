function [ casstr ] = cas(t, index )
%CAS prevede udaj v poli t na retezec HH:MM:SS.FFF
%   t je soucast dat EEG z Motola

casstr = datestr(t(index)/3600/24,'HH:MM:SS.FFF');

end

