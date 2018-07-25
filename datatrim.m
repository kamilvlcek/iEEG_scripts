function [d,tabs,header] = datatrim(d,tabs,ts0,ts1,header)
% redukce dat na vybrany test
% podle ts0 a ts1 - indexy zacatku a konce
% updatuju pole d, tabs a header
% since 1.3.201

%global d tabs; %header nebudu davat globalne, muze se jmenovat pokazde jinak

d = d(ts0:ts1,:);
tabs = tabs(ts0:ts1,:);
if exist('header','var')
    header.starttime = datestr(tabs(1),'HH:MM:SS.FFF');
    header.records = (tabs(end)-tabs(1))*60*24*60; %timestampy jsou ve dnech, jejich rozdil je pocet dnu, ja chci vteriny
    %nemenim anotace, takze ty pak nesedi
else
    header = 0;
end
end