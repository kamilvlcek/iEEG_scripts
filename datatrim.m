function [header] = datatrim(header,ts0,ts1)
% redukce dat na vybrany test
% podle ts0 a ts1 - indexy zacatku a konce
% updatuju pole d, tabs a header
% since 1.3.201

global d tabs; %header nebudu davat globalne, muze se jmenovat pokazde jinak

d = d(ts0:ts1,:);
tabs = tabs(ts0:ts1,:);
header.starttime = datestr(tabs(1),'HH:MM:SS.FFF');
end