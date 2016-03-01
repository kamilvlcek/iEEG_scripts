% redukce dat na vybrany test
% podle ts0 a ts1 - indexy zacatku a konce
% updatuju pole d, tabs a header
% since 1.3.2016

d = d(ts0:ts1,:);
tabs = tabs(ts0:ts1,:);
header.starttime = datestr(tabs(1),'HH:MM:SS.FFF');
