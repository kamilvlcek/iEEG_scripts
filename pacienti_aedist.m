function [ pacienti ] = pacienti_aedist()
%PACIENTI_AEDIST Summary of this function goes here
%   Detailed explanation goes here


pacienti = struct;
p = 1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p073 Pech VT6';
pacienti(p).data = 'VT6_INV Test Vlcek1_X_aedist.mat'; %'VT6_INV Test Vlcek1_X_aedist.mat';
pacienti(p).header = 'p73_header.mat'; % p73_header.mat je ten nejnovejsi od Jirky - 26.5.2017 'p73_header_kamil.mat';
pacienti(p).psychopy = 'p73_aedist.mat';
pacienti(p).rjepoch = 'p073_aedist_RjEpoch.mat'; %2017 
pacienti(p).epievents = 'p73_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p079 Plu VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_aedist.mat'; %'VT8_2015-04-09_09-46_001_concat_X_aedist.mat';
pacienti(p).header = 'P79_header.mat';
pacienti(p).psychopy = 'p79_aedist.mat';
pacienti(p).rjepoch = 'p079_aedist_RjEpoch.mat'; %'aedist RjEpoch Resp.mat';
pacienti(p).epievents = 'p79_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47 64 114]; %#ok<NBR%#ok<MSNU> AK> 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p082 Vov VT9';
pacienti(p).data = 'VT9_2015-04-21_09-46_001_concat_X_aedist.mat';
pacienti(p).header = 'P82_header.mat';
pacienti(p).psychopy = 'p82_aedist.mat';
pacienti(p).rjepoch = 'p082_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p082_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p083 Kol VT10';
pacienti(p).data = 'VT10_2015-05-19_10-00_001_X_aedist.mat';
pacienti(p).header = 'P83_header.mat';
pacienti(p).psychopy = 'p83_aedist.mat';
pacienti(p).rjepoch = 'p083_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p083_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p095 Hav VT11';
pacienti(p).data = 'VT11_2015-12-15_aedist.mat';
pacienti(p).header = 'P95_header.mat';
pacienti(p).psychopy = 'p95_aedist.mat';
pacienti(p).rjepoch = 'p095_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p95_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1;

pacienti(p).todo = 1;
pacienti(p).folder = 'p096 Gro VT12';
pacienti(p).data = 'VT12_2016-01-26_09-16_001_aedist.mat';
pacienti(p).header = 'P96_headerX.mat';
pacienti(p).psychopy = 'p096_aedist.mat';
pacienti(p).rjepoch = ''; %zatim neexistuje
pacienti(p).epievents = 'p096_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

% ve vsech udalostech jsou epi eventy - neda se vubec pouzit
% pacienti(p).todo = 0;
% pacienti(p).folder = 'p097 Nov VT13';
% pacienti(p).data = 'VT13_2016-02-11_09-20_001 aedist.mat';
% pacienti(p).header = 'P97_header.mat';
% pacienti(p).psychopy = 'p97_aedist.mat';
% pacienti(p).rjepoch = '';
% pacienti(p).epievents = 'p097_aedist_epievents.mat'; %2017
% pacienti(p).rjch = [47]; %#ok<NBRAK>
p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p110 Sou VT14';
pacienti(p).data = 'P110_2016-06-08_15-56_001_concat_aedist.mat';
pacienti(p).header = 'p110_header.mat';
pacienti(p).psychopy = 'p110_aedist.mat';
pacienti(p).rjepoch = 'p110_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p110_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p126 Sve VT15';
pacienti(p).data = 'VT15_2016-09-06_09-18_001_concat_aedist.mat';
pacienti(p).header = 'P126_header.mat';
pacienti(p).psychopy = 'p126_aedist.mat';
pacienti(p).rjepoch = 'p126_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p126_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47 50]; 

p=p+1;
% neni hammer header, ani aedist udelane - nedelali jsme 
% pacienti(p).todo = 1;
% pacienti(p).folder = 'p119 Buc VT16';
% pacienti(p).data = 'VT16_2016-10-10_17-11_001_concat_ppa.mat';
% pacienti(p).header = 'p119_header_kamil.mat';
% pacienti(p).psychopy = '';
% pacienti(p).rjepoch = '';
% pacienti(p).epievents = ''; %2017 
% pacienti(p).rjch = [47 50]; 


pacienti(p).todo = 1;
pacienti(p).folder = 'p130 Koc VT17';
pacienti(p).data = 'VT17_2016-10-24_15-41_001_concat_aedist.mat';
pacienti(p).header = 'p130_header.mat';
pacienti(p).psychopy = 'p130_aedist.mat';
pacienti(p).rjepoch = 'p130_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p130_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; 

p=p+1;
% neni hammer header - ma hipp body
% pacienti(p).todo = 1;
% pacienti(p).folder = 'p132 Pol VT18';
% pacienti(p).data = 'VT18_2016-12-06_10-53_001_aedist.mat';
% pacienti(p).header = 'p130_header.mat';
% pacienti(p).psychopy = 'p130_aedist.mat';
% pacienti(p).rjepoch = 'p130_aedist_RjEpoch.mat';
% pacienti(p).epievents = 'p130_aedist_epievents.mat'; %2017 
% pacienti(p).rjch = [47]; 

pacienti(p).todo = 1;
pacienti(p).folder = 'p129 Kuch VT19';
pacienti(p).data = 'VT19_2017-01-16_09-38_001_concat_aedist.mat';
pacienti(p).header = 'p129_header.mat';
pacienti(p).psychopy = 'p129 aedist.mat';
pacienti(p).rjepoch = 'p129_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p129_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [ 34 35 47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p136 Men VT20';
pacienti(p).data = 'VT20_2017-01-30_09-40_001_concat_aedist.mat';
pacienti(p).header = 'p136_header.mat';
pacienti(p).psychopy = 'p136_aedist.mat';
pacienti(p).rjepoch = 'p136_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p136_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; 

% je to grid, nejde tam asi udelat bipolarni reference, 
% p=11;
% pacienti(p).todo = 1; 
% pacienti(p).folder = 'p138 Ven VT21';
% pacienti(p).data = 'VT21_2017-02-28_09-35_001_500hz_concat_aedist.mat';
% pacienti(p).header = 'p138_header.mat';
% pacienti(p).psychopy = 'aedist_p138.mat';
% pacienti(p).rjepoch = 'p138_aedist_RjEpoch.mat';
% pacienti(p).epievents = 'p138_aedist_epievents.mat'; %2017 
% pacienti(p).rjch = []; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p142 Nam VT22';
pacienti(p).data = 'VT22_2017-03-31_19-27_002_concat_aedist.mat';
pacienti(p).header = 'p142_header.mat';
pacienti(p).psychopy = 'aedist_p142.mat';
pacienti(p).rjepoch = 'p142_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p142_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [2 47 81]; 
end

