function [ pacienti ] = pacienti_ppa(  )
%PACIENTI_PPA Summary of this function goes here
%   Detailed explanation goes here

pacienti = struct;
p = 1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p073 Pech VT6'; %'p073_VT6';
pacienti(p).data = 'VT6_INV Test Vlcek1_X_ppa.mat';
pacienti(p).header = 'p73_header.mat';
pacienti(p).psychopy = 'p73_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p73_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p079 Plu VT8'; %'p79_VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_ppa.mat';
pacienti(p).header = 'P79_header.mat';
pacienti(p).psychopy = 'p79_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p79.mat'; 
pacienti(p).epievents = 'p79_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47 64 68]; %#ok<NBR%#ok<MSNU> AK> puvodne 47 64 114

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p082 Vov VT9'; %'VT09';
pacienti(p).data = 'VT9_2015-04-21_09-46_001_concat_X_ppa.mat';
pacienti(p).header = 'P82_header.mat';
pacienti(p).psychopy = 'p82_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p82_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47 68 126]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p083 Kol VT10'; %'VT10';
pacienti(p).data = 'VT10_2015-05-19_10-00_001_X_ppa.mat';
pacienti(p).header = 'P83_header.mat';
pacienti(p).psychopy = 'p83_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p83.mat';
pacienti(p).epievents = 'p83_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47 64]; %  64 ?

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p095 Hav VT11'; %'VT11';
pacienti(p).data = 'VT11_2015-12-15_ppa.mat';
pacienti(p).header = 'P95_header.mat';
pacienti(p).psychopy = 'p95_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p95.mat';
pacienti(p).epievents = 'p95_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK> 64?

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p096 Gro VT12'; %'VT12';
pacienti(p).data = 'VT12_2016-01-26_09-16_001_concat_ppa.mat';
pacienti(p).header = 'P96_header.mat';
pacienti(p).psychopy = 'p96_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p96.mat'; 
pacienti(p).epievents = 'p96_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK> 64?

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p097 Nov VT13'; %'VT13';
pacienti(p).data = 'VT13_2016-02-11_09-20_001_concat_ppa.mat';
pacienti(p).header = 'P97_header.mat';
pacienti(p).psychopy = 'p97_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch-p97.mat';
pacienti(p).epievents = 'p97_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p110 Sou VT14'; %'VT14';
pacienti(p).data = 'P110_2016-06-08_15-56_001_concat_ppa.mat';
pacienti(p).header = 'p110_header.mat';
pacienti(p).psychopy = 'p110_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p110_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47 15]; %#ok<NB%#ok<MSNU> RAK>

p=p+1; % moc zachvatu?
pacienti(p).todo = 1;
pacienti(p).folder = 'p126 Sve VT15'; %'VT15';
pacienti(p).data = 'VT15_2016-09-06_09-18_001_concat_ppa.mat';
pacienti(p).header = 'P126_headerX.mat';
pacienti(p).psychopy = 'p126_ppa.mat';
pacienti(p).rjepoch = ''; %muze byt prazne, pak se nevyrazuji zadne epochy
pacienti(p).epievents = 'p126_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47 50]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 0;
pacienti(p).folder = 'p119 Buc VT16'; %'VT16';
pacienti(p).data = 'VT16_2016-10-10_17-11_001_concat_ppa.mat';
pacienti(p).header = 'p119_header_kamil.mat';%chybi header Hammer
pacienti(p).psychopy = 'p119_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p119_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47 57 64]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p130 Koc VT17'; %'VT17';
pacienti(p).data = 'VT17_2016-10-24_15-41_001_concat_ppa.mat';
pacienti(p).header = 'p130_header.mat';
pacienti(p).psychopy = 'p130_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p126_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p132 Pol VT18'; %'VT18';
pacienti(p).data = 'VT18_2016-12-07_16-20_001_ppa_512hz.mat';
pacienti(p).header = 'p132_header.mat';% lisi se
pacienti(p).psychopy = 'p132_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p132_ppa_epievents.mat'; %2017 
pacienti(p).rjch = []; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p129 Kuch VT19'; %'VT19';
pacienti(p).data = 'VT19_2017-01-16_09-38_001_concat_ppa.mat';
pacienti(p).header = 'p129_header.mat';
pacienti(p).psychopy = 'p129_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p129_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47];  %POZOR na elektrody 34 a 35 od epochy 166-169 + elektroda 67 epocha 182 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p136 Men VT20'; %'VT20';
pacienti(p).data = 'VT20_2017-01-30_09-40_001_concat_ppa.mat';
pacienti(p).header = 'p136_header.mat'; %'P83_64ch_header.mat'
pacienti(p).psychopy = 'p136_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p136_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p138 Ven VT21'; %'VT21';
pacienti(p).data = 'VT21_2017-02-28_09-35_001_500hz_concat_ppa.mat'; %wifi data
pacienti(p).header = 'p138_header.mat'; %'P83_64ch_header.mat'
pacienti(p).psychopy = 'p138_ppa.mat';
pacienti(p).rjepoch = 'ppa_RjEpoch.mat';
pacienti(p).epievents = 'p138_ppa_epievents.mat'; %2017 
pacienti(p).rjch = []; %

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p155 Svi VT27'; %'VT21';
pacienti(p).data = 'VT27_2017-11-22_08-36_001_concat_ppa.matt';
pacienti(p).header = 'p155_header.mat'; %'P83_64ch_header.mat'
pacienti(p).psychopy = 'p155_ppa.mat';
pacienti(p).rjepoch = 'p155_ppa_rjepoch.mat';
pacienti(p).epievents = 'p155_ppa_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %

pacienti(p).todo = 1;
pacienti(p).folder = 'p160 Kor VT28';
pacienti(p).data = 'VT28_2017-12-11_16-05_001_concat_ppa.mat';
pacienti(p).header = 'p160_header.mat';
pacienti(p).psychopy = 'p160_ppa.mat';
pacienti(p).rjepoch = 'p160_ppa_rjepoch.mat';
pacienti(p).epievents = 'p160_ppa_epievents.mat'; 
pacienti(p).rjch = [47,5,79,105, 106, 107, 122];

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p162 Fra VT29';
pacienti(p).data = 'VT29_2018-01-09_10-04_001_ppa.mat';
pacienti(p).header = 'p162_header.mat';
pacienti(p).psychopy = 'p162_ppa.mat';
pacienti(p).rjepoch = 'p162_ppa_rjepoch.mat';
pacienti(p).epievents = 'p162_ppa_epievents.mat'; 
pacienti(p).rjch = [47,72,73,74,75,80,81,82,83,84,85,86,87,91,92,93,94,95,96,104,105,106,107,108,114,115,116,117,118]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p151 Ber VT31';
pacienti(p).data = 'VT31_2018-02-26_10-54_001_concat_ppa.mat';
pacienti(p).header = 'p151_header.mat';
pacienti(p).psychopy = 'p151_ppa.mat';
pacienti(p).rjepoch = 'p151_ppa_rjepoch.mat';
pacienti(p).epievents = 'p151_ppa_epievents.mat'; 
pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p170 Zem VT32';
pacienti(p).data = 'VT32_2018-03-13_10-13_001_concat_ppa.mat';
pacienti(p).header = 'p170_headerX.mat';
pacienti(p).psychopy = 'p170_ppa.mat';
pacienti(p).rjepoch = 'p170_ppa_rjepoch.mat';
pacienti(p).epievents = 'p170_ppa_epievents.mat'; 
pacienti(p).rjch = [1,2,9,10,14,15,16,18,23,24,25,26,27,28,31,34,35,36,37,38,39,40,45,46,47,48,49,50,51,53,54,55,56,82];

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p173 Mar VT34';
pacienti(p).data = 'VT34_2018-04-09_10-39_001_concat_ppa.mat';
pacienti(p).header = 'p173_header.mat';
pacienti(p).psychopy = 'p173_ppa.mat';
pacienti(p).rjepoch = 'p173_ppa_rjepoch.mat';
pacienti(p).epievents = 'p173_ppa_epievents.mat'; 
pacienti(p).rjch = [47]; 

end

