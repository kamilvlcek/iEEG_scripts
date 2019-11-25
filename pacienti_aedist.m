function [ pacienti ] = pacienti_aedist()
%PACIENTI_AEDIST Summary of this function goes here
%   Detailed explanation goes here


pacienti = struct;
p = 1;
pacienti(p).todo = 1; %docasne ho nechci do analyzy, protoza ma jiny pocet epoch a nesedi do CHIlbertMulti
pacienti(p).folder = 'p073 Pech VT6';
pacienti(p).data = 'VT6_INV Test Vlcek1_X_aedist.mat'; %'VT6_INV Test Vlcek1_X_aedist.mat';
pacienti(p).header = 'p73_headerX.mat'; % p73_header.mat je ten nejnovejsi od Jirky - 26.5.2017 'p73_header_kamil.mat';
pacienti(p).psychopy = 'p73_aedist.mat';
pacienti(p).rjepoch = 'p073_aedist_RjEpoch.mat'; %2017 
pacienti(p).epievents = 'p73_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p079 Plu VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_aedist.mat'; %'VT8_2015-04-09_09-46_001_concat_X_aedist.mat';
pacienti(p).header = 'P79_headerK.mat';
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
pacienti(p).header = 'P83_headerK.mat';
pacienti(p).psychopy = 'p83_aedist.mat';
pacienti(p).rjepoch = 'p083_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p083_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p095 Hav VT11';
pacienti(p).data = 'VT11_2015-12-15_aedist.mat';
pacienti(p).header = 'P95_headerK.mat';
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
pacienti(p).rjepoch = 'p096_aedist_rjepoch.mat'; %zatim neexistuje
pacienti(p).epievents = 'p096_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

% p=p+1;
% ve vsech udalostech jsou epi eventy - neda se vubec pouzit
% 2017 - ve vetsine kanalu je vic nez 30% epoch s epi eventy - nepouzivat
% pacienti(p).todo = 0;
% pacienti(p).folder = 'p097 Nov VT13';
% pacienti(p).data = 'VT13_2016-02-11_09-20_001 aedist.mat';
% pacienti(p).header = 'P97_header.mat';
% pacienti(p).psychopy = 'p97_aedist.mat';
% pacienti(p).rjepoch = ''; %zatim neexistuje
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
pacienti(p).header = 'P126_headerX.mat';
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
pacienti(p).header = 'p130_headerK.mat';
pacienti(p).psychopy = 'p130_aedist.mat';
pacienti(p).rjepoch = 'p130_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p130_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; 

p=p+1;
% neni hammer header - ma hipp body
pacienti(p).todo = 1;
pacienti(p).folder = 'p132 Pol VT18';
pacienti(p).data = 'VT18_2016-12-06_10-53_001_aedist_512hz.mat';
pacienti(p).header = 'p132_header.mat';
pacienti(p).psychopy = 'p132_aedist.mat';
pacienti(p).rjepoch = 'p132_aedist_rjepochs.mat';
pacienti(p).epievents = 'p132_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p129 Kuch VT19';
pacienti(p).data = 'VT19_2017-01-16_09-38_001_concat_aedist.mat';
pacienti(p).header = 'p129_headerXK.mat';
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


p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p144 Nav VT23';
pacienti(p).data = 'VT23_2017-05-09_17-12_001_aedist.mat';
pacienti(p).header = 'p144_header2.mat';
pacienti(p).psychopy = 'p144_aedist.mat';
pacienti(p).rjepoch = 'p144_aedist_RjEpoch.mat'; 
pacienti(p).epievents = 'p144_aedist_epievents.mat'; %2017  
pacienti(p).rjch = [ 47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p147 Fli VT24';
pacienti(p).data = 'VT24_2017-06-19_15-28_001_aedist.mat';
pacienti(p).header = 'p147_header.mat';
pacienti(p).psychopy = 'p147_aedist.mat';
pacienti(p).rjepoch = 'p147_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p147_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p153 Ven VT26';
pacienti(p).data = 'VT26_2017-11-08_09-17_001_concat_aedist.mat';
pacienti(p).header = 'p153_header.mat';
pacienti(p).psychopy = 'p153_aedist.mat';
pacienti(p).rjepoch = 'p153_aedist_RjEpoch.mat';
pacienti(p).epievents = 'p153_aedist_epievents.mat'; %2017 
pacienti(p).rjch = [47]; 
 
p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p155 Svi VT27';
pacienti(p).data = 'VT27_2017-11-22_08-36_001_concat_aedist.mat';
pacienti(p).header = 'p155_headerK.mat';
pacienti(p).psychopy = 'p155_aedist.mat';
pacienti(p).rjepoch = 'p155_aedist_RjEpoch.mat'; %2018
pacienti(p).epievents = 'p155_aedist_epievents.mat'; %2018
pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p160 Kor VT28';
pacienti(p).data = 'VT28_2017-12-11_16-05_001_concat_aedist.mat';
pacienti(p).header = 'p160_header.mat';
pacienti(p).psychopy = 'p160_aedist.mat';
pacienti(p).rjepoch = 'p160_aedist_rjepoch.mat'; %2018
pacienti(p).epievents = 'p160_aedist_epievents.mat'; %2018
pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p162 Fra VT29';
pacienti(p).data = 'VT29_2018-01-05_10-21_001_aedist.mat';
pacienti(p).header = 'p162_header2.mat';
pacienti(p).psychopy = 'p162_aedist.mat';
pacienti(p).rjepoch = 'p162_aedist_RjEpoch.mat'; %2018
pacienti(p).epievents = 'p162_aedist_epievents.mat'; %2018
pacienti(p).rjch = [47 114	36	37	38	39	51	52	59	60	61	62	63	64	65	66	67	68	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	112	113	114	115	116	117	118	119	120]; 
% spousta epileptickych kanalu, abych nevyradil vsechny epochy, vyradil jsem ty kanaly, ktere pres 30% epoch k vyrazeni

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p165 Jan VT30';
pacienti(p).data = 'VT30_2018-01-19_09-21_001_aedist.mat';
pacienti(p).header = 'p165_header.mat';
pacienti(p).psychopy = 'p165_aedist.mat';
pacienti(p).rjepoch = 'p165_aedist_RjEpoch.mat'; %2018-05
pacienti(p).epievents = 'p165_aedist_epievents.mat'; %2018
pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p170 Zem VT32';
pacienti(p).data = 'VT32_2018-03-10_13-39_001_concat_aedist.mat';
pacienti(p).header = 'p170_headerX.mat';
pacienti(p).psychopy = 'p170_aedist.mat';
pacienti(p).rjepoch = 'p170_aedist_rjEpoch.mat'; %2018-05
pacienti(p).epievents = 'p170_aedist_epiEvents.mat'; %2018
pacienti(p).rjch = [47 84]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p099 Hos VT33';
pacienti(p).data = 'VT33_2018-03-23_10-33_001_aedist_512hz.mat';
pacienti(p).header = 'p099_header.mat';
pacienti(p).psychopy = 'p099_aedist.mat';
pacienti(p).rjepoch = 'p099_aedist_rjepoch.mat'; %2018-05
pacienti(p).epievents = 'p099_aedist_epievents.mat'; %2018
pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p173 Mar VT34';
pacienti(p).data = 'VT34_2018-04-09_10-39_001_concat_aedist.mat';
pacienti(p).header = 'p173_headerX.mat';
pacienti(p).psychopy = 'p173_aedist.mat';
pacienti(p).rjepoch = 'p173_aedist_rjepoch.mat'; %2018-05
pacienti(p).epievents = 'p173_aedist_epievents.mat'; %2018
pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1; %20.11.2019 Sofiia
pacienti(p).folder = 'p180 Van VT35';
pacienti(p).data = 'VT35_2018-05-22_17-19_001_500hz_5_concat_aedist.mat';
pacienti(p).header = 'p180_64ch_header.mat';
pacienti(p).psychopy = 'p180_aedist.mat';
pacienti(p).rjepoch = 'p180_aedist_rjepoch.mat'; 
pacienti(p).epievents = 'p180_aedist_epievents.mat'; 
pacienti(p).rjch = [33]; 

p=p+1;
pacienti(p).todo = 1; %21.11.2019 Sofiia
pacienti(p).folder = 'p181 Tis VT36';
pacienti(p).data = 'VT36_2018-06-04_18-05_001_512hz_2_concat_aedist.mat';
pacienti(p).header = 'p181_header.mat';
pacienti(p).psychopy = 'p181_aedist.mat';
pacienti(p).rjepoch = 'p181_aedist_rjepoch.mat'; 
pacienti(p).epievents = 'p181_aedist_epievents.mat'; 
pacienti(p).rjch = []; 

p=p+1;
pacienti(p).todo = 1; %22.11.2019 Sofiia
pacienti(p).folder = 'p187 Boh VT39';
pacienti(p).data = 'VT39_2018-10-13_11-12_001_concat_aedist.mat';
pacienti(p).header = 'p187_headerB.mat';
pacienti(p).psychopy = 'p187_aedist.mat';
pacienti(p).rjepoch = 'p187_aedist_rjepoch.mat'; 
pacienti(p).epievents = 'p187_aedist_epievents.mat'; 
pacienti(p).rjch = [];

p=p+1;
pacienti(p).todo = 1; %22.11.2019 Sofiia
pacienti(p).folder = 'p129 Kuch VT40';
pacienti(p).data = 'VT40_2018-11-14_08-54_001_aedist.mat';
pacienti(p).header = 'p129_headerVT40b.mat';
pacienti(p).psychopy = 'p129vt40_aedist.mat';
pacienti(p).rjepoch = 'p129vt40_aedist_rjepoch.mat'; 
pacienti(p).epievents = 'p129vt40_aedist_epievents.mat'; 
pacienti(p).rjch = [];
end

