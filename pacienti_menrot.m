function [ pacienti ] = pacienti_menrot()
%PACIENTI_menrot Summary of this function goes here
%   Detailed explanation goes here


pacienti = struct;
p = 1;
pacienti(p).todo = 0; %rozdeleno na dve casti, je nutne spojit, p073_menrot neni
pacienti(p).folder = 'p073 Pech VT6';
pacienti(p).data = 'VT6_INV Test Vlcek1_X_menrot.mat'; %'VT6_INV Test Vlcek1_X_menrot.mat';
pacienti(p).header = 'p73_headerX.mat'; % p73_header.mat je ten nejnovejsi od Jirky - 26.5.2017 'p73_header_kamil.mat';
pacienti(p).psychopy = 'p073_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat'; %2017 
pacienti(p).epievents = 'p073_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47 68]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p079 Plu VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_menrot.mat'; %'VT8_2015-04-09_09-46_001_concat_X_menrot.mat';
pacienti(p).header = 'P79_headerK.mat';
pacienti(p).psychopy = 'p079_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat'; %'menrot RjEpoch Resp.mat';
pacienti(p).epievents = 'p079_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47 64 114]; 

p=p+1;
pacienti(p).todo = 0; %docasne, ma vic vzorku 512 nez ostatni 256, nejde do CHilbertMulti
pacienti(p).folder = 'p082 Vov VT9';
pacienti(p).data = 'VT9_2015-04-21_09-46_001_concat_X_menrot.mat';
pacienti(p).header = 'P82_header.mat';
pacienti(p).psychopy = 'p082_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p082_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47 68]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p083 Kol VT10';
pacienti(p).data = 'VT10_2015-05-19_10-00_001_X_menrot.mat';
pacienti(p).header = 'P83_headerK.mat';
pacienti(p).psychopy = 'p083_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p083_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>


p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p096 Gro VT12';
pacienti(p).data = 'VT12_2016-01-26_10-16_002_menrot.mat';
pacienti(p).header = 'P96_headerX.mat';
pacienti(p).psychopy = 'p096_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p096_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

% p=p+1;
% % ve vsech udalostech jsou epi eventy - neda se vubec pouzit
% % pacienti(p).todo = 0;
% % pacienti(p).folder = 'p097 Nov VT13';
% % pacienti(p).data = 'VT13_2016-02-11_09-20_001 menrot.mat';
% % pacienti(p).header = 'P97_header.mat';
% % pacienti(p).psychopy = 'p97_menrot.mat';
% % pacienti(p).rjepoch = '';
% % pacienti(p).epievents = 'p097_menrot_epievents.mat'; %2017
% % pacienti(p).rjch = [47]; %#ok<NBRAK>
% 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p110 Sou VT14';
pacienti(p).data = 'P110_2016-06-08_15-56_001_concat_menrot.mat';
pacienti(p).header = 'p110_header.mat';
pacienti(p).psychopy = 'p110_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p110_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p126 Sve VT15';
pacienti(p).data = 'VT15_2016-09-06_09-18_001_concat_menrot.mat';
pacienti(p).header = 'P126_headerX.mat';
pacienti(p).psychopy = 'p126_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p126_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47 50]; 

% p=p+1;
% neni hammer header, ani menrot udelane - nedelali jsme 
% pacienti(p).todo = 1;
% pacienti(p).folder = 'p119 Buc VT16';
% pacienti(p).data = 'VT16_2016-10-10_17-11_001_concat_ppa.mat';
% pacienti(p).header = 'p119_header_kamil.mat';
% pacienti(p).psychopy = '';
% pacienti(p).rjepoch = '';
% pacienti(p).epievents = ''; %2017 
% pacienti(p).rjch = [47 50]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p130 Koc VT17';
pacienti(p).data = 'VT17_2016-10-24_15-41_001_concat_menrot.mat';
pacienti(p).header = 'p130_headerK.mat';
pacienti(p).psychopy = 'p130_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p130_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47,49]; %49 = vice nez 50 procent vyrazenych epoch

 
% p=p+1;
% neni hammer header - ma hipp body
% pacienti(p).todo = 1;
% pacienti(p).folder = 'p132 Pol VT18';
% pacienti(p).data = 'VT18_2016-12-06_10-53_001_menrot.mat';
% pacienti(p).header = 'p132_header.mat';
% pacienti(p).psychopy = '';
% pacienti(p).rjepoch = '';
% pacienti(p).epievents = ''; %2017 
% pacienti(p).rjch = [47]; 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p129 Kuch VT19';
pacienti(p).data = 'VT19_2017-01-16_09-38_001_concat_menrot.mat';
pacienti(p).header = 'p129_headerXK.mat';
pacienti(p).psychopy = 'p129_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p129_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1; %kamil 6.9.2017
pacienti(p).todo = 1;
pacienti(p).folder = 'p136 Men VT20';
pacienti(p).data = 'VT20_2017-01-31_09-12_001_menrot.mat';
pacienti(p).header = 'p136_header.mat';
pacienti(p).psychopy = 'p136_menrot.mat';
pacienti(p).rjepoch = 'p136_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p136_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47,74]; %74 (50 procent vyrazanych epoch)

% % je to grid, nejde tam asi udelat bipolarni reference, 
% % p=p+1;
% % pacienti(p).todo = 1; 
% % pacienti(p).folder = 'p138 Ven VT21';
% % pacienti(p).data = 'VT21_2017-02-28_09-35_001_500hz_concat_menrot.mat';
% % pacienti(p).header = 'p138_header.mat';
% % pacienti(p).psychopy = 'menrot_p138.mat';
% % pacienti(p).rjepoch = 'p138_menrot_RjEpoch.mat';
% % pacienti(p).epievents = 'p138_menrot_epievents.mat'; %2017 
% % pacienti(p).rjch = []; 
% 
p=p+1; %kamil 6.9.2017
pacienti(p).todo = 1; %vyrazeno kvuli WIFI hlavici a frekvenci 500Hz 
pacienti(p).folder = 'p142 Nam VT22';
pacienti(p).data = 'VT22_2017-04-5_10-31_002_concat_500Hz_menrot.mat';
pacienti(p).header = 'p142_64ch_header.mat'; %mereni s wifi
pacienti(p).psychopy = 'p142_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p142_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [2];%#ok<NBRAK>

p=p+1; %kamil 11.1.2018
pacienti(p).todo = 1; 
pacienti(p).folder = 'p153 Ven VT26';
pacienti(p).data = 'VT26_2017-11-08_09-17_001_concat_menrot.mat';
pacienti(p).header = 'p153_header.mat';
pacienti(p).psychopy = 'p153_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p153_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47];%#ok<NBRAK>

p=p+1; %kamil 11.1.2018
pacienti(p).todo = 1; 
pacienti(p).folder = 'p155 Svi VT27';
pacienti(p).data = 'VT27_2017-11-22_08-36_001_concat_menrot.mat';
pacienti(p).header = 'p155_headerK.mat';
pacienti(p).psychopy = 'p155_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p155_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47];%#ok<NBRAK>

p=p+1; %kamil 11.1.2018
pacienti(p).todo = 1; 
pacienti(p).folder = 'p160 Kor VT28';
pacienti(p).data = 'VT28_2017-12-11_09-39_001_menrot.mat';
pacienti(p).header = 'p160_header.mat';
pacienti(p).psychopy = 'p160_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p160_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1; %kamil 25.5.2018
pacienti(p).todo = 1; 
pacienti(p).folder = 'p162 Fra VT29';
pacienti(p).data = 'VT29_2018-01-09_09-39_001_concat_menrot.mat';
pacienti(p).header = 'p162_header2.mat';
pacienti(p).psychopy = 'p162_menrot.mat';
pacienti(p).rjepoch = 'p162_menrot_rjepoch.mat';
pacienti(p).epievents = 'p162_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47  60  61  62  63  64  65  66  71  72  73  74  75  76  77  78  80  81  82  83  84  85  86  87  88  89  91  92  93  94  95  96  97  99 100 101 102 104 105 106 107 108 114 115 116 117 118 119 120]; 

p=p+1; %kamil 15.5.2018
pacienti(p).todo = 1; 
pacienti(p).folder = 'p151 Ber VT31';
pacienti(p).data = 'VT31_2018-02-26_10-54_001_concat_menrot.mat';
pacienti(p).header = 'p151_header.mat';
pacienti(p).psychopy = 'p151_menrot.mat';
pacienti(p).rjepoch = 'p151_menrot_rjepochs.mat';
pacienti(p).epievents = 'p151_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47]; %#ok<NBRAK>


p=p+1; %kamil 15.5.1.2018
pacienti(p).todo = 1; 
pacienti(p).folder = 'p170 Zem VT32';
pacienti(p).data = 'VT32_2018-03-13_10-13_001_concat_menrot.mat';
pacienti(p).header = 'p170_headerX.mat';
pacienti(p).psychopy = 'p170_menrot.mat';
pacienti(p).rjepoch = 'p170_menrot_rjepochs.mat';
pacienti(p).epievents = 'p170_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 34 35 36 37 38 39 40 45 46 48 49 50 51 52 53 54 55 56 57 58 59 60 81 82 83 91];

p=p+1; %kamil 25.5.2018
pacienti(p).todo = 1; %vyrazn protoze ma jinou frekvenci 500Hz
pacienti(p).folder = {'VT33 Hos p099','p099 Hos VT33'};
pacienti(p).data = 'VT33_2018-03-27_10-21_001_500hz_concat_f1Hz_menrot.mat';
pacienti(p).header = 'p099_64ch_header.mat';
pacienti(p).psychopy = 'p099_menrot.mat';
pacienti(p).rjepoch = 'p099_menrot_rjepochs.mat';
pacienti(p).epievents = 'p099_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [7  8  9 10 11 12 13 14 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31 32 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 57 58 59 60];

p=p+1; %kamil 15.5.2018
pacienti(p).todo = 1; 
pacienti(p).folder = 'p173 Mar VT34';
pacienti(p).data = 'VT34_2018-04-10_10-19_001_concat_menrot.mat';
pacienti(p).header = 'p173_headerX.mat';
pacienti(p).psychopy = 'p173_menrot.mat';
pacienti(p).rjepoch = 'p173_menrot_rjepoch.mat';
pacienti(p).epievents = 'p173_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47]; %#ok<NBRAK>

p=p+1; %kamil 22.3.2019, 500Hz
pacienti(p).todo = 1; 
pacienti(p).folder = 'p180 Van VT35';
pacienti(p).data = 'VT35_2018-05-22_19-35_001_500hz_concat_menrot.mat';
pacienti(p).header = 'p180_64ch_header.mat';
pacienti(p).psychopy = 'p180_menrot.mat';
pacienti(p).rjepoch = 'p180_menrot_rjepoch.mat';
pacienti(p).epievents = 'p180_menrot_epievents.mat'; 
pacienti(p).rjch = [33]; %#ok<NBRAK>

p=p+1; %kamil 22.3.2019, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = 'p181 Tis VT36';
pacienti(p).data = 'VT36_2018-06-02_10-34_001_concat_menrot.mat';
pacienti(p).header = 'p181_header.mat';
pacienti(p).psychopy = 'p181_menrot.mat';
pacienti(p).rjepoch = 'p181_menrot_rjepoch.mat';
pacienti(p).epievents = 'p181_menrot_epievents.mat';  
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 4.12.2019, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = {'VT40 Kuch p129','p129 Kuch VT40'};;
pacienti(p).data = 'VT40_2018-11-13_10-10_001_concat_menrot.mat';
pacienti(p).header = 'p129_headerVT40.mat';
pacienti(p).psychopy = 'p129VT40_menrot.mat';
pacienti(p).rjepoch = 'p129VT40_menrot_rjepoch.mat';
pacienti(p).epievents = 'p129VT40_menrot_epievents.mat';  
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 6.12.2019, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = 'p200 Hav VT45';
pacienti(p).data = 'VT45_2019-05-02_15-35_006_34_concat_512hz_menrot.mat';
pacienti(p).header = 'p200_header.mat';
pacienti(p).psychopy = 'p200VT45_menrot.mat';
pacienti(p).rjepoch = 'p200VT45_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p200VT45_menrot_epievents.mat';  
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 9.12.2019, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = 'p190 Hau VT43';
pacienti(p).data = 'VT43_2019-03-11_15-30_001_24_concat_512hz_menrot.mat';
pacienti(p).header = 'p190_header2.mat';
pacienti(p).psychopy = 'p190VT43_menrot.mat';
pacienti(p).rjepoch = 'p190VT43_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p190VT43_menrot_epievents.mat';
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 9.12.2019, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = {'VT47 Tur p183','p183 Tur VT47'};;
pacienti(p).data = 'VT47_2019-05-28_17-35_001_512hz_39_concat_menrot.mat';
pacienti(p).header = 'p183_headerX.mat';
pacienti(p).psychopy = 'p183VT47_menrot.mat';
pacienti(p).rjepoch = 'p183VT47_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p183VT47_menrot_epievents.mat'; 
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 26.11.2020, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = 'p222 Zich VT49';
pacienti(p).data = 'VT49_2019-09-23_09-50_219_32_concat_512hz_menrot.mat';
pacienti(p).header = 'p222_headerEX.mat';
pacienti(p).psychopy = 'p222VT49_menrot.mat';
pacienti(p).rjepoch = 'p222VT49_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p222VT49_menrot_epievents.mat'; 
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 30.11.2020, 512Hz
pacienti(p).todo = 0; % correct header? The number of channels is 180 - our data, 184 data 
pacienti(p).folder = {'VT52 Loj p1855440','p1855440 Loj VT52'};
pacienti(p).data = 'VT52_2019-12-09_08-59_033_18_concat_512hz_menrot.mat';
pacienti(p).header = 'p1855440_header.mat';  % I deleted one empty channel from the header 'p1855440_header.mat' to be the same number as in the data 
%p1855440_headerm.mat - it has 214 channels instead od 184
pacienti(p).psychopy = 'p1855440VT52_menrot.mat';
pacienti(p).rjepoch = 'p1855440VT52_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p1855440VT52_menrot_epievents.mat'; 
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 30.11.2020, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = {'VT53 Kra p1236540','p1236540 Kra VT53'};
pacienti(p).data = 'VT53_2020-01-15_10-03_214_15_concat_512hz_menrot.mat';
pacienti(p).header = 'p1236540_headerX.mat';  
pacienti(p).psychopy = 'p1236540VT53_menrot.mat';
pacienti(p).rjepoch = 'p1236540VT53_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p1236540VT53_menrot_epievents.mat'; 
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 30.11.2020, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = {'VT56 Kus p1967470','p1967470 Kus VT56'};
pacienti(p).data = 'VT56_2020-05-28_11-00_001_34_concat_512hz_menrot.mat';
pacienti(p).header = 'VT56_header.mat';  
pacienti(p).psychopy = 'p1967470VT56_menrot.mat';
pacienti(p).rjepoch = 'p1967470VT56_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p1967470VT56_menrot_epievents.mat'; 
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 1.12.2020, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = {'VT58 Dan p90','p090 Dan VT58'};
pacienti(p).data = 'VT58_2020-09-24_08-06_630_10_concat_512hz_menrot.mat';
pacienti(p).header = 'p090_header.mat';  
pacienti(p).psychopy = 'p090VT58_menrot.mat';
pacienti(p).rjepoch = 'p090VT58_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p090VT58_menrot_epievents.mat'; 
pacienti(p).rjch = []; %#ok<NBRAK>

p=p+1; %anna 1.12.2020, 512Hz
pacienti(p).todo = 1; 
pacienti(p).folder = 'p187 Boh VT39';
pacienti(p).data = 'VT39_2018-10-18_12-06_009_6_concat_500hz_menrot.mat';
pacienti(p).header = 'p187_64ch_header.mat'; % why 152 channels in p187_header64.mat?  
pacienti(p).psychopy = 'p187VT39_menrot.mat';
pacienti(p).rjepoch = 'p187VT39_menrot_RjEpoch.mat';
pacienti(p).epievents = 'p187VT39_menrot_epievents.mat'; 
pacienti(p).rjch = []; %#ok<NBRAK>

pacienti = pacientFolderSelect(pacienti,setup_menrot());

end

