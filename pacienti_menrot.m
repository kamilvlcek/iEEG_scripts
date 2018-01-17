function [ pacienti ] = pacienti_menrot()
%PACIENTI_menrot Summary of this function goes here
%   Detailed explanation goes here


pacienti = struct;
p = 1;
pacienti(p).todo = 0;
pacienti(p).folder = 'p073 Pech VT6';
pacienti(p).data = 'VT6_INV Test Vlcek1_X_menrot.mat'; %'VT6_INV Test Vlcek1_X_menrot.mat';
pacienti(p).header = 'p73_header.mat'; % p73_header.mat je ten nejnovejsi od Jirky - 26.5.2017 'p73_header_kamil.mat';
pacienti(p).psychopy = 'p073_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat'; %2017 
pacienti(p).epievents = 'p73_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p079 Plu VT8';
pacienti(p).data = 'VT8_2015-04-09_09-46_001_concat_X_menrot.mat'; %'VT8_2015-04-09_09-46_001_concat_X_menrot.mat';
pacienti(p).header = 'P79_header.mat';
pacienti(p).psychopy = 'p079_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat'; %'menrot RjEpoch Resp.mat';
pacienti(p).epievents = 'p79_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47 64 114]; %#ok<NBR%#ok<MSNU> AK> 

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p082 Vov VT9';
pacienti(p).data = 'VT9_2015-04-21_09-46_001_concat_X_menrot.mat';
pacienti(p).header = 'P82_header.mat';
pacienti(p).psychopy = 'p082_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p082_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47 68]; %#ok<NB%#ok<MSNU> RAK>

p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p083 Kol VT10';
pacienti(p).data = 'VT10_2015-05-19_10-00_001_X_menrot.mat';
pacienti(p).header = 'P83_header.mat';
pacienti(p).psychopy = 'p083_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p083_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47]; %#ok<NBRAK>


p=p+1;
pacienti(p).todo = 1;
pacienti(p).folder = 'p096 Gro VT12';
pacienti(p).data = 'VT12_2016-01-26_10-16_002_menrot.mat';
pacienti(p).header = 'P96_header.mat';
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
pacienti(p).header = 'P126_header.mat';
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
pacienti(p).header = 'p130_header.mat';
pacienti(p).psychopy = 'p130_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p130_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47,49]; %49 = vice nez 50 procent vyrazenych epoch

 
p=p+1;
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
pacienti(p).header = 'p129_header.mat';
pacienti(p).psychopy = 'p129_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p129_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [47]; 

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
pacienti(p).todo = 1; 
pacienti(p).folder = 'p142 Nam VT22';
pacienti(p).data = 'VT22_2017-04-5_10-31_002_concat_500Hz_menrot.mat';
pacienti(p).header = 'p142_64ch_header.mat';
pacienti(p).psychopy = 'p142_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p142_menrot_epievents.mat'; %2017 
pacienti(p).rjch = [2];

p=p+1; %kamil 11.1.2018
pacienti(p).todo = 0; 
pacienti(p).folder = 'p153 Ven VT26';
pacienti(p).data = 'VT26_2017-11-08_09-17_001_concat_menrot.mat';
pacienti(p).header = 'p153_header.mat';
pacienti(p).psychopy = 'p153_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p153_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47];

p=p+1; %kamil 11.1.2018
pacienti(p).todo = 0; 
pacienti(p).folder = 'p155 Svi VT27';
pacienti(p).data = 'VT27_2017-11-22_08-36_001_concat_menrot.mat';
pacienti(p).header = 'p155_header.mat';
pacienti(p).psychopy = 'p155_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p155_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47];

p=p+1; %kamil 11.1.2018
pacienti(p).todo = 0; 
pacienti(p).folder = 'p160 Kor VT28';
pacienti(p).data = 'VT28_2017-12-11_09-39_001_menrot.mat';
pacienti(p).header = 'p160_header.mat';
pacienti(p).psychopy = 'p160_menrot.mat';
pacienti(p).rjepoch = 'menrot_RjEpoch.mat';
pacienti(p).epievents = 'p160_menrot_epievents.mat'; %2018 
pacienti(p).rjch = [47];

end

