function [ pacienti ] = pacienti_memact()

pacienti = struct;
p = 1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p1883612 And VT59','VT59 And p1883612'};
pacienti(p).data = 'VT59_2020-12-11_10-52_017_512hz_25_concat_memact.mat'; 
pacienti(p).header = 'VT59_p1883612_headerS.mat'; 
pacienti(p).psychopy = 'VT59_p1883612_memact.mat';
pacienti(p).rjepoch = 'VT59_p1883612_memact_rjepoch.mat';
pacienti(p).epievents = 'VT59_p1883612_memact_epievents.mat';
pacienti(p).rjch = [113 128 185 210]; 

p=p+1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p1239007 Sko VT61','VT61 Sko p1239007'};
pacienti(p).data = 'VT61_2021-10-06_15-27_093_512hz_21_concat_memact.mat'; 
pacienti(p).header = 'p1239007_headerS.mat'; 
pacienti(p).psychopy = 'p1239007_memact.mat';
pacienti(p).rjepoch = 'p1239007_memact_rjepoch.mat'; 
pacienti(p).epievents = 'p1239007_memact_epievents.mat';
pacienti(p).rjch = [68 69];

p=p+1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p2179801 Syk VT62','VT62 Syk p2179801'};
pacienti(p).data = 'VT62_2022-05-17_11-27_174_512hz_12_concat_memact.mat'; 
pacienti(p).header = 'p2179801_headerS.mat'; 
pacienti(p).psychopy = 'p2179801_memact.mat';
pacienti(p).rjepoch = 'p2179801_memact_rjepoch.mat'; 
pacienti(p).epievents = 'p2179801_memact_epievents.mat';
pacienti(p).rjch = [13];

p=p+1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p65639 Kam VT63','VT63 Kam p65639'};
pacienti(p).data = 'VT63_2023-01-25_16-00_100_512hz_24_concat_memact.mat'; 
pacienti(p).header = 'p65639_headerS.mat';  
pacienti(p).psychopy = 'p65639_memact.mat';
pacienti(p).rjepoch = 'p65639_memact_rjepoch.mat'; 
pacienti(p).epievents = 'p65639_memact_epievents.mat';
pacienti(p).rjch = [80, 126:128]; % line noise (80) and empty channels

p=p+1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p1554401 Poz VT64','VT64 Poz p1554401'};
pacienti(p).data = 'VT64_2023-02-08_15-24_294_512hz_14_concat_memact.mat'; 
pacienti(p).header = 'p1554401_headerS.mat';  
pacienti(p).psychopy = 'p1554401_memact.mat';
pacienti(p).rjepoch = 'p1554401_memact_rjepoch.mat'; 
pacienti(p).epievents = 'p1554401_memact_epievents.mat';
pacienti(p).rjch = [];

p=p+1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p2185798 Jel VT65','VT65 Jel p2185798'};
pacienti(p).data = 'VT65_2023-03-07_16-42_005_512hz_22_concat_memact.mat'; 
pacienti(p).header = 'p2185798_headerS.mat'; 
pacienti(p).psychopy = 'p2185798_memact.mat';
pacienti(p).rjepoch = 'p2185798_memact_rjepoch.mat'; 
pacienti(p).epievents = 'p2185798_memact_epievents.mat';
pacienti(p).rjch = [129]; % noise

p=p+1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p1883612 And VT66','VT66 And p1883612'};
pacienti(p).data = 'VT66_2023-04-24_16-52_067_512hz_14_concat_memact.mat'; 
pacienti(p).header = 'VT66_p1883612_headerS.mat'; 
pacienti(p).psychopy = 'VT66_p1883612_memact.mat';
pacienti(p).rjepoch = 'VT66_p1883612_memact_rjepoch.mat';
pacienti(p).epievents = 'VT66_p1883612_memact_epievents.mat';
pacienti(p).rjch = []; 

p=p+1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p1396542 Seb VT67','VT67 Seb p1396542'};
pacienti(p).data = 'VT67_2023-09-25_16-25_104_512hz_26_concat_memact.mat'; 
pacienti(p).header = 'p1396542_headerS.mat'; 
pacienti(p).psychopy = 'p1396542_memact.mat';
pacienti(p).rjepoch = 'p1396542_memact_rjepoch.mat';
pacienti(p).epievents = 'p1396542_memact_epievents.mat';
pacienti(p).rjch = []; 

p=p+1;
pacienti(p).todo = 1; 
pacienti(p).folder = {'p2092846 Gru VT68','VT68 Gru p2092846'};
pacienti(p).data = 'VT68_2023-11-07_14-20_084_512hz_20_concat_memact.mat'; 
pacienti(p).header = 'p2092846_headerS.mat'; 
pacienti(p).psychopy = 'p2092846_memact.mat';
pacienti(p).rjepoch = 'p2092846_memact_rjepoch.mat';
pacienti(p).epievents = 'p2092846_memact_epievents.mat';
pacienti(p).rjch = [73 114 127]; % noise

pacienti = pacientFolderSelect(pacienti,setup_memact());
end

