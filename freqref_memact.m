function [ frekvence,reference ] = freqref_memact()
% freqref_memact returns the frequencies for analysis

frekvence = struct;
f=1;
frekvence(f).todo = 0;
frekvence(f).freq = [];
frekvence(f).freqname = 'ERP'; % ERP
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 50:5:150; %20 pasem
frekvence(f).freqname = '50-150Hz'; % broad band gamma
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 8:0.5:13; %16 pasem
frekvence(f).freqname = '8-13Hz'; % alpha
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 50:5:120;
frekvence(f).freqname = '50-120Hz'; %gamma 2
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 30:1:50; %20 pasem
frekvence(f).freqname = '30-50Hz'; % gamma
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 15:1:40; %30 pasem % 26
frekvence(f).freqname = '15-40Hz'; % beta
f=f+1;
frekvence(f).todo = 1;
% frekvence(f).freq = 4:0.5:8; %20 pasem % 8 bands
% frekvence(f).freqname = '4-8HzM'; % theta fast Morlet
frekvence(f).freq = 4:1:13; %for connectivity analysis 4-13hz
frekvence(f).freqname = '4-13HzM'; % theta fast Morlet
frekvence(f).classname = 'Morlet'; % 
f=f+1;
frekvence(f).todo = 0; %hilbertJirka: spatne definovany filtr 1-3.9 Hz - driv to fungovalo?
frekvence(f).freq = 1:3:4; 
frekvence(f).freqname = '1-4Hz'; % theta slow
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 1:0.2:4; %15 pasem
frekvence(f).freqname = '1-4HzM'; % theta slow Morlet
frekvence(f).classname = 'Morlet'; % 
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 2:2:150; 
frekvence(f).freqname = '2-150Hz'; % all range
frekvence(f).prekryv = 0.5; % 50% prekryv sousednich frekvencnich pasem
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 1:0.25:10;
frekvence(f).freqname = '1-10HzM'; % theta slow Morlet
frekvence(f).classname = 'Morlet'; % 



reference = struct;
r=1;
reference(r).todo = 0;
reference(r).name = 'refOrig';
reference(r).char = '';
r=2;
reference(r).todo = 0;
reference(r).name = 'refEle';
reference(r).char = 'e';
r=3;
reference(r).todo = 0;
reference(r).name = 'refHead';
reference(r).char = 'h';
r=4;
reference(r).todo = 1;
reference(r).name = 'refBipo';
reference(r).char = 'b';

end

