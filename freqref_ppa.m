function [ frekvence,reference ] = freqref_ppa()
%FREQ_AEDIST vrati frekvence pro analyzu
%   Detailed explanation goes here

frekvence = struct;
f=1;
frekvence(f).todo = 0;
frekvence(f).freq = [];
frekvence(f).freqname = 'ERP'; % ERP
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 50:5:150;
frekvence(f).freqname = '50-150Hz'; % broad band gamma
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 7:0.5:15;
frekvence(f).freqname = '7-15Hz'; % alpha
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 50:5:120;
frekvence(f).freqname = '50-120Hz'; %gamma 2
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 30:1:50;
frekvence(f).freqname = '30-50Hz'; % gamma
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 15:0.5:30;
frekvence(f).freqname = '15-30Hz'; % beta
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 4:2:8;
frekvence(f).freqname = '4-8Hz'; % theta fast
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 4:0.2:8;
frekvence(f).freqname = '4-8HzM'; % theta fast Morlet
frekvence(f).classname = 'Morlet'; % 
f=f+1;
frekvence(f).todo = 0; %hilbertJirka: spatne definovany filtr 1-3.9 Hz - driv to fungovalo?
frekvence(f).freq = 1:3:4; 
frekvence(f).freqname = '1-4Hz'; % theta slow
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 1:0.2:4;
frekvence(f).freqname = '1-4HzM'; % theta slow Morlet
frekvence(f).classname = 'Morlet'; % 
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 1:0.2:10;
frekvence(f).freqname = '1-10HzM'; % theta slow Morlet
frekvence(f).classname = 'Morlet'; %
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 2:2:150;
frekvence(f).freqname = '2-150Hz'; % all range
frekvence(f).prekryv = 0.5; % 50% prekryv sousednich frekvencnich pasem



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

