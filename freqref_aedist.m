function [ frekvence,reference ] = freqref_aedist()
%FREQ_AEDIST vrati frekvence pro analyzu
%   Detailed explanation goes here

frekvence = struct;
f=1;
frekvence(f).todo = 1;
frekvence(f).freq = [];
frekvence(f).freqname = 'ERP'; % ERP
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 50:10:150;
frekvence(f).freqname = '50-150'; % broad band gamma
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 7:2:15;
frekvence(f).freqname = '7-15'; % alpha
f=f+1;
frekvence(f).todo = 0;
frekvence(f).freq = 50:5:120;
frekvence(f).freqname = '50-120'; %gamma 2
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 30:5:50;
frekvence(f).freqname = '30-50'; % gamma
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 15:3:31;
frekvence(f).freqname = '15-31'; % beta
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 4:2:8;
frekvence(f).freqname = '4-8'; % theta fast
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 4:1:8;
frekvence(f).freqname = '4-8M'; % theta fast Morlet
frekvence(f).classname = 'Morlet'; % 
f=f+1;
frekvence(f).todo = 0; %hilbertJirka: spatne definovany filtr 1-3.9 Hz - driv to fungovalo?
frekvence(f).freq = 1:3:4; 
frekvence(f).freqname = '1-4'; % theta slow
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 1:1:4;
frekvence(f).freqname = '1-4M'; % theta slow Morlet
frekvence(f).classname = 'Morlet'; % 
f=f+1;
frekvence(f).todo = 1;
frekvence(f).freq = 2:2:150;
frekvence(f).freqname = '2-150'; % all range
frekvence(f).prekryv = 0.5; % 50% prekryv sousednich frekvencnich pasem



reference = struct;
r=1;
reference(r).todo = 1;
reference(r).name = 'refOrig';
reference(r).char = '';
r=2;
reference(r).todo = 1;
reference(r).name = 'refEle';
reference(r).char = 'e';
r=3;
reference(r).todo = 1;
reference(r).name = 'refHead';
reference(r).char = 'h';
r=4;
reference(r).todo = 1;
reference(r).name = 'refBipo';
reference(r).char = 'b';

end
