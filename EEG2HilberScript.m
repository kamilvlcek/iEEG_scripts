%24.3.2015

channels = 1:101; %cisla kanalu k prevedeni na obalku
frekvence = 50:10:150; %gamma 50-150 Hz

% prevedu LPF na Hilbertovu obalku u aktualniho datasetu
EEG = EEG2Hilbert(EEG,channels,frekvence);

% necham EEG lab ulozit prevedene EEG do noveho datasetu,
% viz http://sccn.ucsd.edu/eeglab/structtut/datastructold.html
% http://sccn.ucsd.edu/wiki/Chapter_02:_Writing_EEGLAB_Scripts

[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

eeglab redraw; %prekreslim okno eeglab, aby se tam objevil novy dataset