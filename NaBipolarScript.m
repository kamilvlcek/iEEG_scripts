

%els = [ 14 25 37 47 57 70 78 91 104 115 125]; %hranice elektrod pro p48 Sruma -C64+64 - konce
%els = [5 10 20 29 34 40 44 49 54 58 62]; %hranice elektrod pro p48 Sruma -Wifi64 - konce
els = [ 11 21 31 43 53 64 75 92 101]; %hranice elektrod pro p73 - konce
%spocitam bipolarni referenci
%0=bipolarni, -1 zadna
EEG.data = bipolarRef(EEG.data,els,0); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw; %prekreslim okno eeglab,0 aby se tam objevil novy dataset