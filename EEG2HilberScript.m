%24.3.2015

channels = 1:125; %cisla kanalu k prevedeni na obalku
frekvence = 50:10:150; %gamma 50-150 Hz

%nejdriv prevest referenci na bipolarni
%takhle to moc nefunguje - AEdistP68 AEC gamma Allo.fdt

% prevedu LPF na Hilbertovu obalku u aktualniho datasetu
if size(EEG.data,3) < 2 % pokud se nejedna o data rozdelena do epoch
   
    % necham EEG lab ulozit prevedene EEG do noveho datasetu,
    % viz http://sccn.ucsd.edu/eeglab/structtut/datastructold.html
    % http://sccn.ucsd.edu/wiki/Chapter_02:_Writing_EEGLAB_Scripts

    EEG = EEG2Hilbert(EEG,channels,frekvence);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw; %prekreslim okno eeglab, aby se tam objevil novy dataset
else
    %nic neukladam, jen kvuli spektrogramum
    EEG2Hilbert(EEG,channels,frekvence);
end
