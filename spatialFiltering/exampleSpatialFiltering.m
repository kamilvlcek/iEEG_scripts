%% spatial filtering (re-referencing) of iEEG data
% filterName = 'car' (common average reference)
% - option 1 = to re-reference over all channels = 'perHeadbox'
% - option 2 = to re-reference over el. shanks   = 'perElectrode'
% filterName = 'bip' (bipolar reference)
% filterName = 'nan' (no re-reference)

% (c) Jiri, May16

%% load example data (exported from Petr Jezdik)
load('/home/hammer/export/shareData/subjData_carDriving/16_PR1_day1/rawData/amplifierData/P96_2016-01-25_14-30_001.mat');
assert(exist('d','var') == 1);

%% load header 'H'
fileName = ['setHerePath' filesep '16_PR1_day1_header.mat'];  
% for example:
fileName = ['/home/hammer/export/shareData/subjData_carDriving/16_PR1_day1/alignedData' filesep '16_PR1_day1_header.mat'];  % for example
assert(exist(fileName,'file') == 2);
clear H;
load(fileName, 'H');
assert(size(d,2) == size(H.channels,2));

%% selected channels: signal type = iEEG
selCh_H = [];
selSignals = {'SEEG', 'ECoG-Grid', 'ECoG-Strip'};           % select desired channel group
for ch = 1:size(H.channels,2)
    if isfield(H.channels(ch), 'signalType')
        if ismember(H.channels(ch).signalType, selSignals)
            selCh_H = [selCh_H, H.channels(ch).numberOnAmplifier];
        end
    end
end
rawData = d(:,selCh_H);
H.selCh_H = selCh_H;

%% EXAMPLE: design spatial filter: CAR (per headbox)
filterName = 'car';         % options: 'car','bip','nan'
filterSettings.name = filterName;
filterSettings.chGroups = 'perHeadbox';        % options: 'perHeadbox' (=global CAR), OR 'perElectrode' (=local CAR per el. shank)

filterMatrix = createSpatialFilter_kisarg(H, size(rawData,2), filterSettings);
assert(size(rawData,2) == size(filterMatrix,1));

%% apply spatial filter
filtData = rawData * filterMatrix;  
assert(size(filtData,1) == size(rawData,1));



%% EXAMPLE: design spatial filter: CAR (per electrode shank)
filterName = 'car';         % options: 'car','bip','nan'
filterSettings.name = filterName;
filterSettings.chGroups = 'perElectrode';        % options: 'perHeadbox' (=global CAR), OR 'perElectrode' (=local CAR per el. shank)

filterMatrix = createSpatialFilter_kisarg(H, size(rawData,2), filterSettings);
assert(size(rawData,2) == size(filterMatrix,1));

%% apply spatial filter
filtData = rawData * filterMatrix;  
assert(size(filtData,1) == size(rawData,1));



%% EXAMPLE: design spatial filter: BIP (bipolar on electrode shanks)
filterName = 'bip';         % options: 'car','bip','nan'
filterSettings.name = filterName;

filterMatrix = createSpatialFilter_kisarg(H, size(rawData,2), filterSettings);
assert(size(rawData,2) == size(filterMatrix,1));

%% apply spatial filter
filtData = rawData * filterMatrix;  
assert(size(filtData,1) == size(rawData,1));


