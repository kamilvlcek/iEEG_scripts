function filterMatrix = createSpatialFilter_kisarg(H, N_inputCh, filterSettings)
% creates a spatial filter for (intracranial) EEG data
% input vars:
%   H: header structure
%   N_inputCh: number of input channels
%   filterSettings = struct with fields: 
%         name: 'car'
%     chGroups: 'perElectrode' OR 'perHeadbox'
% output var:
%   filterMatrix: N_inputCh x N_inputCh
% example: 
%   S_car = createSpatialFilter(H, 110, filterSettings);
% usage:
%   X_car = X_ref * S_car;   where:
%       X_ref = [samples x N_inputCh] data matrix
%       X_car = [samples x N_inputCh] data matrix
%       S_car = [N_inputCh x N_inputCh] spatial filter (CAR) matrix

% (c) Jiri, May16


%  init
filterFound = false;
           
%% CAR: common average reference
if strfind(filterSettings.name, 'car')        % ~ common average re-reference (CAR)
    filterFound = true;
    
    % define channel groups
    chGroups = getChannelGroups_kisarg(H, filterSettings.chGroups);
    
    % design filter
    filterMatrix = zeros(N_inputCh);                        % init
    for grp = 1:size(chGroups,2)
        selCh = chGroups{grp};
        numCh = size(selCh,2);
        filterMatrix(selCh,selCh) = eye(numCh) - 1/numCh.*ones(numCh);    % set weights for CAR channels
    end    
end

%% BIP: bipolar reference
if strcmp(filterSettings.name, 'bip')
    filterFound = true;
    
    % define channel groups
    chGroups = getChannelGroups_kisarg(H, 'bip');
    
    % design filter
    filterMatrix = zeros(N_inputCh, size(chGroups,2));      % init
    selCh_H = [];
    for grp = 1:size(chGroups,2)
        selCh = chGroups{grp};
        filterMatrix(selCh(1),grp) = 1;                     % set weights for BIP channels
        if size(selCh,2) == 2
            filterMatrix(selCh(2),grp) = -1;                % set weights for BIP channels
        else
            warning(['BIP: only 1 channel on electrode shank, no referencing. Channel = ' num2str(grp)]);
        end
        selCh_H = cat(2, selCh_H, selCh(1));
    end  
end

%% NAN: no spatial filter
if strcmp(filterSettings.name, 'nan')
    filterFound = true;
    filterMatrix = eye(N_inputCh);
end

assert(filterFound);
