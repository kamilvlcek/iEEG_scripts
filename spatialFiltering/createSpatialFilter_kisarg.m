function filterMatrix = createSpatialFilter_kisarg(H, N_inputCh, filterSettings,rejCh)
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
    
    % rejected channels
    %rejCh = [];         % set here indices of rejected channels !!!
    
    % design filter
    filterMatrix = zeros(N_inputCh);                        % init
    for grp = 1:size(chGroups,2)
        selCh = setdiff(chGroups{grp}, rejCh);                              % potential bug: if chGroups{grp} == rejCh (lets hope not!)       
        if(~isempty(selCh)) %pokud vsecny kanaly za elektrodu vyrazene, muzu pokracovat dalsi elektrodou
            numCh = size(selCh,2);
            filterMatrix(selCh,selCh) = eye(numCh) - 1/numCh.*ones(numCh);    % set weights for CAR channels
        end
    end    
end

%% BIP: bipolar reference
if strcmp(filterSettings.name, 'bip')
    filterFound = true;
    
    % define channel groups
    chGroups = getChannelGroups_kisarg(H, 'bip',rejCh);
    
    % rejected channels
    rejChJirka = [];         % set here indices of rejected channels !!!
    % odstranuju RjCh v modifikovane getChannelGroups_kisarg, takze tohle bude vzdy prazdne - kamil 14.6.2016
    
    % design filter
    filterMatrix = zeros(N_inputCh, size(chGroups,2));      % init
    selCh_H = [];
    for grp = 1:size(chGroups,2)
        selCh = setdiff(chGroups{grp}, rejChJirka );              % potential bug: if chGroups{grp} == rejCh (lets hope not!)
        assert(~isempty(selCh));
        filterMatrix(selCh(1),grp) = 1;                     % set weights for BIP channels
        if size(selCh,2) == 2
            filterMatrix(selCh(2),grp) = -1;                % set weights for BIP channels
        else
            warning(['BIP: only 1 channel on electrode shank, no referencing. Channel = ' num2str(grp)]); %zprava kvuli rjch, ted nikdy nenastane - kamil 14.6.2016
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
