function chGroups = getChannelGroups_kisarg(H, groupSpecification,RjCh)
% returns channel groups (cell array) based on group specification

% (c) Jiri, Apr16

chGroups = [];
if ~exist('RjCh','var'), RjCh = []; end
%% selected channels (for example signal type = iEEG)
assert(isfield(H,'selCh_H'));
selCh_H = H.selCh_H;
selCh_raw = 1:size(selCh_H,2);      % corresponds to indices in rawData

%% CAR: channel groups = headboxes of amplifier with different REFs & GNDs
if strcmp(groupSpecification, 'perHeadbox')
    if isfield(H.channels(1), 'headboxNumber')      % re-referencing per headbox
        hdbxAll = [];
        for ch = 1:size(H.channels,2)
            hdbxAll = cat(2, hdbxAll, H.channels(ch).headboxNumber);        % headbox numbers from all channels
        end
        hdbxSel = hdbxAll(selCh_H);                                         % headbox numbers from selected channels (e.g. iEEG channels)
        assert(size(selCh_raw,2) == size(hdbxSel,2));

        hdbxUnq = unique(hdbxSel);                                          % headbox numers
        chGroups = cell(1,size(hdbxUnq,2));
        for grp = 1:size(hdbxUnq,2)
            thisGrp = hdbxUnq(grp);
            for ch = 1:size(selCh_raw,2)
                if hdbxSel(ch) == thisGrp
                    chGroups{grp} = cat(2, chGroups{grp}, ch);
                end
            end
        end
    else                                            % only 1 headbox in recording
        chGroups{1,1} = selCh_raw;
    end
    %v chGroups maji byt dve cell - kazda pro jeden head - 1:64 a 1:61 napriklad
    %v H.channels je headboxNumber 1 nebo 2 
    assert(numel(chGroups)<=2, 'chybny header - vice nez dva headboxy');
end

%% CAR: channel groups = SEEG electrode shanks
if strcmp(groupSpecification, 'perElectrode')
    elsh_all = [];
    for ch = 1:size(H.channels,2)
        elsh_all{ch} = extractFromString(H.channels(ch).name, 'string');    % string part of all channel names
    end    
    elsh_sel = elsh_all(selCh_H);                                           % string part of selected channels
    elsh_unq = unique(elsh_sel);                                            % electrode shank names
    
    chGroups = cell(1,size(elsh_unq,2));
    for grp = 1:size(elsh_unq,2)
        thisGrp = elsh_unq{grp};
        for ch = 1:size(selCh_raw,2)
            if strcmp(extractFromString(H.channels(ch).name, 'string'), thisGrp)
                chGroups{grp} = cat(2, chGroups{grp}, ch);
            end
        end
    end    
end


%% BIP: channel groups = neighboring SEEG channels on same electrode shank
if strcmp(groupSpecification, 'bip')   
    chGroups = [];
    grp = 1;
    chprev = 1; %matlab neumi menit index for smycky za jejiho behu
    for ch = 2:size(H.channels,2)
        if any(ch==RjCh), continue; end %preskocim rejectovane kanaly - kamil 14.6.2016          
        prevCh_shank = extractFromString(H.channels(chprev).name, 'string');
        currCh_shank = extractFromString(H.channels(ch).name, 'string');        
        if strcmp(prevCh_shank, currCh_shank) && strcmp(H.channels(chprev).signalType, 'SEEG') && strcmp(H.channels(ch).signalType, 'SEEG')
            chGroups{grp} = [chprev, ch]; %#ok<AGROW>
            grp = grp+1;
        end
        chprev = ch; %predchozi neni rejectovany kanal
    end   
end
