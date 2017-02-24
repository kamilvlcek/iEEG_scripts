function [ HH ] = HheaderCreate( pacID, nick,channelsfile,triggerCh )
%HHEADERCREATE vytvori header ve formatu Hammer Freiburg
%   pacID = id pacienta, nick - zkracene prijmeni, channelsfile - soubor se jmeny kanalu

HH = struct;
HH.patientTag = pacID;
HH.patientNick = nick;
HH.channels = struct;
assert(exist(channelsfile,'file')==2,'soubor se jmeny kanalu neexistuje');
channelnames = importdata(channelsfile,','); %file, kde na kazdem radku je jedno jmeno kanalu napr A1 
chnum = size(channelnames,1);
if ~exist('triggerCh','var')
    triggerCh = chnum - 2; %defaultni cislo trigerovaciho kanalu
end
for ch = 1:chnum  
    tmp = regexp(channelnames{ch,1},'([^,:]*)','tokens'); %na radce v channels muze byt jako druhe jmeno struktury oddelene carkou
    names = cat(2,tmp{:});
    HH.channels(ch).name=names{1};
    HH.channels(ch).numberOnAmplifier=ch;
    if ch==triggerCh
        signalType = 'triggerCh'; %vetsinou 2 pred koncem, kdyztak opravim rucne
    elseif chnum-ch < 2
        signalType = 'ECG'; %EKG
    else
        signalType = 'SEEG'; %iEEG channel
    end
    HH.channels(ch).signalType=signalType; 
    HH.channels(ch).ass_brainAtlas='n.a.';
    if ch <= 64
        HH.channels(ch).headboxNumber=1; %tohle pripadne muze byt spatne
    else 
        HH.channels(ch).headboxNumber=2;
    end
    if(numel(names)>1)
        HH.channels(ch).neurologyLabel=names{2};
    else
        HH.channels(ch).neurologyLabel='n.a.';
    end
end

end

