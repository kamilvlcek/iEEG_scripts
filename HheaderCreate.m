function [ HH ] = HheaderCreate( pacID, nick,channelsfile,triggerCh )
%HHEADERCREATE vytvori header ve formatu Hammer Freiburg
%   pacID = id pacienta, nick - zkracene prijmeni, channelsfile - soubor se jmeny kanalu

HH = struct;
HH.patientTag = pacID; %tohle jirka v headeru asi nema
HH.patientNick = nick; %tohle jirka v headeru asi nema
HH.subjName = [pacID ' ' nick]; % to je i headeru od Jirky
HH.channels = struct;
assert(exist(channelsfile,'file')==2,'soubor se jmeny kanalu neexistuje');
channelnames = importdata(channelsfile,','); %file, kde na kazdem radku je jedno jmeno kanalu napr A1 
chnum = size(channelnames,1);
if ~exist('triggerCh','var')
    triggerCh = chnum - 2; %defaultni cislo trigerovaciho kanalu
end
fprintf('last channels: ');
for ch = 1:chnum  
    tmp = regexp(channelnames{ch,1},'([^,:]*)','tokens'); %na radce v channels muze byt jako druhe jmeno struktury oddelene carkou
    names = cat(2,tmp{:});
    HH.channels(ch).name=names{1};
    HH.channels(ch).numberOnAmplifier=ch;
    if ch==triggerCh
        signalType = 'triggerCh'; %vetsinou 2 pred koncem, kdyztak opravim rucne
    elseif chnum-ch < 2 && chnum % >= triggerCh %26.4.2017 - trigerovaci kanal nikdy neni az za EKG, to spis neni zadne EKG
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
    if ch >= chnum -5
        fprintf('%i:%s,',ch,names{1});
    end    
end
fprintf('\n');

end

