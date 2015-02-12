function [ M ] = means( ALLEEG,dataset,interval,prestimul,chtodisp )
%MEANS spocita prumery z intervalu hodnot (interval) pres vsechny opakovani
%   Detailed explanation goes here

srate = ALLEEG(dataset).srate;
%prestimul = 500; %ms
interval = interval + prestimul ; %celkovy cas od zacatku
iM = round(interval/1000*srate); % z ms na sekundy a 

iM = iM(1):iM(2);
data = ALLEEG(dataset).data;
M = zeros(size(data,1),1);
for ch = 1:size(data,1)
    M(ch)= mean(mean(squeeze(data(ch,iM,:))));
    if nargin>4 && ch==chtodisp
        X = squeeze(data(ch,:,:));
        XX = mean(X,2);
        figure;
        plot(XX);
        line([prestimul/1000*srate prestimul/1000*srate],[min(XX),max(XX)]);
    end
    
end

end

