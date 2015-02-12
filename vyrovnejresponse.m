function [ data2 ] = vyrovnejresponse( data,latencies,srate,stimultime,figures )
%VYROVNEJRESPONS Summary of this function goes here
%   Detailed explanation goes here
% data: elektrody x time x epochs

if ~exist('figures','var') , figures = 0; end
if figures
    figure('Name','pred upravou')
    M = mean(data,3);
    imagesc(M(1:124,:));
    colorbar;
    caxis([-30,30]);
end
Istim = round(srate * stimultime);
caspred = 0.6; %sekundy pred odpovedi vzit
caspo = 0.9; %sekundy po odpovedi vzit
n_electrodes = size(data,1); %pocet elektrod
n_epochs = size(data,3); %pocet epoch

delkaodpoved = round((caspo+caspred)*srate); %pocet bodu v matici kolem odpovedi
data2 = zeros(n_electrodes,Istim+delkaodpoved+1,n_epochs); %nova matice s daty

for ep = 1:n_epochs %cyklus pres epochy
    Iod =  round((latencies(ep)/1000-caspred + stimultime)*srate);
    Ido =  Iod + delkaodpoved;
    data2(:,:,ep)=data(:,[1:Istim,Iod:Ido],ep);
end

if figures
    %obrazek prumeru podnetu
    %figure();
    %M = mean(data2,3);
    %plot(M(126,:));

    %obrazek prumeru
    figure('Name','po uprave')
    M = mean(data2,3);
    imagesc(M(1:124,:));
    colorbar;
    caxis([-30,30]);
end
end

