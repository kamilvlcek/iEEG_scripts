% vypocet a zobrazeni Hilberovy obalky podle Jirky Hammera
% 18.2.2015
% viz tez d:\prace\programovani\matlab\kareljezek\hilbertovaobalka.m
close all;
pasma = 5:10:145;
castd = 1:10000+10000;

channel = 3;
srate =  1000;
t(2:end) = t(2:end)-t(1);

figure('Name','Spectrogram');
spectrogram(d(castd,channel),256,128,pasma,1000,'yaxis');

figure('Name','Hilbertovy obalky');

hilb = zeros(numel(castd),numel(pasma)); %obalka bude 

for pasmo = 1:numel(pasma)
    loF = pasma(pasmo);
    hiF = loF+10;
    hh  = hilbertJirka(d(castd,channel),loF,hiF,srate);
    hh = (hh./mean(hh)).*100;
    %plot(pasmo,mean(hilb(castd,pasmo)),'o');
    %plot(t(castd,1),hh);
    %ylim([0 1E4]);
    %hold on;
    hilb(castd,pasmo) = hh;
end
imagesc(1:0.001:10,pasma,hilb');
axis xy;