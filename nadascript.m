%% MORLET WAVELETS PACIENT VT10

%clear all;

pacienti = pacienti_aedist();

pacient_cislo = 'p083 Kol VT10';
pacient = pacienti(4);

%channels = [58 59 60 61 66 67];
channels = [58 59];
frequencies = [1:0.5:4];

load(['d:\eeg\motol\pacienti\', pacient_cislo, '\', pacient.data]);

E = CMorlet(d,tabs,fs);
decimatefactor = 8 ;

load(['d:\eeg\motol\pacienti\', pacient_cislo, '\', pacient.header]);

E.GetHHeader(H);
E.RejectChannels(pacient.rjch);
E.ChangeReference('b');

E.PasmoFrekvence(frequencies,channels,0,decimatefactor);

load(['d:\eeg\motol\pacienti\', pacient_cislo, '\aedist\', pacient.psychopy]);
E.ExtractEpochs(aedist,[-0.2 1.2],[-0.5 -0.2]);

load(['d:\eeg\motol\pacienti\', pacient_cislo, '\aedist\', pacient.rjepoch]);
E.RejectEpochs(RjEpoch, RjEpochCh);

E.ResponseSearch(0.1,[0 1 2]);

%E.plotNada.channels = channels;


%% PLOT RESPONSE

for ch = channels
   
   display(['writing for ', num2str(ch)])
   E.PlotResponseCh(ch);
   %saveas(E.plotRCh.fh, fullfile('C:\Users\Nada\Documents\MATLAB\projects\CiEEG\iEEG_scripts\VT10', ['[1-4]Hz_', num2str(ch), 'ch_', E.CH.H.channels(ch).name]),'jpeg')

end

%%
for ch = 1:length(channels)  
   
   correct = aedist.data(find(aedist.data(:,3)),:);

   data = squeeze(E.phases(:,channels(ch), find(aedist.data(:,3)),:));

   fig = figure()
   colormap(hot)

   % PRE 0 = CERVENA
   subplot(1,3,1)

   data0 = data(:,find(correct(:,7)==0),:);
   % result = frekvencie x cas
   result = zeros(size(data0,3),size(data0,1));

   % cez vsetky frekvencie
   for i = 1:size(data0,3)
       result(i,:) = abs(mean(exp(1i*squeeze(data0(:,:,i))),2));      
   end

   T = linspace(E.epochtime(1),E.epochtime(2),size(E.d,1));
   imagesc(result, 'XData', T, 'YData', frequencies); colorbar
   xlabel('Time'); ylabel('Frequency'); title('cervena')

   % PRE 1 = EGO
   subplot(1,3,2)

   data1 = data(:,find(correct(:,7)==1),:);
   result = zeros(size(data1,3),size(data1,1));

   % cez vsetky frekvencie
   for i = 1:size(data1,3)
       result(i,:) = abs(mean(exp(1i*squeeze(data1(:,:,i))),2));  
   end

   imagesc(result, 'XData', T, 'YData', frequencies); colorbar
   xlabel('Time'); ylabel('Frequency'); title('ego')

   % PRE 2 = ALLO
   subplot(1,3,3)

   data2 = data(:,correct(:,7)==2,:);