%% skript na analyzu PPA data z pacienta p73 - Kamil 16.8.2016
%podle manualu CiEEGData matlab scripts.pdf

%tenhle adresar je treba nastavit na kazdem pocitaci jinak
dir = 'D:\eeg\motol\pacienti\p073 Pech VT6\'; %nastaveni adresare, kde data na pocitaci jsou

%% 1 vytvoøení iEEG tøídy ze zdrojových dat
load([dir 'VT6_INV Test Vlcek1_X_ppa.mat']);

%% 1b vytvoøím objekt
E = CiEEGData(d,tabs,fs,[],H);

%% 2 Naètení headeru od Jirky Hammera
load([dir 'p73_header_kamil.mat']);
E.GetHHeader(H);

%% 2b Naètení epileptických událostí 
load([dir 'p73_ppa_epievents.mat']);
E.GetEpiEvents(DE);

%% 3 Zobrazení dat a vyøazení špatných kanálù
E.PlotElectrode();
E.RejectChannels([47 68]);

%% 4 Extrakce epoch a jejich manuální vyøazení
close all; %obrazek zase zavru
load([ dir 'p73_ppa.mat']);
E.ExtractEpochs(ppa,[-0.3 0.8]);
E.PlotElectrode();

%% 4b. TED projdu soubor a vyradim spatne epochy pomoci klavesy Delete
E.Save([dir 'ppa CiEEGData.mat']);
RjEpoch = E.RjEpoch; %ulozim vyrazene epochy pro pozdejsi pouziti
save([dir 'ppa_RjEpoch.mat'],'RjEpoch','-v7.3');

%% 5. Vykreslení úspìšnosti v testu
E.PlotResponses();

%% 6. Vytvoøení tøídy Chilbert ze zdrojových dat
clear E; %uvolnim misto v pameti
E = CHilbert(d,tabs,fs); 

%% 7. Zmìna reference na bipolární
% Nejdøív naètení headeru od Jirky Hammera (který už je v promìnné H),  vyøazení stejných kanálù 
E.GetHHeader(H);
E.RejectChannels([47 68]);
E.ChangeReference('b'); %vlastni zmena reference
E.Save([ dir 'ppa CHilbert.mat']);

%% 8. Vypoèet frekvenèních pásem
%Tohle mì trvalo asi 5 minut, ale u delšího souboru bude dýl
E.PasmoFrekvence([50:10:150]); %frekvencni pasma Hilbertovy obalky 50-150Hz po 10ti Hz

%% 9. Epochování dat a vyøazení tìch døíve vyøazených manuálnì
E.RejectEpochs(RjEpoch); %the epoching should follow loading excluded epochs, to correctly exclude them in CHilbert.HFreq
E.ExtractEpochs(ppa,[-0.3 0.8]); %epocha -0.3 az 0.8s kolem podnetu
E.Save(); %a uložení souboru

%% 10. Výpoèet statistických rozdílù
E.PsyData.Categories(1);
%melo by vypsat 
%0: Ovoce
%1: Scene
%2: Face
%3: Object

E.ResponseSearch(0.1,[2 3 1]); %0.1s window, rozdily vuci baseline a mezi kategoriemi, poradi kategorii Face Object Scene 
E.Save();

%% 11. Vykreslení prùmìrù síly frekvenèních pásem a rozdílù mezi kategoriemi
E.PlotResponseCh(20,[],1); %vykreslim kanal 20, vcetne okamzite hodnoty p

