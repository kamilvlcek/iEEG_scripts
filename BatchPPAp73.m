%% skript na analyzu PPA data z pacienta p73 - Kamil 16.8.2016
%podle manualu CiEEGData matlab scripts.pdf

%tenhle adresar je treba nastavit na kazdem pocitaci jinak
dir = 'D:\eeg\motol\pacienti\p073 Pech VT6\'; %nastaveni adresare, kde data na pocitaci jsou

%% 1 vytvo�en� iEEG t��dy ze zdrojov�ch dat
load([dir 'VT6_INV Test Vlcek1_X_ppa.mat']);

%% 1b vytvo��m objekt
E = CiEEGData(d,tabs,fs,[],H);

%% 2 Na�ten� headeru od Jirky Hammera
load([dir 'p73_header_kamil.mat']);
E.GetHHeader(H);

%% 2b Na�ten� epileptick�ch ud�lost� 
load([dir 'p73_ppa_epievents.mat']);
E.GetEpiEvents(DE);

%% 3 Zobrazen� dat a vy�azen� �patn�ch kan�l�
E.PlotElectrode();
E.RejectChannels([47 68]);

%% 4 Extrakce epoch a jejich manu�ln� vy�azen�
close all; %obrazek zase zavru
load([ dir 'p73_ppa.mat']);
E.ExtractEpochs(ppa,[-0.3 0.8]);
E.PlotElectrode();

%% 4b. TED projdu soubor a vyradim spatne epochy pomoci klavesy Delete
E.Save([dir 'ppa CiEEGData.mat']);
RjEpoch = E.RjEpoch; %ulozim vyrazene epochy pro pozdejsi pouziti
save([dir 'ppa_RjEpoch.mat'],'RjEpoch','-v7.3');

%% 5. Vykreslen� �sp�nosti v testu
E.PlotResponses();

%% 6. Vytvo�en� t��dy Chilbert ze zdrojov�ch dat
clear E; %uvolnim misto v pameti
E = CHilbert(d,tabs,fs); 

%% 7. Zm�na reference na bipol�rn�
% Nejd��v na�ten� headeru od Jirky Hammera (kter� u� je v prom�nn� H),  vy�azen� stejn�ch kan�l� 
E.GetHHeader(H);
E.RejectChannels([47 68]);
E.ChangeReference('b'); %vlastni zmena reference
E.Save([ dir 'ppa CHilbert.mat']);

%% 8. Vypo�et frekven�n�ch p�sem
%Tohle m� trvalo asi 5 minut, ale u del��ho souboru bude d�l
E.PasmoFrekvence([50:10:150]); %frekvencni pasma Hilbertovy obalky 50-150Hz po 10ti Hz

%% 9. Epochov�n� dat a vy�azen� t�ch d��ve vy�azen�ch manu�ln�
E.RejectEpochs(RjEpoch); %the epoching should follow loading excluded epochs, to correctly exclude them in CHilbert.HFreq
E.ExtractEpochs(ppa,[-0.3 0.8]); %epocha -0.3 az 0.8s kolem podnetu
E.Save(); %a ulo�en� souboru

%% 10. V�po�et statistick�ch rozd�l�
E.PsyData.Categories(1);
%melo by vypsat 
%0: Ovoce
%1: Scene
%2: Face
%3: Object

E.ResponseSearch(0.1,[2 3 1]); %0.1s window, rozdily vuci baseline a mezi kategoriemi, poradi kategorii Face Object Scene 
E.Save();

%% 11. Vykreslen� pr�m�r� s�ly frekven�n�ch p�sem a rozd�l� mezi kategoriemi
E.PlotResponseCh(20,[],1); %vykreslim kanal 20, vcetne okamzite hodnoty p

