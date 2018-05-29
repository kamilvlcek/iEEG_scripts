testname = 'aedist';
reference = 'refOrig'; %'refBipo';
%% HLEDANI V HEADERECH
CB = CBrainPlot; %vytvorim tridu
%% najdu v headerech strukturu 
PAC = CB.StructFind({},{'PHG','parah','ent','subi'},testname,iff(strcmp(reference,'refBipo'),'b',[])); %#ok<NASGU> %hledam parahipp a entorhinal gyrus + subiculum
        % pouzivam pouze Martinovy zkratky
        
%ted pole PAC projdu a vymazu radky, ktere tam nepatri  (protoze muze nazev struktury obsahovat napr ent aj.     
%% nacteni drive ulozeneho seznamu
% NEBO - strukturu PAC si muzu zkopirovat do excelu (vcetne nazvu sloupcu) a pak si ji takhle znovu nacist
xlsfile = 'd:\eeg\motol\pacienti\0sumarne\structfind_mat.xlsx';
PAC = CB.StructFindLoad(xlsfile,2); 


%% PRACE S EXTRAKTY 
CM = CHilbertMulti; %vytvorim tridu
frekvence = '50-150';
label = 'PHGent'; % nazev exktraktu, pripoji s filename
datum = '2018-05-16';

%% vytvoreni extraktu s vybranymi kanaly pro kazdeho pacienta
filename = ['AEdist CHilbert ' frekvence ' -0.5-1.2 ' reference ' Ep2018-04_CHilb.mat']; %nazev souboru CHilbert, ze kterych se maji delat extrakty
overwrite = 1; %0= no overwrite - existujici soubory to preskoci 
filenames = CM.ExtractData(PAC,testname,filename,label, overwrite); %#ok<NASGU> 

%% nalezeni existujicich extraktu 
%NEBO pokud vim, ze mam vsechny extrakty vytvorene, muzu pouzit tohle
filenames = CM.FindExtract(testname,label, filename); 

%% zkontroluju extrakty
FILES = CM.TestExtract(filenames); 
%v promenne FILES v druhem sloupci bude pro kazdy soubor jeho velikost aj udaje, ktere se musi shodovat mezi soubory 
% - pocet vzorku, pocet epoch, vzorkovaci frekvence
%ty ktere se neshoduji, je nutne vyradit z promenne filenames

%% naimportuju extrakty
CM.ImportExtract(filenames);

%% smazani tridy
%po neuspesnem importu a pred dalsim importem je treba data tridy smazat
CM.Clear();

%% vypocitam si statistiku
setup = eval(['setup_' testname]); %nactu nastaveni
stat = setup.stat_kats; %resp setup.stat_kats{0} {1} nebo {2} pro menrot
CM.ResponseSearch(0.1,stat);

%% vyslednou sumarni tridu si ulozim, podobne jako kdyz ukladam data tridy CHilbert
CM.Save(['d:\eeg\motol\pacienti\0sumarne\' label ' ' frekvence ' ' reference ' ' datum]);

%% sumarni graf odpovedi
CM.IntervalyResp(); %graf velikosti odpovedi pres vsechny kanaly
CM.PlotResponseCh(); %odpoved pro kazdy kanal zvlast 


%% nactu starsi data
CM.Load(['d:\eeg\motol\pacienti\0sumarne\' label ' ' frekvence ' ' reference ' ' datum '_CHilb.mat']);