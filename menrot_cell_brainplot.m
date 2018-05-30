testname = 'menrot'; %menrot, ppa
reference = 'refOrig'; %'refBipo', refHead, refEle
%% 1. HLEDANI V HEADERECH
CB = CBrainPlot; %vytvorim tridu
%% 1.1a najdu v headerech strukturu 
PAC = CB.StructFind({},{'FG','fug'},testname,iff(strcmp(reference,'refBipo'),'b',[])); %#ok<NASGU> %hledam parahipp a entorhinal gyrus + subiculum
        % pouzivam pouze Martinovy zkratky
        
%ted pole PAC projdu a vymazu radky, ktere tam nepatri  (protoze muze nazev struktury obsahovat napr ent aj.     
%TODO pridelat moznost vyjimek, FG + FuG toho najde strasne moc
%TODO prida skript na hledani preklepu - neexistujicich zkratek
%TODO co dela s pacientama z WIFI, ktery maji fs=500 hz?
%% 1.1b nacteni drive ulozeneho seznamu
% NEBO - strukturu PAC si muzu zkopirovat do excelu (vcetne nazvu sloupcu) a pak si ji takhle znovu nacist
xlsfile = 'd:\eeg\motol\pacienti\0sumarne\menrot_brain.xlsx'; %viz structfind_mat.xlsx na drive
PAC = CB.StructFindLoad(xlsfile,1);  %1 je cislo sheetu


%% 2. PRACE S EXTRAKTY 
CM = CHilbertMulti; %vytvorim tridu
frekvence = '50-150'; %15-31 
label = 'FG'; % nazev exktraktu, pripoji s filename
datum = '2018-05-30'; %dnesni datum - tak se pojmenuju vystupni souhrnny soubor
datumEP = '2018-01'; %datum v nazvu nacitaneho souboru
epochtime = '-1.0-1.5';

%TODO jak to je s vyrazenymi epochami u CHilbertMulti - funguji spravne?
%TODO revize starych rjepochch - u starych pacienti nejsou
%% 2.1a vytvoreni extraktu s vybranymi kanaly pro kazdeho pacienta
filename = ['Menrot CHilbert ' frekvence ' ' epochtime ' ' reference ' Ep' datumEP '_CHilb.mat']; %nazev souboru CHilbert, ze kterych se maji delat extrakty
overwrite = 0; %0= no overwrite - existujici soubory to preskoci 
filenames = CM.ExtractData(PAC,testname,filename,label, overwrite); %#ok<NASGU> 

%% 2.1b nalezeni existujicich extraktu 
%NEBO pokud vim, ze mam vsechny extrakty vytvorene, muzu pouzit tohle
filenames = CM.FindExtract(testname,label, filename); 

%% 2.2 zkontroluju extrakty
FILES = CM.TestExtract(filenames); 
%v promenne FILES v druhem sloupci bude pro kazdy soubor jeho velikost aj udaje, ktere se musi shodovat mezi soubory 
% - pocet vzorku, pocet epoch, vzorkovaci frekvence
%ty ktere se neshoduji, je nutne vyradit z promenne filenames

%% 2.3 naimportuju extrakty
CM.ImportExtract(filenames);

%% 2.3b smazani tridy
%po neuspesnem importu a pred dalsim importem je treba data tridy smazat
CM.Clear();

%% 2.4 vypocitam si statistiku
setup = eval(['setup_' testname]); %nactu nastaveni
kontrast = 1; %3=PT vs Ego 2=vy vs znacka 1=vsechno proti vsemu
stat = setup.stat_kats{kontrast}; %resp setup.stat_kats{1} {2} nebo {3} pro menrot 
CM.ResponseSearch(0.1,stat); 
%TODO - pouziti ruznych kontrastu a jejich ulozeni do jednoho souboru

%% 2.5 vyslednou sumarni tridu si ulozim, podobne jako kdyz ukladam data tridy CHilbert
CM.Save(['d:\eeg\motol\pacienti\0sumarne\CM ' testname ' ' label ' ' frekvence ' ' reference ' ' datum ' ' epochtime ' ' datumEP]);

%% 2.6 sumarni graf odpovedi
CM.IntervalyResp(); %graf velikosti odpovedi pres vsechny kanaly
%TODO - rozdily v kategoriich=podminkach - nesedi velikosti - mozna se pocitaji z celeho casu epochy?
CM.PlotResponseCh(); %odpoved pro kazdy kanal zvlast 
%TODO - zobrazeni reakcnich casu pro kazdeho pacienta i v sumarnim souboru
%% 2.6b nactu starsi data
CM.Load(['d:\eeg\motol\pacienti\0sumarne\CM ' testname ' ' label ' ' frekvence ' ' reference ' ' datum ' ' epochtime ' ' datumEP '_CHilb.mat']);

%% 3 OBRAZEK MOZKU
BPD = CM.ExtractBrainPlotData(); %vytvori data pro import do CBrainPlot
CB.ImportData(BPD); %naimportuje data z CHilbertMulti
CB.PlotBrain3D(); %vykresli obrazek mozku
%TODO funkce na zobrazeni cisteho plotu MNI souradnic - pridat velikost EEG odpovedi