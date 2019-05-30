%% 0. ZAKLADNI NASTAVENI
testname = 'menrot'; %menrot, ppa
reference = 'refBipo'; %'refBipo', refHead, refEle

%pokud uz mam ulozeny CHilbertMulti, staci spustit body 0., 2. a pak v bodu 2.7 soubor nactu
%% 1. HLEDANI V HEADERECH
CB = CBrainPlot; %vytvorim tridu
%% 1.1a najdu v headerech strukturu 
labelstofind = {'PHG','parah','ent','subi'}; %ktere neurology label chci zahrnout - staci cast (entorhinal), bez ohledu na velikost pismen
labelsNOTfind = {}; %ktere neurology label nechci zahrnout 
PacTodo = 1; %jestli se maji brat jen pacienti s todo=1
PAC = CB.StructFind({},labelstofind,testname,iff(strcmp(reference,'refBipo'),'b',[]),labelsNOTfind,PacTodo); %hledam parahipp a entorhinal gyrus + subiculum
        % pouzivam pouze Martinovy zkratky
        
%ted pole PAC projdu a vymazu radky, ktere tam nepatri  (protoze muze nazev struktury obsahovat napr ent aj.     

%% 1.1b nacteni drive ulozeneho seznamu
% NEBO - strukturu PAC si muzu zkopirovat do excelu (vcetne nazvu sloupcu) a pak si ji takhle znovu nacist
%xlsfile = 'd:\eeg\motol\pacienti\0sumarne\menrot_brain.xlsx'; %viz structfind_mat.xlsx na drive
%PAC = CB.StructFindLoad(xlsfile,1);  %1 je cislo sheetu


%% 2. PRACE S EXTRAKTY 
CM = CHilbertMulti; %vytvorim tridu
%% nastavim promenne
datum = '2018-05-30'; %dnesni datum - tak se pojmenuju vystupni souhrnny soubor
label = 'PHG-ent'; % nazev exktraktu, pripoji s filename

%podle nasledujicich tri promennych se sestavi jmeno nacitaneho souboru, napriklad Menrot CHilbert 50-150Hz -1.0-1.5 refBipo Ep2018-01_CHilb.mat
frekvence = '50-150Hz'; %15-31 
datumEP = '2018-01'; %datum v nazvu nacitaneho souboru
epochtime = '-1.0-1.5'; 
M = containers.Map({'ppa','menrot','aedist'},{'PPA','Menrot','Aedist'});
filename = [M(testname) ' CHilbert ' frekvence ' ' epochtime ' ' reference ' Ep' datumEP '_CHilb.mat']; %nazev souboru CHilbert, ze kterych se maji delat extrakty

%% 2.1a vytvoreni extraktu s vybranymi kanaly pro kazdeho pacienta
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
CMlabel = [M(testname) '_' label '_' frekvence]; 
CM.ImportExtract(filenames,CMlabel);

%% 2.3b smazani tridy
%po neuspesnem importu a pred dalsim importem je treba data tridy smazat
%kdyz hlasi import chybu, zkus tohle a pak import znova
CM.Clear();

%% 2.4 vypocitam si statistiku
setup = eval(['setup_' testname]); %nactu nastaveni
signum = 1; %1=chci jen odpovedi vyssi nez baseline, nebo druha kategorie, 0= vsechny signif rozdily
CM.ResponseSearchMulti(0.1,setup.stat_kats);
%nastavim oznaceni kanalu podle signifikance odpovedi
CM.SetStatActive(2); %pridam k prvnimu oznaceni - cislo statistiky podle setup.stat_kats
CM.SelChannelStat({3},1,0,signum); % {'znacka-2D+znacka-3D X vy-2D+vy-3D'} = k
CM.SetStatActive(3); %jeste pridam k oznaceni   - cislo statistiky podle setup.stat_kats
CM.SelChannelStat({3},2,1,signum); % {'vy-3D+znacka-3D X vy-2D+znacka-2D'} = l
CM.SetStatActive(1); %nove oznaceni
CM.SelChannelStat({1, 2, 3, 4},[3 4 5 6],1,signum); % 'vy-2D' 'vy-3D' 'znacka-2D' 'znacka-3D' = fghj

%% 2.5 vyslednou sumarni tridu si ulozim, podobne jako kdyz ukladam data tridy CHilbert
CM.Save([ setup.basedir '0sumarne\CM ' M(testname) ' ' label ' ' frekvence ' ' reference ' ' epochtime ' Ep' datumEP ' '  datum]);
%TODO folder

%% 2.6 sumarni graf odpovedi
CM.IntervalyResp(); %graf velikosti odpovedi pres vsechny kanaly
%TODO - rozdily v kategoriich=podminkach - nesedi velikosti - mozna se pocitaji z celeho casu epochy?
CM.PlotResponseCh(); %odpoved pro kazdy kanal zvlast 
%TODO - zobrazeni reakcnich casu pro kazdeho pacienta i v sumarnim souboru
%% 2.6b nactu starsi data
CM.Load([ setup.basedir '0sumarne\CM ' M(testname) ' ' label ' ' frekvence ' ' reference ' ' epochtime ' Ep' datumEP ' '  datum '_CHilb.mat']);

%% 3 OBRAZEK MOZKU 
%  3.1 nastavim statistiku
statToDo = [2 3]; %ze kterych statistiky chci udelat obrazky, 2=vy vs znacka
stat = statToDo(1); %tohle musim rucne editova - nejdriv 1, pak 2. Muzu take kliknout na cislo a pravym tlacitkem mysi zvolit Increment Value and Run Section
CM.SetStatActive(stat); 
%% 3.2 nactu data
BPD = CM.ExtractBrainPlotData(); %vytvori data pro import do CBrainPlot
CB.ImportData(BPD); %naimportuje data z CHilbertMulti
%% 3.3 - vyrobim obrazky
CB.PlotBrain3DConfig(struct('Names',0,'NoNames',0,'signum',signum,'overwrite',1)); %nahraje defaultni konfiguraci
CB.PlotBrain3D(); %vykresli obrazek mozku
%TODO funkce na zobrazeni cisteho plotu MNI souradnic - pridat velikost EEG odpovedi