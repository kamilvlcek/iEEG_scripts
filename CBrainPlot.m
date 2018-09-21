classdef CBrainPlot < matlab.mixin.Copyable
    %CBRAINPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        VALS; %souhrnne prumery elektrod pres vsechny pacienty - cell(intervaly x kategorie)
        MNI;  %souhrnna MNI data pres vsechny pacienty - cell(intervaly x kategorie)
        NAMES; %souhrnna jmena elektrod pres vsechny pacienty - cell(intervaly x kategorie)
        NLabels; %jmena neurologyLabels z Headeru
        intervals; % intervaly z funkce IntervalyResp
        katstr; %jmena kategorii
        brainsurface; %ulozeny isosurface z main_brainPlot
        testname; %jmeno zpracovavaneho testu
        katstr_pacients;
        numelP; %pocty  signif elektrod pro kazdy pacient x interval x kategorie
        pacients;    
        filename; %jmeno zpracovavanych souboru
        PAC; %struktura stejna jako vraci funkce StructFind, naplni se v IntervalyResp
        iPAC; %index v poli PAC
        reference; %reference
        Hf; %seznam frekvencnich pasem
        selCh; %vyber kanalu, ktere zobrazit
        plotBrain3Dcfg; %struktura s nastavenim plotBrain3D
    end
    
    methods (Access = public)        
        function [obj] = IntervalyResp(obj,testname,intervals,filename,contrast,signum)
            %IntervalyResp(testname,intervals,filename,contrast)
            %vola postupne pro vsechny pacienty E.IntervalyResp a uklada vysledky
            %vyradi vsechny kontakty bez odpovedi nebo se zapornou odpovedi
            %spoji vsechno dohromady
            %vrati vysledky ve formatu pro SEEE-vizualization
            %napr CB.IntervalyResp('aedist',[0.2 0.8],'AEdist CHilbert 50-120 refBipo Ep2017-11_CHilb.mat');
            %signum = jestli chci jen kat1>kat2 (1), nebo obracene (-1), nebo vsechny (0)
            if ~exist('contrast','var'), contrast = 1; end; %defaultni je prvni kontrast            
            if strcmp(testname,'aedist')
                pacienti = pacienti_aedist(); %nactu celou strukturu pacientu    
            elseif strcmp(testname,'ppa')
                pacienti = pacienti_ppa(); %nactu celou strukturu pacientu    
            elseif strcmp(testname,'menrot')
                pacienti = pacienti_menrot(); %nactu celou strukturu pacientu    
            else
                error('neznamy typ testu');
            end
            obj.testname = testname;
            obj.intervals = intervals; 
            obj.filename = filename;
            elcount = []; %jen inicializace            
            P = {}; M = {}; N = {}; %jen inicializace
            obj.PAC = [];
            obj.pacients = cell(numel(pacienti),1); 
            obj.katstr_pacients = []; %musim to smazat, nize testuju, jestil to je prazdne
            obj.numelP = [];  %tam budu ukladat pocty elektrod pro kazdy pacient x interval x kategorie
            for p = 1:numel(pacienti) % cyklus pacienti
                if pacienti(p).todo 
                    disp(['***   ' pacienti(p).folder '   ***']);
                    E = pacient_load(pacienti(p).folder,testname,filename,[],[],[],0); %nejspis objekt CHilbert, pripadne i jiny; loadall = 0
                    if isempty(E)
                        disp('no data');
                        pacienti(p).todo = 0; %nechci ho dal zpracovavat
                        continue;
                    end
                    E.SetStatActive(contrast); %nastavi jeden z ulozenych statistickych kontrastu
                    [prumery, MNI,names,~,katstr,neurologyLabels] = E.IntervalyResp( intervals,[],signum,0);   %#ok<PROPLC> %no figure, funkce z CiEEGData                           
                    obj.pacients{p} = pacienti(p).folder;
                    obj.GetPAC(prumery,E.CH.H,pacienti(p).folder);
                    obj.reference = E.reference;
                    if isprop(E,'Hf')
                        obj.Hf = E.Hf;
                    end
                    clear E;
                    if isempty(obj.katstr_pacients)
                        obj.katstr = [katstr 'AllEl']; %#ok<PROPLC> %katstr se ziskava z IntervalyResp
                        obj.katstr_pacients = cell(numel(pacienti),numel(katstr)); %#ok<PROPLC>
                        obj.numelP = zeros(numel(pacienti),size(intervals,1),numel(katstr)+1); %#ok<PROPLC> %tam budu ukladat pocty elektrod pro kazdy interval a pacienta
                        elcount = zeros(size(prumery,2),size(prumery,3)+1); %pocet elektrod pro kazdy casovy interval a kategorii - interval x kategorie
                        P = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne prumery pro vsechny pacienty: pacient*interval*kategorie
                        M = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne MNI koordinaty pro vsechny pacienty
                        N = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne names pro vsechny pacienty
                        NL = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne neurologyLabels
                            %+1 je pro obrazek vsech elektrod i tech bez odpovedi
                    end
                    obj.katstr_pacients(p,:) = katstr; %#ok<PROPLC> %jsou kategorie u vsech pacientu ve stejnem poradi?
                    for interval = 1:size(prumery,2) % cyklus intervaly
                        for kat = 1:size(prumery,3)+1 % cyklus kategorie podnetu
                            if kat <= size(prumery,3) %obvykle kategorie
                                ip = prumery(:,interval, kat) ~= 0; % chci i zaporny rozdil ; aby tam neco bylo 
                                P{p,interval,kat}=prumery(ip,interval, kat); %#ok<AGROW>
                                M{p,interval,kat}=MNI(ip); %#ok<AGROW,PROPLC>>
                                N{p,interval,kat} = strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],sum(ip),1)),names(ip)); %#ok<AGROW>
                                NL{p,interval,kat}= strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],sum(ip),1)),neurologyLabels(ip)); %#ok<AGROW>
                                elcount(interval,kat) = elcount(interval,kat) + sum(ip); %#ok<AGROW>
                                obj.numelP(p,interval,kat)=sum(ip);
                            else %kategorie jakoby navic pro vykresleni jen pozice elekrod
                                channels = size(prumery,1);
                                P{p,interval,kat}=zeros(channels,1); %#ok<AGROW> % 0 pro kazdy kanal - vsechny stejnou barvou
                                M{p,interval,kat}=MNI; %#ok<AGROW,PROPLC>>
                                N{p,interval,kat} = strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],channels,1)),names); %#ok<AGROW>
                                NL{p,interval,kat}= strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],channels,1)),neurologyLabels); %#ok<AGROW>
                                elcount(interval,kat) = elcount(interval,kat) + channels; %#ok<AGROW>                                
                                obj.numelP(p,interval,kat)=channels;
                            end
                        end                       
                    end 
                end
            end
            %ted z P M a N rozdelenych po pacientech udelam souhrnna data
            obj.VALS = cell(size(elcount)); %souhrnne prumery - interval * kategorie
            obj.MNI = cell(size(elcount)); 
            obj.NAMES = cell(size(elcount)); 
            obj.NLabels = cell(size(elcount)); 
            obj.selCh = cell(size(elcount));  %prazdny vyber elektrod           
            if sum([pacienti.todo])>0 
                for interval = 1:size(prumery,2) 
                    for kat = 1:size(prumery,3)+1                   
                          obj.VALS{interval,kat} = zeros(elcount(interval,kat),1);
                          obj.MNI{interval,kat} = struct('MNI_x',{},'MNI_y',{},'MNI_z',{});
                          obj.NAMES{interval,kat}   = cell(elcount(interval,kat),1);
                          obj.NLabels{interval,kat} = cell(elcount(interval,kat),1);
                          iVALS = 1;
                          for p = 1:numel(pacienti) 
                              if pacienti(p).todo
                                  n = numel(P{p,interval,kat});
                                  obj.VALS{interval,kat} (iVALS:iVALS+n-1)=P{p,interval,kat};
                                  obj.MNI{interval,kat}  (iVALS:iVALS+n-1)=M{p,interval,kat};
                                  obj.NAMES{interval,kat}(iVALS:iVALS+n-1)  =N{p,interval,kat};
                                  obj.NLabels{interval,kat}(iVALS:iVALS+n-1)=NL{p,interval,kat};
                                  iVALS = iVALS + n;
                              end
                          end
                    end
                end             
                disp(''); %prazdna radka
                %disp(['vytvoreny ' num2str(numel(obj.katstr)) ' kategorie: ' cell2str(obj.katstr)]);
                %jeste vypisu pocty elektrod pro kazdou kategorii
                fprintf('\npocty elektrod v %i kategoriich (pro vsechny pacienty):\n',numel(obj.katstr));
                for kat = 1:numel(obj.katstr)
                    fprintf('%s:\t', obj.katstr{kat});
                    for int = 1:size(intervals,1)
                        fprintf(' %i,', sum(obj.numelP(:,int,kat)));
                    end
                    fprintf('\n');
                end               
            else
                disp('zadny soubor nenalezen');
            end
        end
        function [obj] = GetPAC(obj,prumery,H,pac_folder)
            %vytvori a ulozi PAC data do obj.PAC pro jednoho pacienta z prumery a headeru
            %updatuje obj.iPAC jako index v obj.PAC
            if isempty(obj.PAC)
                obj.PAC = cell(size(prumery,2),size(prumery,3));
                obj.iPAC = ones(size(prumery,2),size(prumery,3));
            end
            for interval = 1:size(prumery,2)
                for kat = 1:size(prumery,3)
                    if obj.iPAC(interval,kat) == 1
                        obj.PAC{interval,kat} = {};
                    end                    
                    index = find(prumery(:,interval,kat)~=0);
                    for ii = 1:numel(index)                
                        obj.PAC{interval,kat}(obj.iPAC(interval,kat)).pacient = pac_folder; %#ok<AGROW>
                        obj.PAC{interval,kat}(obj.iPAC(interval,kat)).ch = index(ii); %#ok<AGROW>
                        obj.PAC{interval,kat}(obj.iPAC(interval,kat)).name = H.channels(index(ii)).name; %#ok<AGROW>
                        obj.PAC{interval,kat}(obj.iPAC(interval,kat)).neurologyLabel = H.channels(index(ii)).neurologyLabel; %#ok<AGROW>
                        obj.PAC{interval,kat}(obj.iPAC(interval,kat)).ass_brainAtlas = H.channels(index(ii)).ass_brainAtlas;%#ok<AGROW>
                        obj.PAC{interval,kat}(obj.iPAC(interval,kat)).ass_cytoarchMap = H.channels(index(ii)).ass_cytoarchMap; %#ok<AGROW>
                        obj.iPAC(interval,kat) = obj.iPAC(interval,kat) + 1;
                    end
                end
            end
        end
        function obj = ImportData(obj,BPD)
            %vlozi data, ktera jsem vytvoril pomoci CHilbert.ExtractBrainPlotData
            obj.VALS = BPD.VALS;
            obj.MNI = BPD.MNI;
            obj.NAMES = BPD.NAMES;
            obj.NLabels = BPD.NLABELS;
            obj.katstr = BPD.katstr;
            obj.intervals = BPD.intervals;       
            obj.testname = BPD.testname;
            obj.reference = BPD.reference;
            obj.Hf = BPD.Hf;
            obj.selCh = BPD.selCh; %vyber kanalu z CHilbertMulti aj
            if isfield(BPD,'signum')
                obj.plotBrain3Dcfg.signum = BPD.signum;
            end
        end
        function obj = PlotBrain3DConfig(obj,cfg)
            if ~exist('cfg','var'), cfg = struct; end
            if ~isprop(obj,'plotBrain3Dcfg')
                obj.plotBrain3Dcfg = struct; 
            end
            if isstruct(cfg) && isfield(cfg,'signum')
                obj.plotBrain3Dcfg.signum = cfg.signum;
            elseif ~isfield(obj.plotBrain3Dcfg,'signum') || isempty(obj.plotBrain3Dcfg.signum)
                %vyplnuje se, jen pokud je zatim prazcne
                obj.plotBrain3Dcfg.signum = 0; %defaultni je rozdil kladny i zaporny
            end
            if isstruct(cfg) && isfield(cfg,'outputDir')
                obj.plotBrain3Dcfg.outputDir = cfg.outputDir;
            else
                obj.plotBrain3Dcfg.outputDir = 'd:\eeg\motol\CBrainPlot\';
            end
            if isstruct(cfg) && isfield(cfg,'overwrite')
                obj.plotBrain3Dcfg.overwrite = cfg.overwrite;
            else
                obj.plotBrain3Dcfg.overwrite = 1; %defaultne se vystupni soubory prepisuji
            end
            if isstruct(cfg) && isfield(cfg,'NLabels')
                obj.plotBrain3Dcfg.NLabels = cfg.NLabels;
            else
                obj.plotBrain3Dcfg.NLabels = 0; %defaultne se nevypisuji anatomicke lokalizace
            end
            if isstruct(cfg) && isfield(cfg,'NOnames')
                obj.plotBrain3Dcfg.NOnames = cfg.NOnames;
            else
                obj.plotBrain3Dcfg.NOnames = 0; %defaultne se nedelaji obrazky bez popisu elektrod
            
            end            
        end
        function PlotBrain3D(obj,kategorie)
            %vykresli jpg obrazky jednotlivych kategorii a kontrastu mezi nimi            
            %TODO do jmena vystupniho jpg pridat i frekvence a referenci, aby se to neprepisovalo
            %TODO je mozne ty signif vyexportovat a pak je nacist zase do CHilbertMulti?
            %TODO do vystupni tabulky nejak dostat anatomickou lokalizaci?
            assert(~isempty(obj.VALS),'zadna data VALS');
            plotSetup = {};
            if ~exist('kategorie','var') || isempty(kategorie) , kategorie = 1:size(obj.VALS,2); end %muzu chtit jen nektere kategorie
            signum = obj.plotBrain3Dcfg.signum; 
            plotSetup.outputDir = obj.plotBrain3Dcfg.outputDir;                                              
            
            if ~isempty(obj.brainsurface)
                brainsurface = obj.brainsurface;  %#ok<PROPLC>
            else
                brainsurface = []; %#ok<PROPLC>
            end
            hybernovat = 0; %jestli chci po konci skriptu pocitac uspat - ma prednost
            vypnout = 0;  %jestli chci po konci skriptu pocitac vypnout (a nechci ho hybernovat)             
            plotSetup.figureVisible = 'off';   %nechci zobrazovat obrazek 
            plotSetup.FontSize = 4; 
            plotSetup.myColorMap = iff(signum ~= 0,parula(128) ,jet(128));    %pokud jednostrane rozdily, chci parula
            %barevna skala od Nadi
            plotSetup.customColors.customColor = true; % custom colormap oddeli negativne a pozitivne hodnoty - 29.6.2018
            plotSetup.customColors.flip = 0; %pokud chci prehodit barvy
            plotSetup.customColors.darkneg = [50 145 0]; %tmave zelena
            plotSetup.customColors.lightneg = [212 255 171];
            plotSetup.customColors.lightpos = [246 203 203];
            plotSetup.customColors.darkpos = [162 2 2]; %tmave cervena
            plotSetup.figureNamePrefix = [ obj.testname '_' num2str(obj.Hf([1 end]),'%i-%iHz') '_' obj.reference '_names']; %default name

            tablelog = cell(obj.pocetcykluPlot3D(kategorie,signum)+2,7); % z toho bude vystupni xls tabulka s prehledem vysledku
            tablelog(1,:) = {datestr(now),obj.filename,'','','','',''}; %hlavicky xls tabulky
            tablelog(2,:) = {'interval','kategorie','chname','neurologyLabel','mni','val','selected'}; %hlavicky xls tabulky
            iTL = 2; %index v tablelog
            tic; %zadnu merit cas
            for interval = 1:size(obj.VALS,1) 
                for kat = kategorie
                    if signum > 0 
                        iV = obj.VALS{interval,kat} > 0; %jen kladne rozdily
                    elseif signum <0 
                        iV = obj.VALS{interval,kat} < 0; %jen zaporne rozdily
                    elseif ~isempty(obj.selCh{interval,kat})
                        iV = ismember(1:numel(obj.VALS{interval,kat}),find(any(obj.selCh{interval,kat},2))); %vyber kanalu k zobrazeni, napriklad z CHilbertMulti
                    else
                        iV = true(size(obj.VALS{interval,kat})); %vsechny rozdily
                    end                    
%                
                    katname = obj.katstr{kat};
                    plotSetup.circle_size = iff(strcmp(katname,'all') || strcmp(katname,'AllEl'),28,56); %mensi kulicka pokud vsechny elektrody                
                    
                    if strcmp(plotSetup.figureVisible,'off')
                        disp('figures invisible');
                    end
                    figureNameNames = [ obj.testname '_' num2str(obj.intervals(interval,:),'%.1f-%.1fs')  '_' katname '_' num2str(signum) ...
                            '_' num2str(obj.Hf([1 end]),'%i-%iHz') '_' obj.reference '_names'];
                    figureNameNoNames = [ obj.testname '_' num2str(obj.intervals(interval,:),'%.1f-%.1fs')  '_' katname '_' num2str(signum) ...
                            '_' num2str(obj.Hf([1 end]),'%i-%iHz') '_' obj.reference '_NOnames'];
                    if numel(obj.VALS{interval,kat}(iV)) > 0
                        
                        vals_channels = obj.VALS{interval,kat}(iV); %parametr  main_brainPlot
                        if signum ~= 0
                            vals_channels = vals_channels*signum; %u zapornych hodnot prehodim znamenko
                        end
                        mni_channels = obj.MNI{interval,kat}(iV);                                                                                                 
                        names_channels = iff(obj.plotBrain3Dcfg.NLabels, obj.NLabels{interval,kat}(iV), obj.NAMES{interval,kat}(iV));                        
                        
                        if ~strcmp(obj.katstr{kat},'AllEl') %nechci to pro kategorii vsech elektrod
                            for iVal = 1:numel(vals_channels)
                                tablelog(iVal + iTL,:) = { sprintf('[%.1f %.1f]',obj.intervals(interval,:)),obj.katstr{kat}, obj.NAMES{interval,kat}{iVal}, obj.NLabels{interval,kat}{iVal}, ...
                                    sprintf('[%.1f,%.1f,%.1f]',mni_channels(iVal).MNI_x, mni_channels(iVal).MNI_y, mni_channels(iVal).MNI_z), vals_channels(iVal),int8(iV(iVal))};
                            end
                            iTL = iTL + numel(vals_channels);
                        end
                        
                        %nejdriv vykreslim bez popisku elektrod
                        if obj.plotBrain3Dcfg.NOnames 
                            if  isempty(dir([ plotSetup.outputDir '3D_model\' figureNameNoNames '*'])) || obj.plotBrain3Dcfg.overwrite==1 
                                plotSetup.figureNamePrefix = figureNameNoNames;
                                disp(plotSetup.figureNamePrefix);                            
                                brainsurface = main_brainPlot(vals_channels,mni_channels,[],brainsurface,plotSetup);  %#ok<PROPLC>
                                %volam Jirkuv skript, vsechny ty promenne predtim jsou do nej
                                if isempty(obj.brainsurface)
                                    obj.brainsurface = brainsurface; %#ok<PROPLC> %ulozim si ho pro dalsi volani
                                end
                            else
                                disp(['soubor uz existuje ' figureNameNoNames ' - neprepisuju ']);
                            end
                        else
                             disp([figureNameNoNames ': obrazky noNames negeneruju' ]);
                        end
                        
                        %a pak jeste s popisy elektrod                        
                        if isempty(dir([ plotSetup.outputDir '3D_model\' figureNameNames '*'])) || obj.plotBrain3Dcfg.overwrite==1 
                            plotSetup.figureNamePrefix = figureNameNames;
                            disp(plotSetup.figureNamePrefix);                                                     
                            brainsurface = main_brainPlot(vals_channels,mni_channels,names_channels,brainsurface,plotSetup);    %#ok<PROPLC>  
                            if isempty(obj.brainsurface)
                                obj.brainsurface = brainsurface; %#ok<PROPLC> %ulozim si ho pro dalsi volani
                            end
                        else
                            disp(['soubor uz existuje ' figureNameNames ' - neprepisuju ']);
                        end
                    else  
                        disp(['zadne hodnoty pro ' plotSetup.figureNamePrefix ' - neukladam ']);                                         
                    end
                end
            end
            toc; %ukoncim mereni casu a vypisu
            logfilename = ['logs\PlotBrain3D_' figureNameNames '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') ];
            
            xlswrite([plotSetup.outputDir logfilename '.xls'],tablelog); %zapisu do xls tabulky
            if hybernovat
                system('shutdown -h')  %#ok<UNRCH>
            elseif vypnout            
                system('shutdown -s') %#ok<UNRCH>
            end
        end
    end
    methods (Static,Access = public)
        function PAC = StructFind(struktura,label,testname,reference,labelnot)
            %najde pacienty, jejich headery obsahuji mozkovou strukturu
            %struktura je nazev struktury podle atlas napriklad hippo, label je kratky nazev podle martina, napriklad hi
            if ~exist('label','var'),    label = struktura; end %defaultni test
            if ~exist('testname','var') || isempty(testname), testname = 'aedist'; end %defaultni test
            if ~exist('reference','var') || isempty(reference), reference = []; end %defaultni test
            if ~exist('labelnot','var'),    labelnot = {}; end %defaultni test
            if ischar(struktura), struktura = {struktura}; end %prevedu na cell array
            if ischar(label), label = {label}; end %prevedu na cell array
            [ pacienti, setup ] = pacienti_setup_load( testname );
            PAC = {};
            iPAC = 1;
            for p = 1:numel(pacienti)
                disp(['* ' pacienti(p).folder ' - ' pacienti(p).header ' *']);
                hfilename = [setup.basedir pacienti(p).folder '\' pacienti(p).header];                
                if exist(hfilename,'file')==2
                    load(hfilename);
                else
                    disp(['header ' hfilename ' neexistuje']);
                    continue; %zkusim dalsiho pacienta, abych vypsal, ktere vsechny headery neexistujou
                end               
                if ~isempty(reference)
                    CH = CHHeader(H);
                    CH.RejectChannels( pacienti(p).rjch); %musim vyradit vyrazene kanaly, protoze ty se vyrazuji v bipolarni referenci
                    CH.ChangeReference(reference); %nove od 18.1.2018
                    H = CH.H;
                end
                ii = ~cellfun(@isempty,{H.channels.neurologyLabel}); %neprazdne cells
                if ~isempty(struktura) || ~isempty(label)                    
                    labels = lower({H.channels(ii).neurologyLabel});
                    iLabels = contains(labels,lower(label)); %najde vsechny label naraz                    
                    if ~isempty(labelnot)
                        iLabelsNOT = contains(labels,lower(labelnot)); %najde vsechny label naraz  
                        iLabels = iLabels & ~iLabelsNOT;
                    end
                    index = find(iLabels);%bude obsahovat cisla vybranych kanalu - jeden radek
                    iiBA = ~cellfun(@isempty,{H.channels.ass_brainAtlas}); %neprazdne cells
                    iiCM = ~cellfun(@isempty,{H.channels.ass_cytoarchMap}); %neprazdne cells
                    for jj = 1:size(struktura,2)
                        index = [ index find(~cellfun('isempty',strfind(lower({H.channels(iiCM).ass_cytoarchMap}),lower(struktura{jj}))))]; %#ok<AGROW>
                        index = [ index find(~cellfun('isempty',strfind(lower({H.channels(iiBA).ass_brainAtlas}),lower(struktura{jj}))))]; %#ok<AGROW>
                    end
                    index = union(index,[]); %vsechny tri dohromady - seradi podle velikost a odstrani duplikaty
                else
                    index = 1:numel(H.channels); %pokud neuvedu nic k hledani, beru vsechny kanaly
                end
                if isempty(reference) || reference ~= 'b' %pokud jsem kanaly nevyradil uz pri zmene reference - vyrazuji se jen pri bipolarni
                    indexvyradit = ismember(index, pacienti(p).rjch); %vyrazene kanaly tady nechci
                    index(indexvyradit)=[]; 
                end
                
                %vrati indexy radku ze struct array, ktere obsahuji v sloupci neurologyLabel substring struktura
                for ii = 1:numel(index)                
                    PAC(iPAC).pacient = pacienti(p).folder; %#ok<AGROW>
                    PAC(iPAC).ch = index(ii); %#ok<AGROW>
                    PAC(iPAC).name = H.channels(index(ii)).name; %#ok<AGROW>
                    PAC(iPAC).neurologyLabel = H.channels(index(ii)).neurologyLabel; %#ok<AGROW>
                    PAC(iPAC).ass_brainAtlas = H.channels(index(ii)).ass_brainAtlas;%#ok<AGROW>
                    PAC(iPAC).ass_cytoarchMap = H.channels(index(ii)).ass_cytoarchMap; %#ok<AGROW>
                    PAC(iPAC).MNI_x = H.channels(index(ii)).MNI_x; %#ok<AGROW>
                    PAC(iPAC).MNI_y = H.channels(index(ii)).MNI_y; %#ok<AGROW>
                    PAC(iPAC).MNI_z = H.channels(index(ii)).MNI_z; %#ok<AGROW>
                    iPAC = iPAC + 1;
                end
            end 
            xlsname = ['./logs/StructFind PAC_' testname '_' cell2str(struktura,1) '_' cell2str(label,1) '.xlsx'];
            writetable(struct2table(PAC), xlsname); %ulozimvysledek taky do xls
        end
        function PAC = StructFindLoad(xlsfile,sheet)
            %nacteni struktury PAC z existujiciho xls souboru, napr po editaci radku            
             if ~exist('sheet','var'), sheet = 1; end %defaultni je prvni list
             [~ ,~ , raw]=xlsread(xlsfile,sheet); 
             for iraw = 1:numel(raw)
                 if(~isnumeric(raw{iraw}))
                     raw{iraw} = strrep(raw{iraw},'''',''); %neprisel jsem na zpusob, jak o udelat hromadne, isnumeric nefunguje na cely cellarray
                     %mozna by to slo po sloupcich, to ted neresim
                 end
             end
             PAC = cell2struct(raw(2:end,:),raw(1,:),2)';  %originalni PAC struktura z StructFind ma rozmer 1 x N, takze transponuju z excelu
             disp( [ basename(xlsfile) ': soubor nacten']);
        end
%         function PAC2Bipolar(PAC)
%             %kdyz nemam strukturu PAC bipolarne a potrebuju ji
%             PACpac = unique({PAC.pacient});
%             [ pacienti, setup ] = pacienti_setup_load( testname );            
%             for iPAC = 1:numel(PACpac)
%                 p = pacient_find(pacienti,PACpac{iPAC});
%                 if p < 0
%                     disp(['pacient nenalezen: ' nick]);                        
%                     continue;
%                 end
%                 disp(['* ' pacienti(p).folder ' - ' pacienti(p).header ' *']);                
%                 hfilename = [setup.basedir pacienti(p).folder '\' pacienti(p).header];                
%                 if exist(hfilename,'file')==2
%                     load(hfilename);
%                 else
%                     disp(['header ' hfilename ' neexistuje']);
%                     continue; %zkusim dalsiho pacienta, abych vypsal, ktere vsechny headery neexistujou
%                 end  
%                 
%             end
%         end
        function MIS = StructFindErr(testname)
            [ pacienti, setup ] = pacienti_setup_load( testname );
            load('BrainAtlas_zkratky.mat');
            MIS = {}; %pacient, ch, zkratka-  z toho bude vystupni xls tabulka s prehledem vysledku            
            iMIS = 1;
            for p = 1:numel(pacienti)
                disp(['* ' pacienti(p).folder ' - ' pacienti(p).header ' *']);
                hfilename = [setup.basedir pacienti(p).folder '\' pacienti(p).header];                
                if exist(hfilename,'file')==2
                    load(hfilename);
                else
                    disp(['header ' hfilename ' neexistuje']);
                    continue; %zkusim dalsiho pacienta, abych vypsal, ktere vsechny headery neexistujou
                end  
                for ch = 1:numel(H.channels)
                    z = strsplit(H.channels(ch).neurologyLabel,{'/','(',')'});
                    for iz = 1:numel(z)
                        if isempty(find(~cellfun('isempty',strfind(lower(BrainAtlas_zkratky(:,1)),lower(z{iz}))), 1)) %#ok<NODEF>
                           MIS(iMIS).pac = pacienti(p).folder; %#ok<AGROW>
                           MIS(iMIS).ch = ch; %#ok<AGROW>
                           MIS(iMIS).neurologyLabel = H.channels(ch).neurologyLabel; %#ok<AGROW>
                           MIS(iMIS).label = z{iz}; %#ok<AGROW>
                           MIS(iMIS).brainAtlas = H.channels(ch).ass_brainAtlas; %#ok<AGROW>
                           MIS(iMIS).cytoarchMap = H.channels(ch).ass_cytoarchMap; %#ok<AGROW>
                           iMIS = iMIS+ 1;
                        end
                    end
                end
            end
        end
    end
    methods (Access=private)
        function n = pocetcykluPlot3D(obj,kategorie,signum)
            %spocita kolik kanalu celkem vykresli PlotBrain3D pro tyto parametry
            n = 0; 
            for interval = 1:size(obj.VALS,1)  
                for kat = kategorie
                    if ~strcmp(obj.katstr{kat},'AllEl') %nechci to pro kategorii vsech elektrod
                        if signum > 0 
                            iV = obj.VALS{interval,kat} > 0; %jen kladne rozdily
                        elseif signum <0 
                            iV = obj.VALS{interval,kat} < 0; %jen zaporne rozdily
                        else
                            iV = true(size(obj.VALS{interval,kat})); %vsechny rozdily
                        end
                        n = n + numel(obj.VALS{interval,kat}(iV));
                    end
                end
            end
        end
    end
    
end

