classdef CBrainPlot < matlab.mixin.Copyable
    %CBRAINPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        VALS; %souhrnne prumery elektrod pres vsechny pacienty - cell(intervaly x kategorie)
        MNI;  %souhrnna MNI data pres vsechny pacienty - cell(intervaly x kategorie)
        NAMES; %souhrnna jmena elektrod pres vsechny pacienty - cell(intervaly x kategorie)
        NLabels; %jmena neurologyLabels z Headeru
        EpiInfo; %info o epilepsii z Headeru
        PVALS; %summary p values from statistics across all patients - cell(intervaly x kategorie)
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
        label; %label importovan z BPD dat
        signum; %kopie signum z IntervalyResp
    end
    properties (Constant)
        plevel = 0.05; % used significance level
    end
    methods (Access = public)        
        function [obj] = IntervalyResp(obj,testname,intervals,filename,contrast,signum,pvals_time)
            %IntervalyResp(testname,intervals,filename,contrast)
            %vola postupne pro vsechny pacienty E.IntervalyResp a uklada vysledky
            %vyradi vsechny kontakty bez odpovedi nebo se zapornou odpovedi
            %spoji vsechno dohromady
            %vrati vysledky ve formatu pro SEEE-vizualization
            %napr CB.IntervalyResp('aedist',[0.2 0.8],'AEdist CHilbert 50-120 refBipo Ep2017-11_CHilb.mat');
            %contrast - cislo statistiky, kterou chci pouzit
            %signum = jestli chci jen kat1>kat2 (1), nebo obracene (-1), nebo vsechny (0)
            
            if ~exist('contrast','var'), contrast = 1; end %defaultni je prvni kontrast            
            if ~exist('signum','var'), signum = 0; end %defaultne chci oba smery rozdilu mezi kategoriemi
            if ~exist('pvals_time','var'), pvals_time = 1; end %defaultne sbiram vsechny hodnoty ze vsech kanalu i casu 
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
            assert(size(intervals,2)==2,'musi byt prave dve hodnoty v kazdem intervalu');
            obj.intervals = intervals;             
            obj.filename = filename;
            elcount = []; %jen inicializace            
            P = {}; M = {}; N = {}; %jen inicializace
            obj.PAC = [];
            obj.pacients = cell(numel(pacienti),1); 
            obj.katstr_pacients = []; %musim to smazat, nize testuju, jestil to je prazdne
            obj.numelP = [];  %tam budu ukladat pocty elektrod pro kazdy pacient x interval x kategorie
            obj.signum = signum; %zavadim kvuli CMlabel
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
                    [prumery, MNI,names,~,katstr,neurologyLabels,pvalues] = E.IntervalyResp( intervals,[],signum,0,pvals_time);   %#ok<PROPLC> %no figure, funkce z CiEEGData                           
                    epiInfo = E.CH.GetChEpiInfo(); %uses all channels of E
                    obj.pacients{p} = pacienti(p).folder;
                    obj.GetPAC(prumery,E.CH.H,pacienti(p).folder); %gets and stores PAC structure for all intervals and categories to obj.PAC
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
                        EPI = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne neurologyLabels
                        PV = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne p values pro vsechny pacienty: pacient*interval*kategorie
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
                                NL{p,interval,kat}= strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],sum(ip),1)),neurologyLabels(ip)); 
                                EPI{p,interval,kat}= epiInfo(ip); 
                                PV{p,interval,kat}=squeeze(pvalues(interval, kat,:,:)); %ntervaly x kategorie x channels (x time) - pvalues from all channels
                                PVtime = size(pvalues,4); %n time samples , same for all patietns, to use later
                                elcount(interval,kat) = elcount(interval,kat) + sum(ip); %#ok<AGROW>
                                obj.numelP(p,interval,kat)=sum(ip);
                            else %kategorie jakoby navic pro vykresleni jen pozice elekrod
                                channels = size(prumery,1);
                                P{p,interval,kat}=zeros(channels,1); %#ok<AGROW> % 0 pro kazdy kanal - vsechny stejnou barvou
                                M{p,interval,kat}=MNI; %#ok<AGROW,PROPLC>>
                                N{p,interval,kat} = strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],channels,1)),names); %#ok<AGROW>
                                NL{p,interval,kat}= strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],channels,1)),neurologyLabels); 
                                EPI{p,interval,kat}= epiInfo; 
                                PV{p,interval,kat}=iff(pvals_time==1,ones(channels,PVtime),ones(channels,1)); %1 like nothing significant
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
            obj.EpiInfo = cell(size(elcount)); 
            obj.PVALS = cell(size(elcount)); %souhrnne p values - interval * kategorie
            obj.selCh = cell(size(elcount));  %prazdny vyber elektrod           
            if sum([pacienti.todo])>0 
                for interval = 1:size(prumery,2) 
                    for kat = 1:size(prumery,3)+1                   
                          obj.VALS{interval,kat} = zeros(elcount(interval,kat),1);
                          obj.MNI{interval,kat} = struct('MNI_x',{},'MNI_y',{},'MNI_z',{});
                          obj.NAMES{interval,kat}   = cell(elcount(interval,kat),1);
                          obj.NLabels{interval,kat} = cell(elcount(interval,kat),1);
                          obj.EpiInfo{interval,kat} = zeros(elcount(interval,kat),1);
                          obj.PVALS{interval,kat} = ones(elcount(interval,end),PVtime); %all channels in the last category x time (assuming the same for all patients
                          iVALS = 1; %index in VALS - only significant channels
                          iVALSall = 1; %index i PVALS - all channels
                          for p = 1:numel(pacienti) 
                              if pacienti(p).todo
                                  n = numel(P{p,interval,kat});
                                  n_all = size(PV{p,interval,kat},1); %PV contains all pvalues for all channels, channels x time
                                  obj.VALS{interval,kat} (iVALS:iVALS+n-1)=P{p,interval,kat};
                                  obj.MNI{interval,kat}  (iVALS:iVALS+n-1)=M{p,interval,kat};
                                  obj.NAMES{interval,kat}(iVALS:iVALS+n-1)  =N{p,interval,kat};
                                  obj.NLabels{interval,kat}(iVALS:iVALS+n-1)=NL{p,interval,kat};
                                  obj.EpiInfo{interval,kat}(iVALS:iVALS+n-1)=EPI{p,interval,kat};
                                  obj.PVALS{interval,kat} (iVALSall:iVALSall+n_all-1,:)=PV{p,interval,kat};
                                  iVALS = iVALS + n;
                                  iVALSall = iVALSall + n_all; 
                              end
                          end
                    end
                end             
                fprintf('\n'); %2 empty lines
                                
                %jeste vypisu pocty elektrod pro kazdou kategorii
                disp(['pacients processed: ' num2str(sum(~cellfun(@isempty,obj.pacients))) ]);
                fprintf('number of channels in %i categories (across all patients):\n',numel(obj.katstr));
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
        function [obj] = FdrVALS(obj,fdr,label)
            %performs fdr correction on obj.PVALS and corrects the channels in obj.VALS according to results
            if ~exist('fdr','var') || isempty(fdr), fdr = 1; end %default is less strict
            if ~exist('label','var'), label = ''; end %label of the output xls file
            if fdr == 2, fdrmethod = 'dep'; else fdrmethod = 'pdep'; end %#ok<SEPEX> %dep je striktnejsi nez pdep
            chremoved = zeros(size(obj.PVALS));
            %data for output xls table
            out = zeros(size(obj.PVALS{1,1},1),size(obj.PVALS,1) * numel(obj.iPAC) * 2 + 1); %channels x (kats * intervals * 2)
            out(:,1) = 1:size(obj.PVALS{1,1},1); %channel numbers
            outnames = cell(size(obj.PVALS{1,1},1),size(obj.PVALS,1) * numel(obj.iPAC)); %names of significant channels
            colnames = cell(1,size(out,2)+size(outnames,2)); %column names
            colnames{1} = 'ch';
            
            for interval = 1:size(obj.PVALS,1) %over time intervals
                for kat = 1:size(obj.PVALS,2) %over categories                          
                    ich0 = min(obj.PVALS{interval,kat},[],2)<obj.plevel; %index of originally significant channels
                    [~, ~, adj_p]=fdr_bh(obj.PVALS{interval,kat},obj.plevel,fdrmethod,'no'); 
                    if kat <= numel(obj.iPAC) %not for last category with all channels
                        coln = ((interval-1)*size(obj.PVALS,2) + (kat-1))*2+2;   %column index in out      2 columns per interval and kat + 1 first                                  
                        out(:,coln)=min(obj.PVALS{interval,kat},[],2); %minimum from each row, which is time
                        out(:,coln+1)=min(adj_p,[],2);  %minimum from each row, which is time                                              
                        coln = (interval-1)*size(obj.PVALS,2) + kat; %column index in outnames             1 column per interval and kat
                        outnames(ich0,coln) = obj.NAMES{interval,kat};
                        coln = ((interval-1)*size(obj.PVALS,2) + (kat-1))*3+2;  %column index in colnames - 3 columns per interval and kat + 1 first
                        colnames{coln} = ['int' num2str(interval) '_' obj.katstr{kat} '_name'];
                        colnames{coln+1} = ['int' num2str(interval) '_' obj.katstr{kat} '_p'];
                        colnames{coln+2} = ['int' num2str(interval) '_' obj.katstr{kat} '_padj'];                        
                    end
                    if sum(sum(adj_p<1))>0 %if some numbers smaller than one
                        obj.PVALS{interval,kat} = adj_p; %replace original pvalues with corrected ones - to correspond to obj.VALS etc                    
                        adj_p = adj_p(ich0,:); %leave only those originally significant
                        ich1 = min(adj_p,[],2)<obj.plevel; %index of now significant from the original ones                    
                        %filter out original values according to fdr corretion
                        obj.VALS{interval,kat} = obj.VALS{interval,kat}(ich1);
                        obj.MNI{interval,kat}  = obj.MNI{interval,kat}(ich1);
                        obj.NAMES{interval,kat} = obj.NAMES{interval,kat}(ich1);
                        obj.NLabels{interval,kat} = obj.NLabels{interval,kat}(ich1);
                        obj.EpiInfo{interval,kat} = obj.EpiInfo{interval,kat}(ich1);
                        obj.PAC{interval,kat} = obj.PAC{interval,kat}(ich1); %PAC table
                        obj.iPAC(interval,kat) = numel(obj.PAC{interval,kat});
                        chremoved(interval,kat) = sum(min(adj_p,[],2)>=obj.plevel); 
                    else
                        %for kat AllEl with all channels and all pvalues = 1
                        %do not filter out any channels
                        chremoved(interval,kat) = 0; 
                    end
                end
            end
            xlsfilename = ['logs\FdrVALS_' label '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.xls'];
            outcell = num2cell(out);
            for c = 1:size(outnames,2)
                ic = (c-1)*3+2;
                outcell = [outcell(:,1:(ic-1)),outnames(:,c),outcell(:,ic:end)]; %insert  column with channel names                
            end            
            xlswrite(xlsfilename,vertcat(colnames,outcell)); %zapisu do xls tabulky
            disp(['p values orig + fdr corrected written to ' xlsfilename]);
            fprintf('number of removed/left channels in %i categories (across all patients):\n',numel(obj.katstr));
            for kat = 1:numel(obj.katstr)
                fprintf('%s:\t', obj.katstr{kat});
                for interval = 1:size(obj.PVALS,1)
                    fprintf(' %i/%i,', chremoved(interval,kat), numel(obj.VALS{interval,kat}));
                end
            fprintf('\n');
            end  
        end
        function label = CMLabel(obj,intv,kat)
            %vrati standardni label pro CM.ExtractData, ve tvaru katstr_(intervalstring)_sigX, format podle BatchExtracts
            %zatim pouze po CB.IntervalyResp
            intvstr = sprintf('(%1.1f-%1.1f)',obj.intervals(intv,:)); %pojmenovani intervalu
            label = [ obj.katstr{intv,kat} '_' intvstr '_sig' num2str(obj.signum)];
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
            obj.EpiInfo = cell(size(BPD.EPI));
            for int = 1:size(BPD.EPI,1) %intervaly
                for kat = 1:size(BPD.EPI,2) %kategorie   
                    obj.EpiInfo{int,kat} = nan(numel(BPD.EPI{int,kat}),1); %nan zustanou u tech, ktere nemaji epiinfo v headeru
                    %po kanalech musim kvuli tomu, ze u nekterych kanalu neni ve struct ani 0 ani 1 ale [];
                    for ch = 1:numel(BPD.EPI{int,kat})
                        if ~isempty(BPD.EPI{int,kat}(ch).seizureOnset) % kanaly u kterych je epiinfo (nechybi v header)
                            obj.EpiInfo{int,kat}(ch) = double(BPD.EPI{int,kat}(ch).seizureOnset | BPD.EPI{int,kat}(ch).interictalOften); %vrati 1 pokud je jedno nebo druhe 1   
                        end
                    end                   
                end
            end
            obj.filename = BPD.filename;
            obj.label = BPD.label; %label z CHilbertMulti, napr PPA_PHGent_50-150Hz
            disp(['data importovana ze souboru "' BPD.filename '"']);
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
            if isstruct(cfg) && isfield(cfg,'NoNames')
                obj.plotBrain3Dcfg.NoNames = cfg.NoNames;
            else
                obj.plotBrain3Dcfg.NoNames = 0; %defaultne se nedelaji obrazky bez popisu elektrod            
            end            
            if isstruct(cfg) && isfield(cfg,'Names')
                obj.plotBrain3Dcfg.Names = cfg.Names;
            else
                obj.plotBrain3Dcfg.Names = 1; %defaultne se delaji obrazky s popisem elektrod            
            end   
            if isstruct(cfg) && isfield(cfg,'colorScale')
                obj.plotBrain3Dcfg.colorScale = cfg.colorScale;            
            end   
           %barevna skala od Nadi
           %default colors
            redscale = [162 2 2; 246 203 203 ];  %dark red %light red
            greenscale = [0 153 0;212 255 171];  %dark green %light green
            bluescale = [0 0 153; 0 153 255]; %dark blue light blue
            obj.plotBrain3Dcfg.customColors.customColor = true; % custom colormap oddeli negativne a pozitivne hodnoty - 29.6.2018
            obj.plotBrain3Dcfg.customColors.flip = 0; %pokud chci prehodit barvy
            obj.plotBrain3Dcfg.customColors.darkneg = greenscale(1,:); %dark green
            obj.plotBrain3Dcfg.customColors.lightneg = greenscale(2,:); %light green  
            
            %obj.plotBrain3Dcfg.customColors.lightpos = redscale(1,:); %dark red
            %obj.plotBrain3Dcfg.customColors.darkpos = redscale(2,:); %light red  
            
            obj.plotBrain3Dcfg.customColors.darkpos = bluescale(1,:); %dark blue
            obj.plotBrain3Dcfg.customColors.lightpos = bluescale(2,:); %light red 
            obj.plotBrain3Dcfg.customColors.zeroclr = [0 0 0]; %the color of zero values 
            obj.plotBrain3Dcfg.customColors.supermaxcolor = [192 192 192]; %the color of custom value
            
            %we can optionally overwrite the default values

            if isstruct(cfg) && isfield(cfg,'Colors')
                if isfield(cfg.Colors,'pos')
                    switch cfg.Colors.pos
                        case 'r'
                            obj.plotBrain3Dcfg.customColors.darkpos = redscale(1,:); %dark red
                            obj.plotBrain3Dcfg.customColors.lightpos = redscale(2,:); %light red 
                        case 'g'
                            obj.plotBrain3Dcfg.customColors.darkpos = greenscale(1,:); %dark green
                            obj.plotBrain3Dcfg.customColors.lightpos = greenscale(2,:); %light green 
                        case 'b'
                            obj.plotBrain3Dcfg.customColors.darkpos = bluescale(1,:); %dark blue
                            obj.plotBrain3Dcfg.customColors.lightpos = bluescale(2,:); %light red 
                    end
                end
                if isfield(cfg.Colors,'neg')
                    switch cfg.Colors.neg
                        case 'r'
                            obj.plotBrain3Dcfg.customColors.darkneg = redscale(1,:); %dark blue
                            obj.plotBrain3Dcfg.customColors.lightneg = redscale(2,:); %light red 
                        case 'g'
                            obj.plotBrain3Dcfg.customColors.darkneg = greenscale(1,:); %dark blue
                            obj.plotBrain3Dcfg.customColors.lightneg = greenscale(2,:); %light red 
                        case 'b'
                            obj.plotBrain3Dcfg.customColors.darkneg = bluescale(1,:); %dark blue
                            obj.plotBrain3Dcfg.customColors.lightneg = bluescale(2,:); %light red 
                    end
                end 
                if isfield(cfg.Colors,'zero')
                    obj.plotBrain3Dcfg.customColors.zeroclr = cfg.Colors.zero; %the color of zero values 
                end 
            end
            
        end
        function PlotBrain3D(obj,kategorie)
            %vykresli jpg obrazky jednotlivych kategorii a kontrastu mezi nimi            
            %TODO do jmena vystupniho jpg pridat i frekvence a referenci, aby se to neprepisovalo
            %TODO je mozne ty signif vyexportovat a pak je nacist zase do CHilbertMulti?
            %TODO do vystupni tabulky nejak dostat anatomickou lokalizaci?
            assert(~isempty(obj.VALS),'zadna data VALS');
            assert(~isempty(obj.plotBrain3Dcfg),'zadna konfigurace - je nutne volat PlotBrain3DConfig');
            plotSetup = {};
            if ~exist('kategorie','var') || isempty(kategorie) , kategorie = 1:size(obj.VALS,2); end %muzu chtit jen nektere kategorie
            signum = obj.plotBrain3Dcfg.signum; %#ok<PROPLC>
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
            plotSetup.myColorMap = iff(signum ~= 0,parula(128) ,jet(128));  %#ok<PROPLC>  %pokud jednostrane rozdily, chci parula
            
            plotSetup.customColors = obj.plotBrain3Dcfg.customColors; %barevna skala od Nadi       
            plotSetup.colorScale = obj.plotBrain3Dcfg.colorScale; %custom range to be plotted
            plotSetup.figureNamePrefix = [ obj.testname '_' num2str(obj.Hf([1 end]),'%i-%iHz') '_' obj.reference '_names']; %default name
            if strcmp(plotSetup.figureVisible,'off')
                disp('figures invisible');
            end
            if ~obj.plotBrain3Dcfg.NoNames 
               disp('obrazky NoNames negeneruju');
            end
            if ~obj.plotBrain3Dcfg.Names 
               disp('obrazky Names negeneruju');
            end
            [tablerows,VALS_rows,iV_All] = obj.pocetcykluPlot3D(kategorie,signum); %#ok<PROPLC>
            tablelog = cell(tablerows+2,7);% z toho bude vystupni xls tabulka s prehledem vysledku
            tablelog(1,:) = {datestr(now),obj.filename,'','','','',''}; %hlavicky xls tabulky
            tablelog(2,:) = {'interval','kategorie','chname','neurologyLabel','mni','val','selected'}; %hlavicky xls tabulky
            [ChMap,ChNames] = obj.ChannelMap();
            
            iTL = 2; %index v tablelog
            tic; %zadnu merit cas            
            for interval = 1:size(obj.VALS,1) 
                for ikat = 1:numel(kategorie) %cycle over categories plotted to different jpgs
                    VSum = VALS_rows(interval,ikat);
                    if VSum >0 %if anny channels to plot across all categories
                        
                        %1. fill the data for one JPG pictures
                        VALS_channels = zeros(VSum,1); %will containt all vals_channels together
                        MNI_channels = struct('MNI_x',{},'MNI_y',{},'MNI_z',{});
                        MNI_channels(VSum,1).MNI_x = 0; %inicializa the whole struct
                        
                        NAMES_channels = cell(VSum,1);
                        iVALS = 1;
                        katikat = cellval(kategorie,ikat); %index of internal category - plotted to the same figure
                        if isnumeric(katikat), katikat = num2cell(katikat); end                        
                        for ikat2=1:numel(katikat) %cycle over categories plotted to the same jpg
                            kat = cellval(katikat,ikat2);                              
                            ValsNegative = false;                            
                            plotSetup.customColors.SuperMax = false; 
                            if isstring(kat) || ischar(kat) %if kat is string it should be plotted with all same values in a specific color : customColors.supermaxcolor
                                kat = str2double(kat); % the string kat has to be the last one in kategorie
                                katikat{ikat2} = kat;
                                plotSetup.customColors.SuperMax = true; %transfer info to use the supermax value to main_brainPlot
                            elseif kat < 0
                                kat = abs(kat);
                                if iscell(katikat) %for supermax kat, the kats numbers hav to be in cellarray
                                    katikat{ikat2} = kat;
                                else
                                    katikat(ikat2) = kat;
                                end
                                ValsNegative = true;                            
                            end
                            iV = iV_All{interval,ikat}{ikat2}; 
                            if ikat2==1 %this will be set from the first category
                                katname = obj.katstr{kat};                                
                                plotSetup.circle_size = iff(strcmp(katname,'all') || strcmp(katname,'AllEl'),28,56); %mensi kulicka pokud vsechny elektrody                
                                if numel(katikat)>1
                                    katname = cell2str(obj.katstr(abs(cell2double(katikat))),1); %kats can be negative, and could be also string in cell array if SuperMax is used
                                end
                                brainlabel = obj.GetBrainLabel(); %pokud v label na druhe pozici je nazev mozkove oblasti 
                                figureNameNames = [ obj.testname brainlabel '_' num2str(obj.intervals(interval,:),'%.1f-%.1fs')  '_' katname '_' num2str(signum) ... %#ok<PROPLC>
                                        '_' num2str(obj.Hf([1 end]),'%i-%iHz') '_' obj.reference '_names'];
                                figureNameNoNames = [ obj.testname brainlabel '_' num2str(obj.intervals(interval,:),'%.1f-%.1fs')  '_' katname '_' num2str(signum) ... %#ok<PROPLC>
                                        '_' num2str(obj.Hf([1 end]),'%i-%iHz') '_' obj.reference '_NOnames'];
                            end
                            if numel(obj.VALS{interval,kat}(iV)) > 0

                                vals_channels = obj.VALS{interval,kat}(iV); %parametr  main_brainPlot
                                if signum ~= 0 %#ok<PROPLC>
                                    vals_channels = vals_channels*signum; %#ok<PROPLC> %u zapornych hodnot prehodim znamenko
                                end
                                if ValsNegative % this category should be plotted as negative;
                                    vals_channels = vals_channels*-1; 
                                elseif plotSetup.customColors.SuperMax
                                    vals_channels = ones(size(vals_channels))*max(max(vals_channels),max(VALS_channels))*1.01; %1 percent larger value, % the string kat has to be the last one in kategorie
                                    plotSetup.colorScale = [min(VALS_channels) max(VALS_channels)]; %range of values without this supermax value - all values until now
                                end 
                                mni_channels = obj.MNI{interval,kat}(iV);                                                                                                 
                                names_channels = iff(obj.plotBrain3Dcfg.NLabels, obj.NLabels{interval,kat}(iV), obj.NAMES{interval,kat}(iV));                        

                                if ~strcmp(obj.katstr{kat},'AllEl') %do not export all-channels category to xls
                                    for iVal = 1:numel(vals_channels)
                                        tablelog(iVal + iTL,:) = { sprintf('[%.1f %.1f]',obj.intervals(interval,:)),obj.katstr{kat}, obj.NAMES{interval,kat}{iVal}, obj.NLabels{interval,kat}{iVal}, ...
                                            sprintf('[%.1f,%.1f,%.1f]',mni_channels(iVal).MNI_x, mni_channels(iVal).MNI_y, mni_channels(iVal).MNI_z), vals_channels(iVal),int8(iV(iVal))};
                                        iChNames = contains(ChNames(:,1),obj.NAMES{interval,kat}{iVal});
                                        ChMap(iChNames,kat,interval) = vals_channels(iVal);
                                    end
                                    iTL = iTL + numel(vals_channels);
                                end
                                VALS_channels(iVALS:iVALS+sum(iV)-1,1) = vals_channels;
                                [MNI_channels(iVALS:iVALS+sum(iV)-1,1).MNI_x] = mni_channels(:).MNI_x;
                                [MNI_channels(iVALS:iVALS+sum(iV)-1,1).MNI_y] = mni_channels(:).MNI_y;
                                [MNI_channels(iVALS:iVALS+sum(iV)-1,1).MNI_z] = mni_channels(:).MNI_z;
                                NAMES_channels(iVALS:iVALS+sum(iV)-1) = names_channels;
                                iVALS = iVALS + sum(iV);                                                              
                            end
                        end
                        
                        %2. plot the JPGs
                            %without channel description
                            if ~isfield(plotSetup,'colorScale') || isempty(plotSetup.colorScale) %force range of values to main_brainPlot
                                plotSetup.colorScale = [min(VALS_channels) max(VALS_channels)];
                            elseif ~plotSetup.customColors.SuperMax
                                VALS_channels(VALS_channels < plotSetup.colorScale(1)) = plotSetup.colorScale(1);
                                VALS_channels(VALS_channels > plotSetup.colorScale(2)) = plotSetup.colorScale(2);                                    
                            end
                            if obj.plotBrain3Dcfg.NoNames 
                                if  isempty(dir([ plotSetup.outputDir '3D_model\' figureNameNoNames '*'])) || obj.plotBrain3Dcfg.overwrite==1 
                                    plotSetup.figureNamePrefix = figureNameNoNames;
                                    disp(plotSetup.figureNamePrefix);                            
                                    brainsurface = main_brainPlot(VALS_channels,MNI_channels,[],brainsurface,plotSetup);  %#ok<PROPLC>
                                    %volam Jirkuv skript, vsechny ty promenne predtim jsou do nej
                                    if isempty(obj.brainsurface)
                                        obj.brainsurface = brainsurface; %#ok<PROPLC> %ulozim si ho pro dalsi volani
                                    end
                                else
                                    disp(['the file already exists ' figureNameNoNames ' - not ovetwriting ']);
                                end                       
                            end
                            %with channel names
                            if obj.plotBrain3Dcfg.Names                                                  
                                if isempty(dir([ plotSetup.outputDir '3D_model\' figureNameNames '*'])) || obj.plotBrain3Dcfg.overwrite==1 
                                    plotSetup.figureNamePrefix = figureNameNames;
                                    disp(plotSetup.figureNamePrefix);                                                     
                                    brainsurface = main_brainPlot(VALS_channels,MNI_channels,NAMES_channels,brainsurface,plotSetup);    %#ok<PROPLC>  
                                    if isempty(obj.brainsurface)
                                        obj.brainsurface = brainsurface; %#ok<PROPLC> %ulozim si ho pro dalsi volani
                                    end
                                else
                                    disp(['the file already exists ' figureNameNames ' - not ovetwriting ']);
                                end                        
                            end
                        
                    else  
                            disp(['no values for ' plotSetup.figureNamePrefix ' - not saving jpg ']);  
                    end
                    hotovoProc = floor(((interval-1)*(numel(kategorie)) + kat) / numel(obj.VALS) * 100) ;                    
                    disp(['hotovo '  num2str(hotovoProc) '% za ' num2str(floor(toc/6)/10) ' minut']);
                end
            end
                      
            logfilename = [num2str([obj.intervals(1,1) obj.intervals(end,2)],'%.1f-%.1fs') '_sig' num2str(signum) '_' basename(obj.filename) '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') ];            
            xlsfilename = ['logs\PlotBrain3D_' logfilename '.xls'];
            xlswrite(xlsfilename,tablelog); %zapisu do xls tabulky
            disp([ 'xls tables saved: ' xlsfilename]);
            xlsfilename = ['logs\PlotBrain3D_MAP_' logfilename '.xls'];
            obj.ChannelMap2Xls(ChNames,ChMap,xlsfilename);
            disp([ 'xls tables saved: ' xlsfilename ]);
            
            toc; %ukoncim mereni casu a vypisu 
            if hybernovat
                system('shutdown -h')  %#ok<UNRCH>
            elseif vypnout            
                system('shutdown -s') %#ok<UNRCH>
            end
        end
        function [ChMap,ChNames] = ChannelMap(obj)
            %chci ziskat prehled ze vsech kanalu, ktere vykazuji nejake odpovedi na jednotlive kategorie
            if strcmp(obj.katstr(end),'AllEl')
                ChNames = obj.NAMES{1,end}; %seznam vsech kanalu pres vsechny elektrody                
                ChLabels = obj.NLabels{1,end}; %seznam vsech Neurol lokalizaci pres vsechny elektrody  
                ChMNI = obj.MNI{1,end};
                ChEpiInfo = obj.EpiInfo{1,end}; %potrebuju double, protoze cast za else mi taky poskytuje double
            else
                ChNames = {}; ChLabels = {}; ChMNI = struct('MNI_x',{},'MNI_y',{},'MNI_z',{}); ChEpiInfo = [];
                for k = 1:numel(obj.NAMES) %cyklus pres vsechny intervaly a kategorie
                    ChNames = cat(1,ChNames,obj.NAMES{k}); % vynecha duplikaty, seradi
                    ChLabels = cat(1,ChLabels,obj.NLabels{k});
                    ChMNI = cat(1,ChMNI, obj.MNI{k}'); %#ok<AGROW>
                    ChEpiInfo = cat(1,ChEpiInfo,obj.EpiInfo{k}); %epiinfo je double
                end
                [ChNames,iChNames] = unique(ChNames); % chci polozky seradit
                ChLabels = ChLabels(iChNames);
                ChMNI = ChMNI(iChNames);
                ChEpiInfo = ChEpiInfo(iChNames);
            end
            ChEpiInfo = num2cell(ChEpiInfo); %kvuli exportu do excelu potrebuju cell array
            ChEpiInfo(isnan(cell2mat(ChEpiInfo))) = {'NaN'}; %excel neumi zapsat nan hodnoty, musi to byt string
            ChNames = cat(2,ChNames,ChLabels,ChEpiInfo,{ChMNI.MNI_x}',{ChMNI.MNI_y}',{ChMNI.MNI_z}'); %budu mit v jednom cellarray vic sloupcu
            ChMap = zeros(numel(ChNames),numel(obj.katstr),size(obj.intervals,1)); %tam budu  ukladat odpovedi pro jednotlive kanaly, kategorie a intervaly          
        end  
        function ChannelMap2Xls(obj,ChNames,ChMap,xlsfilename)
            %zapise ChMap do xls tabulky spolu se jmeny a nlabely kanalu
            %pro kazdou kategorii bere maximalni absolutni hodnotu pres vsechny intervaly
            tablexls = cell(size(ChNames,1)+1,size(ChMap,2)+size(ChNames,2)+1); %sloupce navic - ChNames (+ ChLabels+MINxyz) + AnyKat
            tablexls(1,:) = cat(2,{'channel','nlabel','epileptic','MNI_x','MNI_y','MNI_z','anykat'},obj.katstr);
            for iCh = 1:size(ChNames,1)
                [~,im]=max(abs(ChMap(iCh,:,:)),[],3); %indexy maximalnich absolutnich hodnot pres cas
                M = ChMap(sub2ind(size(ChMap),iCh*ones(1,numel(im)),1:numel(im),im)); %vyzvednu ty maximalni hodnoty pomoci absolutniho indexovani - pro kazdou polozku im jeden index
                tablexls(iCh+1,:)=cat(2,ChNames(iCh,:),num2cell(int8(any(M~=0))),num2cell(M)); %ukladam maximalni hodnoty ze vsech casu pro kazdou kategorii zvlast
            end
            xlswrite(xlsfilename ,tablexls); %zapisu do xls tabulky
        end
    end
    
    methods (Static,Access = public)
        function [PAC,ChNum] = StructFind(struktura,label,testname,reference,labelnot,pactodo)
            %najde pacienty, jejich headery obsahuji mozkovou strukturu
            %struktura je nazev struktury podle atlas napriklad hippo, label je kratky nazev podle martina, napriklad hi
            if ~exist('label','var'),    label = struktura; end %defaultni test
            if ~exist('testname','var') || isempty(testname), testname = 'aedist'; end %defaultni test
            if ~exist('reference','var') || isempty(reference), reference = []; end %defaultni test
            if ~exist('labelnot','var'),    labelnot = {}; end %defaultni test
            if ~exist('pactodo','var'), pactodo = 0; end %jestli se maji brat jen pacienti s todo=1
            if ischar(struktura), struktura = {struktura}; end %prevedu na cell array
            if ischar(label), label = {label}; end %prevedu na cell array
            [ pacienti, setup ] = pacienti_setup_load( testname );
            PAC = {};
            iPAC = 1;
            ChNum = nan(numel(pacienti),2); %number of channels for each pacient, columns pacient no, num of channels             
            for p = 1:numel(pacienti)
                index = [];
                if pacienti(p).todo == 1 || ~pactodo %pokud je u pacienta todo, nebo se nema pouzivat
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
                ChNum(p,:) = [p,numel(index)];
            end 
            if numel(PAC) > 0
                xlsname = ['./logs/StructFind PAC_' testname '_' cell2str(struktura,1) '_' cell2str(label,1) '.xlsx'];
                writetable(struct2table(PAC), xlsname); %ulozimvysledek taky do xls
                disp(['found ' num2str(nansum(ChNum(:,2))) ' channels in ' num2str(sum(~isnan(ChNum(:,1)))) ' pacients']);
            else
                disp('no channels found');
            end
        end
        function PAC = MNIFind(XYZ,distance,testname,reference,pactodo,distNames)
            %finds channels across all pacients with a maximal distance to given XYZ MNI coordinates
            %XYZ and distance can contain multiple (but equal) number of values-positions
            %reference is the bipolar or other reference
            if ~exist('testname','var') || isempty(testname), testname = 'aedist'; end %defaultni test            
            if ~exist('reference','var') || isempty(reference), reference = []; end %default is original reference 
            if ~exist('pactodo','var') || isempty(pactodo), pactodo = 0; end %jestli se maji brat jen pacienti s todo=1
            if ~exist('distNames','var'), distNames = num2cell(1:numel(distance)); end %optional names for XYZ-distance sets
            
            assert(size(XYZ,1)==numel(distance),'the number of xyz coordinates and distances must to be the same');
            [ pacienti, setup ] = pacienti_setup_load( testname );
            PAC = {};
            iPAC = 1;
            for p = 1:numel(pacienti)
                if pacienti(p).todo == 1 || ~pactodo %pokud je u pacienta todo, nebo se nema pouzivat
                    disp(['* ' pacienti(p).folder ' - ' pacienti(p).header ' *']);
                    hfilename = [setup.basedir pacienti(p).folder '\' pacienti(p).header];                
                    if exist(hfilename,'file')==2
                        load(hfilename); %#ok<LOAD> %load the header
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

                    % MAIN SEARCH OF CHANNES FOR THIS PACIENT 
                    MNI = [H.channels.MNI_x; H.channels.MNI_y; H.channels.MNI_z]';                    
                    [D1,I1] = pdist2(XYZ,MNI,'euclidean','Smallest',1); %vector of smallest distance of each MNI to any of XYZ
                    XYZ2 = [-XYZ(:,1) XYZ(:,2)  XYZ(:,3) ]; %coordinates on the left brain side 
                    [D2,I2] = pdist2(XYZ2,MNI,'euclidean','Smallest',1); %vector of distance 

                    iCh = false(1,size(MNI,1)); %index nalezenych kanalu v ramci tohoto pacienta
                    iChannelDist = cell(1,size(MNI,1));
                    for iDist = 1:numel(distance)
                       iIn = find(I1==iDist); %numbers of channels closest to this XYZ-distance set
                       iCh0 = D1(iIn)<distance(iDist);                       
                       if(sum(iCh0)>0) %some channels found
                         iCh(iIn)=iCh(iIn) | iCh0;                        
                         iChannelDist(iIn(iCh0))=repmat(distNames(iDist),1,numel(iIn(iCh0)));
                       end
                       iIn = find(I2==iDist); %numbers of channels closest to this XYZ-distance set
                       iCh0 = D2(iIn)<distance(iDist);
                       if(sum(iCh0)>0) %some channels found
                         iCh(iIn)=iCh(iIn) | iCh0; 
                         iChannelDist(iIn(iCh0))=repmat(distNames(iDist),1,numel(iIn(iCh0)));
                       end
                    end
                    iCh = find(iCh); 
                    %pokud jsem kanaly nevyradil uz pri zmene reference - vyrazuji se jen pri bipolarni
                    if isempty(reference) || reference ~= 'b' 
                        indexvyradit = ismember(iCh, pacienti(p).rjch); %vyrazene kanaly tady nechci
                        iCh(indexvyradit)=[]; 
                    end

                    %vrati indexy radku ze struct array, ktere obsahuji v sloupci neurologyLabel substring struktura
                    for ii = 1:numel(iCh)                
                        PAC(iPAC).pacient = pacienti(p).folder; %#ok<AGROW>
                        PAC(iPAC).ch = iCh(ii); %#ok<AGROW>
                        PAC(iPAC).name = H.channels(iCh(ii)).name; %#ok<AGROW>
                        PAC(iPAC).neurologyLabel = H.channels(iCh(ii)).neurologyLabel; %#ok<AGROW>
                        PAC(iPAC).ass_brainAtlas = H.channels(iCh(ii)).ass_brainAtlas;%#ok<AGROW>
                        PAC(iPAC).ass_cytoarchMap = H.channels(iCh(ii)).ass_cytoarchMap; %#ok<AGROW>
                        PAC(iPAC).MNI_x = H.channels(iCh(ii)).MNI_x; %#ok<AGROW>
                        PAC(iPAC).MNI_y = H.channels(iCh(ii)).MNI_y; %#ok<AGROW>
                        PAC(iPAC).MNI_z = H.channels(iCh(ii)).MNI_z; %#ok<AGROW>
                        PAC(iPAC).MNIdist = min(D1(iCh(ii)),D2(iCh(ii))); %#ok<AGROW>
                        PAC(iPAC).DistN = iChannelDist(iCh(ii)); %#ok<AGROW>
                        iPAC = iPAC + 1;
                    end
                end
            end
            if numel(PAC) > 0
                if size(XYZ,1) == 1                    
                    xlsname = ['./logs/MNIFind PAC_' testname '_mni' num2str(XYZ,'(%3.1f %3.1f %3.1f)') '_dist' num2str(distance) '.xlsx'];
                else
                    d = num2str(round(distance),'%x'); %some hexadecimal code for the distance array
                    d(d==' ')=[]; %strip all spaces
                    xlsname = ['./logs/MNIFind PAC_' testname '_mni' num2str(size(XYZ,1)) 'x_dist' d '.xlsx'];
                end
                writetable(struct2table(PAC), xlsname); %ulozimvysledek taky do xls
                disp(['table with ' num2str(numel(PAC)) ' channels written to ' xlsname]);
            else
                disp('no channels found');
            end
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
        function PAC = MergePAC(PAC1,PAC2,varargin)
            %funkce spoji nekolik PAC struktur dohromady, variabilni seznam argumentu >= 2
            PAC = [PAC1,PAC2];
            for n = 1:numel(varargin)
                PAC = [PAC,varargin{n}]; %#ok<AGROW>
            end
            T = struct2table(PAC); %prevedu na table, aby fungovalo uniue            
            T = unique(T,'rows'); %cele unikatni radky
            PAC = table2struct(T)'; %chci mit prvni rozmer 1 a radky ve druhem rozmeru, jako jsou vsechny PAC
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
                        if isempty(find(~cellfun('isempty',strfind(lower(BrainAtlas_zkratky(:,2)),lower(z{iz}))), 1)) %#ok<NODEF>
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
        function [tablerows,VALS_rows,iV_All] = pocetcykluPlot3D(obj,kategorie,signum)
            %counts how many channels to plot for each category and interval
            %also creates logical indexing of channels for each category and interval
            %spocita kolik kanalu celkem vykresli PlotBrain3D pro tyto parametry
            tablerows = 0; %number of rows in the output xls table
            VALS_rows = zeros(size(obj.VALS,1),numel(kategorie)); %number of channels (rows) to be plotted
            iV_All = cell(size(obj.VALS,1),numel(kategorie)); %intervals x kategories %to store indexes of channels for each category 
            for interval = 1:size(obj.VALS,1)  %inervals to be plotted
                for ikat = 1:numel(kategorie) %kategories to be plotted 
                    katikat = cellval(kategorie,ikat);
                    iV = cell(numel(katikat),1); %logical index of channels to be plotted in this internal category
                    for ikat2=1:numel(katikat) %INTERNAL CATEGORY: we can have multiple categories in one picture. In this order                                                          
                        kat = cellval(katikat,ikat2);
                        if isstring(kat) || ischar(kat), kat = str2double(kat); end %because SuperMax values, which are as chars
                        kat = abs(kat); %kategory number can be negative
                        if signum > 0 
                            iV{ikat2} = obj.VALS{interval,kat} > 0; %jen kladne rozdily
                        elseif signum <0 
                            iV{ikat2} = obj.VALS{interval,kat} < 0; %jen zaporne rozdily
                        elseif ~isempty(obj.selCh{interval,kat})
                            iV{ikat2} = ismember(1:numel(obj.VALS{interval,kat}),find(any(obj.selCh{interval,kat},2))); %vyber kanalu k zobrazeni, napriklad z CHilbertMulti
                        else
                            iV{ikat2} = true(size(obj.VALS{interval,kat})); %vsechny rozdily
                        end
                        if ~strcmp(obj.katstr{kat},'AllEl') %nechci xls rows pro kategorii vsech elektrod
                            tablerows = tablerows +  sum(iV{ikat2});
                        end
                        VALS_rows(interval,ikat) = VALS_rows(interval,ikat) + sum(iV{ikat2});
                    end
                    iV_All{interval,ikat} = iV;
                end
            end
        end  
        function bl = GetBrainLabel(obj)
            bl = '';
            if ~isempty(obj.label)
                [~,interval,~] = CHilbertMulti.GetLabelInfo(obj.label);
                if ischar(interval)
                    bl = ['_' interval];  %pokud druha polozka labelu z CM je string, nejspis je to zkratka mozkove oblasti, napriklad PHGent              
                    %podtrzitko kvuli pozdejsimu vlozeni do nazvu souboru 
                end            
            end
        end
    end
    
end

