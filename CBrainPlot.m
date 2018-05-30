classdef CBrainPlot < handle
    %CBRAINPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        VALS; %souhrnne prumery elektrod pres vsechny pacienty - cell(intervaly x kategorie)
        MNI;  %souhrnna MNI data pres vsechny pacienty - cell(intervaly x kategorie)
        NAMES; %souhrnna jmena elektrod pres vsechny pacienty - cell(intervaly x kategorie)
        intervals; % intervaly z funkce IntervalyResp
        katstr; %jmena kategorii
        brainsurface; %ulozeny isosurface z main_brainPlot
        testname; %jmeno zpracovavaneho testu
        katstr_pacients;
        numelP; %pocty  signif elektrod pro kazdy pacient x interval x kategorie
        pacients;         
    end
    
    methods (Access = public)        
        function [obj] = IntervalyResp(obj,testname,intervals,filename,contrast)
            %IntervalyResp(testname,intervals,filename,contrast)
            %vola postupne pro vsechny pacienty E.IntervalyResp a uklada vysledky
            %vyradi vsechny kontakty bez odpovedi nebo se zapornou odpovedi
            %spoji vsechno dohromady
            %vrati vysledky ve formatu pro SEEE-vizualization
            %napr CB.IntervalyResp('aedist',[0.2 0.8],'AEdist CHilbert 50-120 refBipo Ep2017-11_CHilb.mat');
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
            elcount = []; %jen inicializace            
            P = {}; M = {}; N = {}; %jen inicializace
            obj.pacients = cell(numel(pacienti),1); 
            obj.katstr_pacients = []; %musim to smazat, nize testuju, jestil to je prazdne
            obj.numelP = [];  %tam budu ukladat pocty elektrod pro kazdy pacient x interval x kategorie
            for p = 1:numel(pacienti) % cyklus pacienti
                if pacienti(p).todo 
                    disp(['***   ' pacienti(p).folder '   ***']);
                    E = pacient_load(pacienti(p).folder,testname,filename); %nejspis objekt CHilbert, pripadne i jiny
                    if isempty(E)
                        disp('no data');
                        pacienti(p).todo = 0; %nechci ho dal zpracovavat
                        continue;
                    end
                    E.SetStatActive(contrast); %nastavi jeden z ulozenych statistickych kontrastu
                    [prumery, MNI,names,~,katstr] = E.IntervalyResp( intervals,[],0);   %#ok<PROPLC> %no figure, funkce z CiEEGData                           
                    obj.pacients{p} = pacienti(p).folder;
                    clear E;
                    if isempty(obj.katstr_pacients)
                        obj.katstr = [katstr 'AllEl']; %#ok<PROPLC> 
                        obj.katstr_pacients = cell(numel(pacienti),numel(katstr)); %#ok<PROPLC>
                        obj.numelP = zeros(numel(pacienti),size(intervals,1),numel(katstr)+1); %#ok<PROPLC> %tam budu ukladat pocty elektrod pro kazdy interval a pacienta
                        elcount = zeros(size(prumery,2),size(prumery,3)+1); %pocet elektrod pro kazdy casovy interval a kategorii - interval x kategorie
                        P = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne prumery pro vsechny pacienty: pacient*interval*kategorie
                        M = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne MNI koordinaty pro vsechny pacienty
                        N = cell([numel(pacienti),size(prumery,2),size(prumery,3)+1]); % souhrnne names pro vsechny pacienty
                            %+1 je pro obrazek vsech elektrod i tech bez odpovedi
                    end
                    obj.katstr_pacients(p,:) = katstr; %#ok<PROPLC> %jsou kategorie u vsech pacientu ve stejnem poradi?
                    for interval = 1:size(prumery,2) % cyklus intervaly
                        for kat = 1:size(prumery,3)+1 % cyklus kategorie podnetu
                            if kat <= size(prumery,3) %obvykle kategorie
                                ip = prumery(:,interval, kat) ~= 0; % chci i zaporny rozdil ; aby tam neco bylo 
                                P{p,interval,kat}=prumery(ip,interval, kat); %#ok<AGROW>
                                M{p,interval,kat}=MNI(ip); %#ok<AGROW,PROPLC>>
                                N{p,interval,kat}= strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],sum(ip),1)),names(ip)); %#ok<AGROW>
                                elcount(interval,kat) = elcount(interval,kat) + sum(ip); %#ok<AGROW>
                                obj.numelP(p,interval,kat)=sum(ip);
                            else %kategorie jakoby navic pro vykresleni jen pozice elekrod
                                channels = size(prumery,1);
                                P{p,interval,kat}=zeros(channels,1); %#ok<AGROW> % 0 pro kazdy kanal - vsechny stejnou barvou
                                M{p,interval,kat}=MNI; %#ok<AGROW,PROPLC>>
                                N{p,interval,kat}= strcat(cellstr(repmat([pacienti(p).folder(1:4) '_'],channels,1)),names); %#ok<AGROW>
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
            if sum([pacienti.todo])>0 
                for interval = 1:size(prumery,2) 
                    for kat = 1:size(prumery,3)+1                   
                          obj.VALS{interval,kat} = zeros(elcount(interval,kat),1);
                          obj.MNI{interval,kat} = struct('MNI_x',{},'MNI_y',{},'MNI_z',{});
                          obj.NAMES{interval,kat} = cell(elcount(interval,kat),1);
                          iVALS = 1;
                          for p = 1:numel(pacienti) 
                              if pacienti(p).todo
                                  n = numel(P{p,interval,kat});
                                  obj.VALS{interval,kat} (iVALS:iVALS+n-1)=P{p,interval,kat};
                                  obj.MNI{interval,kat}  (iVALS:iVALS+n-1)=M{p,interval,kat};
                                  obj.NAMES{interval,kat}(iVALS:iVALS+n-1)=N{p,interval,kat};
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
        function obj = ImportData(obj,BPD)
            %vlozi data, ktera jsem vytvoril pomoci CHilbert.ExtractBrainPlotData
            obj.VALS = BPD.VALS;
            obj.MNI = BPD.MNI;
            obj.NAMES = BPD.NAMES;
            obj.katstr = BPD.katstr;
            obj.intervals = BPD.intervals;       
            obj.testname = BPD.testname; 
        end
        function PlotBrain3D(obj,kategorie,signum,outputDir)
            assert(~isempty(obj.VALS),'zadna data VALS');
            plotSetup = {};
            if ~exist('kategorie','var') || isempty(kategorie) , kategorie = 1:size(obj.VALS,2); end %muzu chtit jen nektere kategorie
            if ~exist('outputDir','var')
                plotSetup.outputDir = 'd:\eeg\motol\CBrainPlot\';    
            else
                plotSetup.outputDir = outputDir;
            end
            if ~exist('signum','var'), signum = 0; end; %defaultni je rozdil kladny i zaporny
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
            
            tic; %zadnu meric cas
            for interval = 1:size(obj.VALS,1) 
                for kat = kategorie
                    if signum > 0 
                        iV = obj.VALS{interval,kat} > 0; %jen kladne rozdily
                    elseif signum <0 
                        iV = obj.VALS{interval,kat} < 0; %jen zaporne rozdily
                    else
                        iV = true(size(obj.VALS{interval,kat})); %vsechny rozdily
                    end
                    
%                     if kat <= numel(obj.katstr) %pokud puvodni kategorie z intervalyResp
                    katname = obj.katstr{kat};
%                     elseif kat==size(obj.VALS,2) %posledni kategorie se vsemi kanaly
%                         katname = 'AllEl';                    
%                     end
                    
                    if strcmp(plotSetup.figureVisible,'off')
                        disp('figures invisible');
                    end
                    plotSetup.figureNamePrefix = [ obj.testname '_' mat2str(obj.intervals(interval,:))  '_' katname '_' num2str(signum) '_NOnames'];
                    if numel(obj.VALS{interval,kat}(iV)) > 0                        
                        disp(plotSetup.figureNamePrefix);
                        vals_channels = obj.VALS{interval,kat}(iV); %parametr  main_brainPlot
                        if signum ~= 0
                            vals_channels = vals_channels*signum; %u zapornych hodnot prehodim znamenko
                        end
                        mni_channels = obj.MNI{interval,kat}(iV);                                                                         
                        names_channels = []; 
                          
                        
                        %nejdriv vykreslim bez popisku elektrod
                        brainsurface = main_brainPlot(vals_channels,mni_channels,names_channels,brainsurface,plotSetup);  %#ok<PROPLC>
                        %volam Jirkuv skript, vsechny ty promenne predtim jsou do nej
                        if isempty(obj.brainsurface)
                            obj.brainsurface = brainsurface; %#ok<PROPLC> %ulozim si ho pro dalsi volani
                        end
                        
                        %a pak jeste s popisy elektrod
                        plotSetup.figureNamePrefix = [obj.testname '_' mat2str(obj.intervals(interval,:)) '_' katname '_' num2str(signum) '_names'];
                        disp(plotSetup.figureNamePrefix);
                        names_channels = obj.NAMES{interval,kat}; 
                        brainsurface = main_brainPlot(vals_channels,mni_channels,names_channels,brainsurface,plotSetup);    %#ok<PROPLC>  
                        
                    else
                        disp(['zadne hodnoty pro ' plotSetup.figureNamePrefix ' - neukladam ']);
                    end
                end
            end
            toc; %ukoncim mereni casu a vypisu
            if hybernovat
                system('shutdown -h')  %#ok<UNRCH>
            elseif vypnout            
                system('shutdown -s') %#ok<UNRCH>
            end
        end
    end
    methods (Static,Access = public)
        function PAC = StructFind(struktura,label,testname,reference)
            %najde pacienty, jejich headery obsahuji mozkovou strukturu
            %struktura je nazev struktury podle atlas napriklad hippo, label je kratky nazev podle martina, napriklad hi
            if ~exist('label','var'),    label = struktura; end %defaultni test
            if ~exist('testname','var') || isempty(testname), testname = 'aedist'; end %defaultni test
            if ~exist('reference','var') || isempty(reference), reference = []; end %defaultni test
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
                index = [];
                labels = lower({H.channels(ii).neurologyLabel}');
                for jj = 1:size(label,2)                    
                    indexjj =  find(~cellfun('isempty',strfind(labels,lower(label{jj}))))'; %rozepsal jsem, aby se to dalo lepe debugovat
                    index = [index indexjj];  %#ok<AGROW>
                    % 3.5.2018 nejakym zahadnym zpusobem funguje hledani pomoci strfind ve sloupci a ne v radku. 
                    % Proto nejdriv prehodim pomoci ' na sloupec a pak zase na radek
                end
                iiBA = ~cellfun(@isempty,{H.channels.ass_brainAtlas}); %neprazdne cells
                iiCM = ~cellfun(@isempty,{H.channels.ass_cytoarchMap}); %neprazdne cells
                for jj = 1:size(struktura,2)
                    index = [ index find(~cellfun('isempty',strfind(lower({H.channels(iiCM).ass_cytoarchMap}),lower(struktura{jj}))))]; %#ok<AGROW>
                    index = [ index find(~cellfun('isempty',strfind(lower({H.channels(iiBA).ass_brainAtlas}),lower(struktura{jj}))))]; %#ok<AGROW>
                end
                index = union(index,[]); %vsechny tri dohromady
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
                    iPAC = iPAC + 1;
                end
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
    end
    
end

