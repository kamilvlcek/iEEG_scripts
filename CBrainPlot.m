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
    end
    
    methods (Access = public)        
        function [obj] = IntervalyResp(obj,testname,intervals,filename)
            %vola postupne pro vsechny pacienty E.IntervalyResp a uklada vysledky
            %vyradi vsechny kontakty bez odpovedi nebo se zapornou odpovedi
            %spoji vsechno dohromady
            %vrati vysledky ve formatu pro SEEE-vizualization
            %napr CB.IntervalyResp('aedist',[0.2 0.8],'AEdist CHilbert 50-120 refBipo Ep2017-11_CHilb.mat');
            if strcmp(testname,'aedist')
                pacienti = pacienti_aedist(); %nactu celou strukturu pacientu    
            end
            obj.intervals = intervals;            
            elcount = []; %jen inicializace            
            P = {}; M = {}; N = {}; %jen inicializace
            for p = 1:numel(pacienti) % cyklus pacienti
                disp(['***   ' pacienti(p).folder '   ***']);
                E = pacient_load(pacienti(p).folder,'aedist',filename); %nejspis objekt CHilbert, pripadne i jiny
                if isempty(E)
                    disp('no data');
                    continue;
                end
                [prumery, MNI,names,~,katstr] = E.IntervalyResp( intervals,[],0); %#ok<PROP>    %no figure, funkce z CiEEGData       
                clear E;
                if p==1
                    obj.katstr = katstr; %#ok<PROP>
                    elcount = zeros(size(prumery,2),size(prumery,3)); %pocet elektrod pro kazdy casovy interval a kategorii - interval x kategorie
                    P = cell([numel(pacienti),size(prumery,2),size(prumery,3)]); % souhrnne prumery pro vsechny pacienty + + jejich kombinace
                    M = cell([numel(pacienti),size(prumery,2),size(prumery,3)]); % souhrnne MNI koordinaty pro vsechny pacienty
                    N = cell([numel(pacienti),size(prumery,2),size(prumery,3)]); % souhrnne names pro vsechny pacienty
                end
                for interval = 1:size(prumery,2) % cyklus intervaly
                    for kat = 1:size(prumery,3) % cyklus kategorie podnetu
                        if kat <= numel(obj.katstr)
                            ip = prumery(:,interval, kat) > 0; % index pro prumery, MNI i names, chci jen kladne odpovedi
                        else
                            ip = prumery(:,interval, kat) ~= 0; % pro vyssi kategorie, ktere jsou rozdily odpovedi, chci i zaporny rozdil ; aby tam neco bylo 
                        end
                        P{p,interval,kat}=prumery(ip,interval, kat); %#ok<AGROW>
                        M{p,interval,kat}=MNI(ip); %#ok<AGROW,PROP>
                        N{p,interval,kat}= strcat(cellstr(repmat([pacienti(p).folder '_'],sum(ip),1)),names(ip)); %#ok<AGROW>
                        elcount(interval,kat) = elcount(interval,kat) + sum(ip); %#ok<AGROW>
                    end                       
                end                
            end
            %ted z P M a N rozdelenych po pacientech udelam souhrnna data
            obj.VALS = cell(size(elcount)); %souhrnne prumery 
            obj.MNI = cell(size(elcount)); 
            obj.NAMES = cell(size(elcount));       
            for interval = 1:size(prumery,2) 
                for kat = 1:size(prumery,3)                    
                      obj.VALS{interval,kat} = zeros(elcount(interval,kat),1);
                      obj.MNI{interval,kat} = struct('MNI_x',{},'MNI_y',{},'MNI_z',{});
                      obj.NAMES{interval,kat} = cell(elcount(interval,kat),1);
                      iVALS = 1;
                      for p = 1:numel(pacienti) 
                          n = numel(P{p,interval,kat});
                          obj.VALS{interval,kat} (iVALS:iVALS+n-1)=P{p,interval,kat};
                          obj.MNI{interval,kat}  (iVALS:iVALS+n-1)=M{p,interval,kat};
                          obj.NAMES{interval,kat}(iVALS:iVALS+n-1)=N{p,interval,kat};
                          iVALS = iVALS + n;
                      end
                end
            end            
        end
        function obj = ImportData(obj,CB)
            %vlozi data, ktera jsem vytvoril pomoci CHilbert.ExtractBrainPlotData
            obj.VALS = CB.VALS;
            obj.MNI = CB.MNI;
            obj.NAMES = CB.NAMES;
            obj.katstr = CB.katstr;
            obj.intervals = CB.intervals;       
        end
        function PlotBrain3D(obj,kategorie)
            assert(~isempty(obj.VALS),'zadna data VALS');
            if ~exist('kategorie','var'), kategorie = 1:size(obj.VALS,2); end %muzu chtit jen nektere kategorie
            if ~isempty(obj.brainsurface)
                brainsurface = obj.brainsurface;  %#ok<PROPLC>
            end
            hybernovat = 0; %jestli chci po konci skriptu pocitac uspat - ma prednost
            vypnout = 0;  %jestli chci po konci skriptu pocitac vypnout (a nechci ho hybernovat) 
            kombinace = [2 1; 3 1 ; 3 2 ]; %kombinace kategorii
            figureVisible = 'off';   %nechci zobrazovat obrazek 
            tic; %zadnu meric cas
            for interval = 1:size(obj.VALS,1) 
                for kat = kategorie                   
                    if kat <= numel(obj.katstr)
                        katname = obj.katstr{kat};
                    else
                        katname = [obj.katstr{kombinace(kat-3,1)} '-' obj.katstr{kombinace(kat-3,2)}];
                    end
                    figureNamePrefix = ['AEdist_' mat2str(obj.intervals(interval,:)) '_' katname '_NOnames'];
                    if strcmp(figureVisible,'off')
                        disp('figures invisible');
                    end
                    if numel(obj.VALS{interval,kat}) > 0
                        disp(figureNamePrefix);
                        vals_channels = obj.VALS{interval,kat}; %#ok<NASGU>
                        mni_channels = obj.MNI{interval,kat};   %#ok<NASGU>                    
                        outputDir = 'd:\eeg\motol\CBrainPlot\';    %#ok<NASGU>                           
                        names_channels = []; %#ok<NASGU> 
                        %nejdriv vykreslim bez popisku elektrod
                        main_brainPlot; %volam Jirkuv skript, vsechny ty promenne predtim jsou do nej
                        if isempty(obj.brainsurface)
                            obj.brainsurface = brainsurface; %#ok<PROPLC> %ulozim si ho pro dalsi volani
                        end
                        
                        %a pak jeste s popisy elektrod
                        figureNamePrefix = ['AEdist_' mat2str(obj.intervals(interval,:)) '_' katname '_names'];
                        disp(figureNamePrefix);
                        names_channels = obj.NAMES{interval,kat}; %#ok<NASGU>
                        main_brainPlot;                     else
                        disp(['zadne hodnoty pro ' figureNamePrefix ' - neukladam ']);
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
            if ~exist('reference','var'), reference = []; end %defaultni test
            if ischar(struktura), struktura = {struktura}; end %prevedu na cell array
            if ischar(label), label = {label}; end %prevedu na cell array
            if strcmp(testname,'aedist')
                pacienti = pacienti_aedist(); %nactu celou strukturu pacientu    
            elseif strcmp(testname,'menrot')
                pacienti = pacienti_menrot(); %nactu celou strukturu pacientu    
            else
                error('nezname jmeno testu');
            end
            PAC = {};
            iPAC = 1;
            for p = 1:numel(pacienti)
                disp(['* ' pacienti(p).folder ' - ' pacienti(p).header ' *']);
                hfilename = ['D:\eeg\motol\pacienti\' pacienti(p).folder '\' pacienti(p).header];
                if exist(hfilename,'file')==2
                    load(hfilename);
                else
                    disp(['header ' hfilename ' neexistuje']);
                end               
                if ~isempty(reference)
                    CH = CHHeader(H);
                    CH.RejectChannels( pacienti(p).rjch); %musim vyradit vyrazene kanaly, protoze ty se vyrazuji v bipolarni referenci
                    CH.ChangeReference(reference); %nove od 18.1.2018
                    H = CH.H;
                end
                ii = ~cellfun(@isempty,{H.channels.neurologyLabel}); %neprazdne cells
                index = [];
                for jj = 1:size(label,2)
                    index = [index find(~cellfun('isempty',strfind(lower({H.channels(ii).neurologyLabel}),lower(label{jj}))))];  %#ok<AGROW>
                end
                iiBA = ~cellfun(@isempty,{H.channels.ass_brainAtlas}); %neprazdne cells
                iiCM = ~cellfun(@isempty,{H.channels.ass_cytoarchMap}); %neprazdne cells
                for jj = 1:size(struktura,2)
                    index = [ index find(~cellfun('isempty',strfind(lower({H.channels(iiCM).ass_cytoarchMap}),lower(struktura{jj}))))]; %#ok<AGROW>
                    index = [ index find(~cellfun('isempty',strfind(lower({H.channels(iiBA).ass_brainAtlas}),lower(struktura{jj}))))]; %#ok<AGROW>
                end
                index = union(index,[]); %vsechny tri dohromady
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
    end
    
end

