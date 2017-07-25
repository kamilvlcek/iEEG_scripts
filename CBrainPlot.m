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
            if strcmp(testname,'aedist')
                pacienti = pacienti_aedist(); %nactu celou strukturu pacientu    
            end
            obj.intervals = intervals;            
            elcount = []; %jen inicializace            
            P = {}; M = {}; N = {}; %jen inicializace
            for p = 1:numel(pacienti)
                disp(['***   ' pacienti(p).folder '   ***']);
                E = pacient_load(pacienti(p).folder,'aedist',filename);
                [prumery, MNI,names,~,katstr] = E.IntervalyResp( intervals,[],0); %#ok<PROP>    %no figure       
                clear E;
                if p==1
                    obj.katstr = katstr; %#ok<PROP>
                    elcount = zeros(size(prumery,2),size(prumery,3)); %pocet elektrod pro kazdy casovy interval a kategorii - interval x kategorie
                    P = cell([numel(pacienti),size(prumery,2),size(prumery,3)]); % souhrnne prumery pro vsechny pacienty + + jejich kombinace
                    M = cell([numel(pacienti),size(prumery,2),size(prumery,3)]); % souhrnne MNI koordinaty pro vsechny pacienty
                    N = cell([numel(pacienti),size(prumery,2),size(prumery,3)]); % souhrnne names pro vsechny pacienty
                end
                for interval = 1:size(prumery,2)
                    for kat = 1:size(prumery,3)
                        if kat <= numel(obj.katstr)
                            ip = prumery(:,interval, kat) > 0; % index pro prumery MNI i names, chci jen kladne odpovedi
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
        function PlotBrain3D(obj,kategorie)
            assert(~isempty(obj.VALS),'zadna data z IntervalyResp');
            if ~exist('kategorie','var'), kategorie = 1:size(obj.VALS,2); end %muzu chtit jen nektere kategorie
            if ~isempty(obj.brainsurface)
                brainsurface = obj.brainsurface; %#ok<PROP>
            end
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
                            obj.brainsurface = brainsurface; %#ok<PROP> %ulozim si ho pro dalsi volani
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
        end
    end
    
end

