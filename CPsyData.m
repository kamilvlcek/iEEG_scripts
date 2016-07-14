classdef CPsyData < handle
    %CPSYDATA Trida na praci s behavioralnimi daty z psychopy
    %   vytvorena pomoci ppa_data.m, aedist_data.m aj
    % Kamil Vlcek, FGU AVCR, since 2016 04
    
    properties  (Access = public)
        P; %psychopy behavioural data
        fhR; %figure handle from PlotResponses
    end
    
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS
        function obj = CPsyData(psy)
            %konstruktor
            obj.P = psy;
        end
        function rt = ReactionTime(obj)
            %vrati matici vsech reakcnich casu roztridenych do sloupcu podle kategorii podnetu
            %kde je min hodnot v kategorii, zarovnano pomoci NaN, takze se musi pouzit nanmean
            kat = unique(obj.P.data(:,obj.P.sloupce.kategorie)); %ciselne vyjadreni kategorie podnetu 0-n
            rt = nan(size(obj.P.data,1),numel(kat)); %pole kam budu ukladat reakcni casy v sekundach
            
            for k = kat' %funguje jen pro radky, kat je sloupec
                ikat = obj.P.data(:,obj.P.sloupce.kategorie)==k;
                rt(1:sum(ikat),k+1) = 24*3600*(obj.P.data(ikat,obj.P.sloupce.ts_odpoved) - obj.P.data(ikat,obj.P.sloupce.ts_podnet));
            end
            rt = rt(any(~isnan(rt),2),:); % necha jen radky, kde je nejake ~NaN cislo           
        end
        
        function isi = InterStimulusInterval(obj)
            %vraci matici vsech interstimulus intervalu v sekundach
            isi = 24*3600*(obj.P.data(2:end,obj.P.sloupce.ts_podnet) - obj.P.data(1:end-1,obj.P.sloupce.ts_podnet));
        end
        
        function ts_podnety = TimeStimuli(obj)
            %vraci matici timestampu vsech podnetu
            ts_podnety = obj.P.data(:,obj.P.sloupce.ts_podnet);
        end
        
        function [kat katnum] = Category(obj,event)
            %vrati 1 retezec s popisem kategorie eventu a 2. cislo kategorie
            katnum = obj.P.data(event,obj.P.sloupce.kategorie); %cislo kategorie
            kat = obj.P.strings.podminka{katnum+1};            
        end 
        
        function kat = Categories(obj)
            %vraci cisla vsech kategorii
            kat = cell2mat(obj.P.strings.podminka(:,2))'; %kategorie v radku
        end
        function [kat] = CategoryName(obj,katnum)
            %vraci jmeno kategorie z cisla
            kat = obj.P.strings.podminka{katnum+1};
        end
        
        function [blocks, srate, test, kategorie]= GetBlocks(obj)
            %vraci promenne pro bloky o stejne kategorii
            b1 = [find(obj.P.data(2:end,obj.P.sloupce.kategorie) ~= obj.P.data(1:end-1,7)); size(obj.P.data,1)]; %konce bloku
            b0 = [ 1; (b1(1:end-1)+1)]; %zacatky bloku
            kategorie = obj.P.data(b0,obj.P.sloupce.kategorie);
            srate = zeros(numel(b0),1); % prumerna uspesnot za blok
            test = ones(numel(b0),1);
            for block = 1:numel(b0)
                srate(block,1)=mean(obj.P.data(b0(block) : b1(block),obj.P.sloupce.spravne));    
                test(block) = test(block) - mean(obj.P.data(b0(block) : b1(block),obj.P.sloupce.zpetnavazba));    
            end
            blocks = [b0 b1];
        end
        
        function [resp,rt,kat,test] = GetResponses(obj)
            %vraci odpovedi cloveka - spravne/spatne, reakcni cas, kategorie 
            S = obj.P.sloupce;
            resp = obj.P.data(:,S.spravne); %vratim i reakcni casy
            rt = obj.P.data(:,S.rt); %vratim i reakcni casy
            kat = obj.P.data(:,S.kategorie); %vratim i reakcni casy
            test = obj.P.data(:,S.zpetnavazba)==0; %vratim index testovych trialu
        end
        
        function chyby = GetErrorTrials(obj)
            %vrati pole indikujici chybu/vyrazeni pro kazdy trial=radky - podle sloupcu  chyby, chybne bloky a treningovy trial
            %chybny blok ma < 75% uspesnost
            S = obj.P.sloupce;
            chyby = zeros(size(obj.P.data,1),3); %tri sloupce - chybne trials a chybne bloky, trening
            chyby(:,1) = obj.P.data(:,S.spravne)==0;
            [blocks,srate,blocktest]=obj.GetBlocks();
            for b = 1:size(blocks,1)
                if srate(b) < 0.75  %chybny blok
                    chyby(blocks(b,1) : blocks(b,2) , 2) = ones( blocks(b,2) - blocks(b,1) +1,1); %vyplnim jednickami
                end 
                if  blocktest(b)==0 %treningovy blok
                    chyby(blocks(b,1) : blocks(b,2) , 3) = ones( blocks(b,2) - blocks(b,1) +1,1); %vyplnim jednickami
                end
            end
        end
        
        %% PLOT FUNCTIONS
        function [chyby] = PlotResponses(obj)
            %nakresli graf rychlosti vsech odpovedi a bloku, vcetne chyb a uspesnosti za blok
            S = obj.P.sloupce;
            test = obj.P.data(:,:); %vyberu vsechny trialy
            if isempty(obj.fhR)
                obj.fhR = figure('Name','Responses');
            else
                figure(obj.fhR); %pokud uz graf existuje, nebudu tvorit znova
                clf; %smazu aktualni figure
            end
            plot(test(:,S.rt),'-o');
            hold on;
            treningtrials = find(obj.P.data(:,S.zpetnavazba)==1);            
            bar(treningtrials,test(treningtrials,S.rt),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]); %oznacim treningove pokusy
            
            chyby = find(test(:,S.spravne)==0);
            plot(chyby,test(chyby,S.kategorie),'sr','MarkerSize',10,'MarkerFaceColor','y'); %vykreslim chyby na vysce kategorie
            plot(chyby,test(chyby,S.rt),'or','MarkerFaceColor','y'); %vykreslim reakcni cas cervene
            
            plot(test(:,S.kategorie),'g','LineWidth',2); %vykreslim kategorie podnetu, jako zelenou caru
            for j = 1:size(obj.P.strings.podminka,1)
                text (10,obj.P.strings.podminka{j,2}+0.05,obj.P.strings.podminka{j,1},'Color','r','FontSize',15);
            end
                        
            [blocks,srate,blocktest]=obj.GetBlocks();
            for b = 1:size(blocks,1)
                if blocktest(b) == 1 %pokud se jedna o testovy blok
                    procent = round(srate(b)*100);
                    if procent < 75, color = 'r'; else color = 'b'; end %cervenou barvou, pokud je uspesnost pod 75%
                    text( blocks(b,1) , -0.1, [ num2str(procent) '%'],'FontSize',8,'Color',color);
                else
                    text( blocks(b,1) , -0.1, '*','FontSize',8,'Color','k'); %treningove bloky jako hvezdicky
                end
            end
        end
    end
    
end

