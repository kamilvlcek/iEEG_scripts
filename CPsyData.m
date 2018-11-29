classdef CPsyData < matlab.mixin.Copyable %je mozne kopirovat pomoci E.copy();
    %CPSYDATA Trida na praci s behavioralnimi daty z psychopy
    %   vytvorena pomoci ppa_data.m, aedist_data.m aj
    % Kamil Vlcek, FGU AVCR, since 2016 04
    
    properties  (Access = public)
        P; %psychopy behavioural data
        fhR; %figure handle from PlotResponses
        warning_rt=false; %jestli uz byl warning o reakcnich casech
        testname; %jmeno testu, ze ktereho jsou data
        blocks; %struktura, ktera uklada vysledky GetBlocks, kvuli uspore casu
    end
    
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS
        function obj = CPsyData(psy)
            %konstruktor
            obj.P = psy;
            if ~isfield(psy,'pacientid'), obj.P.pacientid = ''; end
            if ~isfield(psy,'eegfile'), obj.P.eegfile = ''; end
            obj.DoplnZpetnavazba();
            obj.testname = obj.GetTestName(inputname(1));
        end
        function [obj] = SubjectChange(obj,~)
            %jen prazdna funkce, ale pretezuju ji v CPsyDataMulti, a naplnuju obsahem
            % v tehle trida mam jen jeden subjekt
        end
        function [testname2] = GetTestName(obj,testname)  
            %zjisti jmeno testu, ze ktereho jsou data            
            if strcmp(testname,'aedist') || strcmp(testname,'menrot') || strcmp(testname,'ppa')
                testname2 = testname;
            elseif strcmp(obj.P.strings.podminka{1,1},'cervena')
                testname2 = 'aedist';
            elseif strcmp(obj.P.strings.podminka{1,1},'vy-2D')
                testname2 = 'menrot';
            elseif strcmp(obj.P.strings.podminka{1,1},'Ovoce')
                testname2 = 'ppa';
            else
                testname2 = '';
            end
        end
        function rt = ReactionTime(obj,nokategories)
            %vrati matici vsech reakcnich casu roztridenych do sloupcu podle kategorii podnetu            
            %kde je min hodnot v kategorii, zarovnano pomoci NaN, takze se musi pouzit nanmean
            %reakcni casy se berou ze synchronizacnich pulzu
            % pokud nokategories= 1, tak jen jeden sloupec reakcnich casu, bez rozliseni na kategorie
            if ~exist('nokategories','var') || nokategories == 0
                kat = unique(obj.P.data(:,obj.P.sloupce.kategorie)); %ciselne vyjadreni kategorie podnetu 0-n
                rt = nan(size(obj.P.data,1),numel(kat)); %pole kam budu ukladat reakcni casy v sekundach

                for k = kat' %funguje jen pro radky, kat je sloupec
                    ikat = obj.P.data(:,obj.P.sloupce.kategorie)==k;
                    rt(1:sum(ikat),k+1) = 24*3600*(obj.P.data(ikat,obj.P.sloupce.ts_odpoved) - obj.P.data(ikat,obj.P.sloupce.ts_podnet));
                end
                rt = rt(any(~isnan(rt),2),:); % necha jen radky, kde je nejake ~NaN cislo           
            else
                rt = 24*3600*(obj.P.data(:,obj.P.sloupce.ts_odpoved) - obj.P.data(:,obj.P.sloupce.ts_podnet));
            end
            if any(sum(rt<0)>=1) && isempty(obj.warning_rt)
                warning('Pozor: %d zaporne reakcni casy',sum(rt<0));
                obj.warning_rt = true;
            end
        end
        
        function isi = InterStimulusInterval(obj)
            %vraci matici vsech interstimulus intervalu v sekundach
            isi = 24*3600*(obj.P.data(2:end,obj.P.sloupce.ts_podnet) - obj.P.data(1:end-1,obj.P.sloupce.ts_podnet));
        end
        
        function ts_podnety = TimeStimuli(obj,response)
            %vraci matici timestampu vsech podnetu/odpovedi, pokud response=1 vracim odpovedi
            if exist('response','var') && response==1
                ts_podnety = obj.P.data(:,obj.P.sloupce.ts_odpoved); %vratim timestampy odpovedi
            else
                ts_podnety = obj.P.data(:,obj.P.sloupce.ts_podnet); %vratim timestampy podnetu
            end
        end
                       
        function [kat, katnum] = Category(obj,event)
            %vrati 1 retezec s popisem kategorie eventu a 2. cislo kategorie
            katnum = obj.P.data(event,obj.P.sloupce.kategorie); %cislo kategorie
            kat = obj.P.strings.podminka{katnum+1};            
        end 
        
        function [katnum, katstr] = Categories(obj,tisk)
            %vraci cisla vsech kategorii
            if ~exist('tisk','var'), tisk = 0; end %dfaultne netisknu kategorie
            katnum = cell2mat(obj.P.strings.podminka(:,2))'; %kategorie v radku
            katstr = cell(size(katnum));
            for k = 1:numel(katnum)
                katstr{k} = obj.P.strings.podminka{k};
                if tisk
                    disp([ num2str(katnum(k)) ': ' katstr{k}]);
                end
            end                
        end
        function [kat] = CategoryName(obj,katnum,concat)
            %vraci jmeno kategorie z cisla, jmeno kategorie se pocita od 0
            %muze byt i vic cisel kategorii
            % pak se jmena spoji do jednoho retezce pomoci parametru concat (default +)
            % pokud je concat prazdne, vrati se cell array
            
            if ~exist('concat','var'), concat = '+'; end %defaulte se vice jmen kategorii spojuje pomoci +
            if numel(katnum) == 1
                kat = obj.P.strings.podminka{katnum+1};
            elseif ~isempty(concat)
                kat = '';
                for k = katnum
                    if numel(kat) == 0
                        kat = obj.P.strings.podminka{k+1};
                    else
                        kat = [ kat '+' obj.P.strings.podminka{k+1}]; %#ok<AGROW>
                    end
                end
            else %chci jmena kategorii vratit jako cell array
                kat = cell(1,numel(katnum));
                for k  = 1:numel(katnum)
                    if numel(cellval(katnum,k))>1
                        katkat = cell(numel(katnum{k}),1);
                        for kk = 1:numel(katnum{k})
                            katkat{kk} = obj.P.strings.podminka{katnum{k}(kk)+1};
                        end
                        kat{k} = cell2str(katkat);
                    else
                        kat{k} = obj.P.strings.podminka{cellval(katnum,k)+1};
                    end
                end
            end
        end
        function [katnum] = CategoryNum(obj,katname)
            %vraci cislo kategorie podle jmena
            %muze byt i nekolik jmen kategorii v cell array
            katnum = nan(size(katname));
            for k = 1:numel(katname)
                ikat = not(cellfun('isempty',strfind(obj.P.strings.podminka(:,1),katname{k}))); %logical index, ale s jen jednim true
                katnum(k) = obj.P.strings.podminka{ikat,2};
            end
        end
        function [blocks, srate, test, kategorie]= GetBlocks(obj)
            %vraci promenne pro bloky o stejne kategorii
            b1 = [find(obj.P.data(2:end,obj.P.sloupce.kategorie) ~= obj.P.data(1:end-1,obj.P.sloupce.kategorie)); size(obj.P.data,1)]; %konce bloku
            if min(diff(b1)) == 1
                b1 = (1:size(obj.P.data,1))'; %pro PPA test - pokud jsou nektere bloky velikosti 1, tak je udelam vsechy tak mal
            end %jinak to dela problemy pri prohazeni a pocitam srate
            b0 = [ 1; (b1(1:end-1)+1)]; %zacatky bloku
            kategorie = obj.P.data(b0,obj.P.sloupce.kategorie);            
            blocks = [b0 b1];
            if isempty(obj.blocks) || size(blocks,1) ~= numel(obj.blocks.srate)
                obj.blocks.srate = zeros(numel(b0),1); % prumerna uspesnost za blok
                obj.blocks.test = ones(numel(b0),1); %jestli by blok testovy, tady nastavim ze vsechny
                for block = 1:numel(b0)
                    obj.blocks.srate(block,1)=mean(obj.P.data(b0(block) : b1(block),obj.P.sloupce.spravne));                       
                    obj.blocks.test(block) = obj.blocks.test(block) - max(obj.P.data(b0(block) : b1(block),obj.P.sloupce.zpetnavazba));                        
                        %test=1 vsechny epochy v bloku testove, test=0 alespon jedna epocha treningova
                end                
            end
            srate = obj.blocks.srate;
            test = obj.blocks.test;
            
        end
        
        function [resp,rt,kat,test] = GetResponses(obj)
            %vraci odpovedi cloveka - spravne/spatne, reakcni cas, kategorie 
            S = obj.P.sloupce;
            resp = obj.P.data(:,S.spravne); % spravnost odpovedi
            rt = obj.P.data(:,S.rt); % reakcni casy
            kat = obj.P.data(:,S.kategorie); %kategorie
            test = obj.P.data(:,S.zpetnavazba)==0; %vratim index testovych trialu
        end
        
        function chyby = GetErrorTrials(obj)
            %vrati pole indikujici chybu/vyrazeni pro kazdy trial=radky - sloupce: jednotlive chyby, chybne bloky, treningovy trial, prilis kratky reakcni cas
            %chybny blok ma < 75% uspesnost
            S = obj.P.sloupce;          
            
            chyby = zeros(size(obj.P.data,1),4); %ctyri sloupce - chybne trials a chybne bloky, trening, prilis rychle reakcni casy
            chyby(:,1) = obj.P.data(:,S.spravne)==0; %pro PPA jsou vsechny ovoce spatne. Sloupec spravne je u ovoce vzdy 0, chyba v PHP asi
            rt = obj.ReactionTime(1); %reakcni casy podle Sychropulsu
            rtPsy = obj.P.data(:,S.rt); %reakcni cas podle psychopy
            chyby(:,4) = rt(:,1) < 0.1 | (rtPsy(:,1) < 0.1 & rtPsy(:,1) > 0);  %v PPA clovek nereaguje spravne, takze 0 jako cas odpovedi me nezajima 
                    %chyba, pokud je reakcni cas prilis kratky (0 v PsychoPy znamena, ze nereagoval, to je taky chyba)
            [blocks,srate,blocktest]=obj.GetBlocks();            %#ok<PROP>
            for b = 1:size(blocks,1) %#ok<PROP>
                if srate(b) < 0.75  %chybny blok
                    chyby(blocks(b,1) : blocks(b,2) , 2) = ones( blocks(b,2) - blocks(b,1) +1,1);%#ok<PROP> %vyplnim jednickami
                end 
                if  blocktest(b)==0 %treningovy blok
                    chyby(blocks(b,1) : blocks(b,2) , 3) = ones( blocks(b,2) - blocks(b,1) +1,1); %#ok<PROP> %vyplnim jednickami
                end
            end
        end
        function opak = GetOpakovani(obj)
            %vraci cislo opakovani obrazku pro kazdou epochu
            S = obj.P.sloupce;
            if isfield(S,'opakovani_obrazku')
                opak = obj.P.data(:,S.opakovani_obrazku); % ppa test
            elseif isfield(S,'opakovani_obrazku')
                opak = obj.P.data(:,S.opakovani);  %aedist test
            else
                opak = zeros(size(obj.P.data,1),1); %opakovani 0
            end
        end
        function obj = Cond2Epochs(obj)
            %predela conditions na epochy, kvuli ITPC
            katnum = obj.Categories();
            data2 = zeros(numel(katnum),size(obj.P.data,2));
            for k = 1:numel(katnum)
                idata = obj.P.data(:,obj.P.sloupce.kategorie)==katnum(k) & obj.P.data(:,obj.P.sloupce.zpetnavazba)==0 ;
                data2(k,:) = [0 0 1 mean(obj.P.data(idata,4)) 0 0 katnum(k) min(obj.P.data(idata,8)) min(obj.P.data(idata,9))];
            end
            obj.P.data = data2;
        end
        %% PLOT FUNCTIONS
        function [obj, chyby] = PlotResponses(obj)
            %nakresli graf rychlosti vsech odpovedi a bloku, vcetne chyb a uspesnosti za blok
            S = obj.P.sloupce;
            test = obj.P.data(:,:); %vyberu vsechny trialy
            if isfield(obj,'fh') && ishandle(obj.fhR) 
                figure(obj.fhR); %pokud uz graf existuje, nebudu tvorit znova
                clf; %smazu aktualni figure
            else                
                obj.fhR = figure('Name','Responses');
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
                    if procent < 75 || size(blocks,1)<100 %vypisuju uspesnost u vsech bloku jen pokud jich je min nez 100; jinak jen chybne
                        text( blocks(b,1) , -0.1, [ num2str(procent) '%'],'FontSize',8,'Color',color);
                    end
                elseif size(blocks,1) < 100
                    text( blocks(b,1) , -0.1, '*','FontSize',8,'Color','k'); %treningove bloky jako hvezdicky
                end
            end
            resp = obj.P.data(:,S.spravne); %vratim i reakcni casy            
            text(0,-0.2, ['pocet chyb: ' num2str(numel(resp)-sum(resp)) ]);
        end
        function [obj] = PlotITI(obj)
            %graf intervalu mezi obrazky
            %TODO korelace mezi RT podle Synchropulsu a PsychoPy
            figure('Name','ITI');
            subplot(1,2,1);
            isi = obj.InterStimulusInterval();
            plot(isi,'.');
            ylim([1 prctile(isi,99)]); %rozsah do 99  percentilu - odstranim odlehle hodnoty
            title('ITI');
            ylabel('sec');
            
            subplot(1,2,2);
            rt = obj.ReactionTime(); %ctyri sloupce reakcnich casu podle synchropulsu
            numkatvals = size(rt,1);
            numkats = size(rt,2);
            rt = rt(:);      % jeden sloupec, v poradi ovoce, 
            plot(rt,'o');
            for j = 1:numkats
                line([numkatvals*j numkatvals*j],[0.4 1.4]);
                text(numkatvals*(j-0.5),0.1,obj.P.strings.podminka{j},'Color','blue');
            end
            title('RT');
            ylabel('sec');
            hold on;
            rtPsy = obj.P.data(:,obj.P.sloupce.rt);
            plot(rtPsy,'xm'); %reakcni casy podle PsychoPy - tohle neni serazene podle kategorii je tam 0 pokud zadna reakce
            
            %druhy graf - 21.6.2017
            figure('Name','ITI - Synchropuls vs PsychoPy');
            rt = obj.ReactionTime(1);
            plot(rt,rtPsy,'o');
            xlabel('RT Synchropulse');
            ylabel('RT PsychoPy');
       
        end
        function [id] = PacientID(obj)
            %funkce ktera bude vracet id pacienta
            id = obj.P.pacientid;
            if length(id) <= 5 %pokud id obsahuje jen napr p132
                podtrzitko = strfind(obj.P.eegfile,'_');
                id = [id '-' obj.P.eegfile(1 : podtrzitko-1)]; %pritam jeste zacatke jmena eeg souboru, naprilad VT18
             end
        end
        function [obj]= RemoveTraining(obj)
            %odstrani treningove trialy - ty kde byla zpetna vazba
            %since 15.12.2017
            test = obj.P.data(:,obj.P.sloupce.zpetnavazba)==0;
            obj.P.data = obj.P.data(test,:);
        end
    end
    methods  (Access = protected)
        function [obj] = DoplnZpetnavazba(obj)
            %doplni sloupec zpetnavazba, pokud neexistuje a naplni ho nulama
            if isstruct(obj.P) && ~isfield(obj.P.sloupce,'zpetnavazba')                
                obj.P.data = [ obj.P.data zeros(size(obj.P.data,1),1)]; %doplnim dalsi sloupec do dat
                obj.P.sloupce.zpetnavazba = size(obj.P.data,2); %pojmenuju ho zpetnavazba
            end
        end
    end
    methods (Static,Access = public)
        function [opakname] = OpakovaniName(opakovani)
            %vrati jmeno opakovani. predpoklada jako parametr array
            opakname = ['Opak' num2str(opakovani)];
        end
    end
    
end

