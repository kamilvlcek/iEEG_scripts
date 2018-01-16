classdef CHilbertMulti < CHilbert 
    %CHILBERTMULTI zpracovava data z nekolika CHilbert vyexportovanych data
    %   nacita data vyprodukovana pomoci CHilbert.ExtractData
    
    properties
        orig; %struktura s puvodnimi daty 
        filenames; %cell array se jmeny souboru
        blokyprehazej; %matrix s originalnimi cisly bloku a cisly bloku k prehazeni z noveho souboru. Meni se pro kazdy pridany soubor
           %plni se v GetEpochData
    end
    
    methods
        function obj = CHilbertMulti(filenames)  
            obj.ImportExtract(filenames);
        end
        function obj = ImportExtract(obj,filenames)
            for fileno = 1:numel(filenames)
                filename = filenames{fileno}; %cell array, zatim to musi byt plna cesta                
                if exist(filename,'file') 
                    disp(obj.basename(filename)); %zobrazim jmeno souboru s pouze koncem 
                    obj.filenames{fileno} = filename;
                    load(filename); %nacte vsechny promenne
                    test = ~P.data(:,P.sloupce.zpetnavazba); %index testovych epoch
                    
                    %identita jednotlivych epoch - epochData 
                    %predpokladam epochovana data, pro jina jsem ani nezkousel
                    obj.GetEpochData(epochData,fileno,test); %soucasne naplni obj.blokyprehazej, ktere muzu pouzivat dal                    
                    [d,tabs,RjEpochCh,P]=obj.PrehazejEpochy(d,tabs,RjEpochCh,P,test); %vyradi treningove trialy a prehazi epochy pokud je to treba
                    
                    %d hlavni data
                    obj.GetD(d); %ulozi nova data d, pripadne prehazi epochy
                    
                    %tabs a tabs_orig
                    obj.GetTabs(tabs,tabs_orig,fileno);                    
                    
                    %fs - vzorkovaci prekvence - musi byt pro vsechny soubory stejna
                    if ~isempty(obj.fs), assert(obj.fs == fs,'fs musi byt stejne');  else, obj.fs = fs; end
                    
                    %epochtime - cas epochy, napriklad -0.2 - 1.2, taky musi byt pro vsechny soubory stejne
                    if ~isempty(obj.epochtime)
                        assert(isequal(obj.epochtime,epochtime(1:2)),'epochtime musi byt stejne');
                    else
                        obj.epochtime = epochtime(1:2);
                    end
                    
                    %vyrazene epochy x kanaly
                    obj.RjEpochCh = cat(1,obj.RjEpochCh,RjEpochCh); %spojim pres kanaly, pocet epoch musi by stejny                    
                    
                    %PsychoPy data
                    obj.GetPsyData(P,fileno);                    
                    
                    %Hammer header
                    obj.GetHHeader(H,fileno);
                    
                    %frekvencni data
                    obj.GetHfreq(Hf,Hfmean,HFreq);
                    
                    %jen kopie - zatim nezpracovavam                                                           
                    obj.orig(fileno).DatumCas = DatumCas;
                    obj.orig(fileno).filename = filename;                                   
                end
                
            end
        end
        function [d,tabs,RjEpochCh,P]= PrehazejEpochy(obj,d,tabs,RjEpochCh,P,test) 
            %vyradi treningove epochy ze vsech dat 
            %pokud je potreba, prehazi epochy podle blokyprehazej
            d = d(:,:,test); %jen testove epochy, trening vymazu
            if ~isempty(obj.d) %pokud uz existuji data v obj.d - nejedna se o prvni soubor
                assert(size(obj.d,1)==size(d,1),['pocet vzorku musi byt stejny: d:' num2str(size(d,1)) ' x predchozi:' num2str(size(obj.d,1))]);
                assert(size(obj.d,3)==size(d,3),['pocet epoch musi byt stejny: d:' num2str(size(d,3)) ' x predchozi:' num2str(size(obj.d,3))]);
                assert(size(obj.d,3)==size(obj.epochData,1), 'pocet epoch v epochData a d musi byt stejny');
            end              
            
            tabs = tabs(:,test); %velikost tabs se odviji od velikosti d, takze nemusim overovat            
            RjEpochCh = RjEpochCh(:,test);
            Psy = CPsyData(P);%vytvorim objekt CPsyData 
            Psy.RemoveTraining(); %odstranim treninkove trialy
            P = Psy.P; %zase vytahnu strukturu P ven
            
            if ~isempty(obj.blokyprehazej) %pokud je treba prehazet epochy - jine poradi bloku nez v predchozim souboru
                %tam budu ukladat kopie dat s prehazenymi epochami
                d2=zeros(size(d)); 
                tabs2 = zeros(size(tabs));
                RjEpochCh2 = zeros(size(RjEpochCh));
                Pdata2 = zeros(size(P.data));
                epocha2=1; %nove cislo epochy 1-n
                for blok2 = 1:size(obj.blokyprehazej,1) %blok 16ti epoch (u AEdist) - nova cisla bloku 1-n
                    blok1 = obj.blokyprehazej(blok2,2); %puvodni cislo bloku, ktere chci zmenit na nove                
                    for epocha1 = obj.blokyprehazej(blok1,3) :  obj.blokyprehazej(blok1,4) %jednotlive epochy - puvodni cisla epoch                 
                       
                       %vlozim epochu na novou pozici
                       d2(:,:,epocha2)= d(:,:,epocha1); %time x channels x epochs
                       tabs2(:,epocha2) = tabs(:,epocha1); %time x epochs
                       RjEpochCh2(:,epocha2) = RjEpochCh(:,epocha1); %channels x epochs
                       Pdata2(epocha2,:) = P.data(epocha1,:);
                       
                       epocha2 = epocha2 + 1;
                    end
                end
                d = d2; %clear d2;
                tabs = tabs2; %clear tabs2;
                RjEpochCh = RjEpochCh2;
                P.data = Pdata2; %clear Pdata2;
            end
            
        end
        function GetD(obj,d)
            %ulozi nova eeg data k predchazejicim - prida je do spolecneho pole obj.d                     
            obj.d = cat(2,obj.d,d); %spojim pres channels - jen data z testovych epoch            
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
        end
        function GetHfreq(obj,Hf,Hfmean,HFreq)
            %spoji frekvencni data z predchozich a noveho souboru
            if isempty(obj.Hf)
                obj.Hf = Hf; %jen prvni soubor
            else
                assert(isequal(obj.Hf,Hf),'frekvencni pasma musi byt stejna');
            end
            if isempty(obj.Hfmean)
                obj.Hfmean = Hfmean; %jen prvni soubor
            end
            if isempty(obj.HFreq)
                obj.HFreq = HFreq;%time x channel x freq (x kategorie)
            else
                obj.HFreq = cat(2,obj.HFreq,HFreq);
            end
        end
        function obj = GetTabs(obj,tabs,tabs_orig,f)
            %spoji tabs z predchozich a noveho souboru
            if isempty(obj.tabs) %budu pouzivat tabs z prvniho souboru
                for ep = 1 : size(tabs,2) %pres vsechny epochy
                    if ep==1
                        tabs(:,1) = tabs(:,1) - tabs(1,1); %aby to zacinalo od 0
                    else
                        tabs(:,ep) = (tabs(:,ep) - tabs(1,ep)) + tabs(end,ep-1) + (tabs(end,ep-1)-tabs(end-1,ep-1));  %odectu prvni a pak prictu posledni z predchozi epochy a rozdil oproti predposlenimu
                    end
                end
                obj.tabs = tabs;                        
            end
            obj.orig(f).tabs = tabs;   %ale pro kazdy soubor ulozim originalni tabs a tabs_orig                     
            obj.orig(f).tabs_orig = tabs_orig; 
        end
        function obj = GetPsyData(obj,P,fileno)
            %pokud prvni soubor, ulozi P data
            %pokud dalsi soubor v poradi, ulozi P data z noveho souboru do zalohy v obj.ori            
            obj.orig(fileno).P = P; %psychopy data  
            if isempty(obj.PsyData)
                obj.PsyData = CPsyData(P); %vytvorim objekt CPsyData - z prvniho souboru
                obj.PsyData.RemoveTraining(); %odstranim treninkove trialy
            else
                obj.PsyData.P.pacientid = [ obj.PsyData.P.pacientid ',' P.pacientid]; %spojim pacient ID od vice pacientu
                %ostatni psydata zatim nepouzivam, reakcni uspesnost a rychlost bude u kazdeno jina
            end
        end
        function obj = GetEpochData(obj,epochData,fileno,test)
            %porovna poradi bloku v puvodnich souborej a novem souboru
            %pokud je jine, vytvori pole blokyprehazej, ktere se pak pouzije v PrehazejEpochy
            epochData =  epochData(test,:); %vyradim trening;
            if isempty(obj.epochData)                                                
                %bloky = obj.GetBlocks(epochData); 
                obj.epochData = epochData; %pouziju epochData prvniho souboru, ulozim bez treningu
                obj.blokyprehazej = []; %nebudu nic prehazovat
            else
                assert(size(obj.epochData,1)==size(epochData,1), 'pocet podminek v epochData musi byt stejny');                                                 
                bloky0 = obj.GetBlocks(obj.epochData); 
                bloky1 = obj.GetBlocks(epochData); 
                assert(size(bloky0,1)==size(bloky1,1),'pocet bloku musi byt stejny');
                prehazet = []; %v kterych blocich jsou jine podminky
                for b = 1:size(bloky0,1) %pres vsechny bloky
                   if bloky0(b,2) ~= bloky1(b,2)                               
                       assert(false,['ruzna velikost bloku ' b ]);
                   end
                   if bloky0(b,1) ~= bloky1(b,1), prehazet = [prehazet b]; end %#ok<AGROW> %budu muset poradi bloku prehazet, podminky v jinem poradi 
                end
                if numel(prehazet) > 0
                    obj.blokyprehazej = obj.ShuffleBlocks(bloky0, bloky1);                    
                else
                    obj.blokyprehazej = []; %nebudu nic prehazovat
                end
            end           
            obj.orig(fileno).epochData = epochData; % ulozim original taky                        
        end
        function obj = GetHHeader(obj,H,f)
            %spoji header (kanaly) v puvodnich souborech a tom novem
            if isempty(obj.CH)
                GetHHeader@CHilbert(obj,H);  
                obj.CH.chgroups = {1:numel(H.channels)};
                obj.CH.els = numel(H.channels);
                obj.els = obj.CH.els;
            else
                obj.CH.H.subjName = [obj.CH.H.subjName ',' H.subjName]; %spojim jmena subjektu
                obj.CH.H.channels = cat(2,obj.CH.H.channels,H.channels); %spojim udaje o kanalech
                obj.CH.chgroups{f} = (1:numel(H.channels)) + obj.CH.els(f-1);
                obj.CH.els(f) = numel(H.channels)+obj.CH.els(f-1);
                obj.orig(f).H = H; %orignalni header ulozim
            end
        end
        
    end
    methods (Static,Access = private)
        function bloky = GetBlocks(epochData)
            %ziska info o blocich 
            %data uz bez treningu
            typ = cell2mat(epochData(:,2)) ;
            b1 = [find(typ(2:end) ~= typ(1:end-1)); size(typ,1)]; %konce bloku
            b0 = [ 1; (b1(1:end-1)+1)]; %zacatky bloku
            bloky = [typ(b0),b1-b0+1,b0,b1];
        end
        function blokyprehazej = ShuffleBlocks(bloky0, bloky1)
            %zjisti indexy k prehazeni bloku, v bloky1=Cil podle bloky0=Zdroj
            blokyprehazej = [bloky0(:,1) zeros(size(bloky0,1),1) bloky0(:,3:4)]; %vzorove poradi bloku , k tomu budu pridavat cisla 
            lastkat = [unique(bloky1(:,1)) zeros(numel(unique(bloky1(:,1))),1)]; %tam si budu ukladat podledni nalezena cisla bloku
            for b = 1:size(blokyprehazej,1)
                
                blok = blokyprehazej(b,1); %cislo podminky z bloku 1 v Zdroj
                ilastkat = find(lastkat(:,1)==blok,1); %index podminky v kastkat
                temp = find(bloky1(lastkat(ilastkat,2)+1:end,1)==blok,1)+lastkat(ilastkat,2); %najde dalsi blok s touto podminkou v Cil
                %find hleda jen v te casti bloky1, takze musim pricist index, odkud jsem hledal
                lastkat(ilastkat,2) = temp;               
                blokyprehazej(b,2)=temp;
            end
            
        end       
        function [str]= basename(filename)
            % vraci filename s koncem path pro identifikaci pacienta
            %[path,basename,ext] = fileparts(filename); %takhle to nechci
            fslash = strfind(filename,'\');
            str = filename(fslash(end-2)+1:end); %dve casti path pred basename
        end
        
    end

    
end

