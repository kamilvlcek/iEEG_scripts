classdef CHilbertMulti < CHilbert 
    %CHILBERTMULTI zpracovava data z nekolika CHilbert vyexportovanych data
    %   nacita data vyprodukovana pomoci CHilbert.ExtractData
    
    properties
        orig; %struktura s puvodnimi daty 
        filenames; %cell array se jmeny souboru
    end
    
    methods
        function obj = CHilbertMulti(filenames)  
            obj.ImportExtract(filenames);
        end
        function obj = ImportExtract(obj,filenames)
            for f = 1:numel(filenames)
                filename = filenames{f}; %cell array, zatim to musi byt plna cesta                
                if exist(filename,'file')  
                    obj.filenames{f} = filename;
                    load(filename); %nacte vsechny promenne
                    %d hlavni data
                    if ~isempty(obj.d)
                        assert(size(obj.d,1)==size(d,1),['pocet vzorku musi byt stejny: d:' num2str(size(d,1)) ' x predchozi:' num2str(size(obj.d,1))]);
                        assert(size(obj.d,3)==size(d,3),['pocet epoch musi byt stejny: d:' num2str(size(d,3)) ' x predchozi:' num2str(size(obj.d,3))]);
                    end
                    obj.d = cat(2,obj.d,d); %spojim pres channels
                    [obj.samples,obj.channels, obj.epochs] = obj.DSize();
                    %tabs a tabs_orig
                    obj.GetTabs(tabs,tabs_orig,f);   
                    
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
                    obj.GetPsyData(P,f);
                    
                    %identita jednotlivych epoch - epochData 
                    obj.GetEpochData(epochData,f);
                    %Hammer header
                    obj.GetHHeader(H,f);
                    
                    %frekvencni data
                    obj.GetHfreq(Hf,Hfmean,HFreq);
                    
                    %jen kopie - zatim nezpracovavam                                                           
                    obj.orig(f).DatumCas = DatumCas;
                    obj.orig(f).filename = filename;
                    %budu mit spoustu chybejicich promennych
                    %nejake funkce mi urcite nebudou fungovat 
                    %budu doplnovat podle prvniho souboru v poradi?                    
                    
                end
                
            end
        end
        function GetHfreq(obj,Hf,Hfmean,HFreq)
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
        function obj = GetPsyData(obj,P,f)
            obj.orig(f).P = P; %psychopy data  
            if isempty(obj.PsyData)
                obj.PsyData = CPsyData(P); %vytvorim objekt CPsyData - z prvniho souboru
            else
                obj.PsyData.P.pacientid = [ obj.PsyData.P.pacientid ',' P.pacientid]; %spojim pacient ID od vice pacientu
            end
        end
        function obj = GetEpochData(obj,epochData,f)
            if isempty(obj.epochData)
                assert(size(epochData,1)==size(obj.d,3), 'pocet podminek v epochData a epoch musi byt stejny');
                obj.epochData = epochData; %pouziju epochData prvniho souboru 
            else
                assert(size(obj.epochData,1)==size(epochData,1), 'pocet podminek v epochData musi byt stejny');                                                 
                for ep = 1:size(epochData,1) %pres vsechny podminky
                   if strcmp(obj.epochData{ep,1},epochData{ep,1})==0                               
                       assert(false,['ruzne podminky pro epochu ' ep ]);
                   end
                end
            end
            obj.orig(f).epochData = epochData; % ulozim original taky                        
        end
        function obj = GetHHeader(obj,H,f)
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
    
end

