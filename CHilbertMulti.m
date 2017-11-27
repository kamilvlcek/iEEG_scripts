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
                    load(filename);
                    %d hlavni data
                    if ~isempty(obj.d)
                        assert(size(obj.d,1)==size(d,1),['pocet vzorku musi byt stejny: d:' num2str(size(d,1)) ' x predchozi:' num2str(size(obj.d,1))]);
                        assert(size(obj.d,3)==size(d,3),['pocet epoch musi byt stejny: d:' num2str(size(d,1)) ' x predchozi:' num2str(size(obj.d,1))]);
                    end
                    obj.d = cat(2,obj.d,d); %spojim pres channels
                    [obj.samples,obj.channels, obj.epochs] = obj.DSize();
                    %tabs a tabs_orig
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
                    
                    %jen kopie - zatim nezpracovavam                                                           
                    obj.orig(f).DatumCas = DatumCas;
                    
                    %budu mit spoustu chybejicich promennych
                    %nejake funkce mi urcite nebudou fungovat 
                    %budu doplnovat podle prvniho souboru v poradi?                    
                    
                end
                
            end
        end
        function obj = GetPsyData(obj,P,f)
            obj.orig(f).P = P; %psychopy data  
            obj.PsyData = CPsyData(P); %vytvorim objekt CPsyData - z prvniho souboru
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
                obj.H.subjName = [obj.H.subjName ',' H.subjName]; %spojim jmena subjektu
                obj.H.channels = cat(1,obj.H.channels,H.channels); %spojim udaje o kanalech
                obj.orig(f).H = H; %orignalni header ulozim
            end
        end
    end
    
end

