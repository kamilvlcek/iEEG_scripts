classdef CHilbertMulti < CHilbert 
    %CHILBERTMULTI zpracovava data z nekolika CHilbert vyexportovanych data
    %   nacita data vyprodukovana pomoci CHilbert.ExtractData
    
    properties
        orig; %struktura s puvodnimi daty 
        filenames; %cell array se jmeny souboru
    end
    
    methods
        function obj = ImportExtract(filenames)
            for f = 1:numel(filenames)
                filename = filenames{f}; %cell array, zatim to musi byt plna cesta                
                if exist(CHilbert.filenameH(filename),'file')  
                    obj.filenames{f} = filename;
                    load(CHilbert.filenameH(filename));
                    %d hlavni data
                    if ~isempty(obj.d)
                        assert(size(obj.d,1)==size(d,1),['pocet vzorku musi byt stejny: d:' num2str(size(d,1)) ' x predchozi:' num2str(size(obj.d,1))]);
                        assert(size(obj.d,3)==size(d,3),['pocet epoch musi byt stejny: d:' num2str(size(d,1)) ' x predchozi:' num2str(size(obj.d,1))]);
                    end
                    obj.d = cat(2,obj.d,d); %spojim pres channels
                    %tabs a tabs_orig
                    if isempty(obj.tabs)
                        for ep = 1 : size(tabs,2) %pres vsechny epochy
                            if ep==1
                                tabs(:,1) = tabs(:,1) - tabs(1,1); %aby to zacinalo od 0
                            else
                                tabs(:,ep) = (tabs(:,ep) - tabs(1,ep)) + tabs(end,ep-1) + (tabs(end,ep-1)-tabs(end-1,ep-1));  %odectu prvni a pak prictu posledni z predchozi epochy a rozdil oproti predposlenimu
                            end
                        end
                        obj.tabs = tabs;                        
                    end
                    obj.orig(f).tabs = tabs;                        
                    obj.orig(f).tabs_orig = tabs_orig;                        
                   %fs
                    if ~isempty(obj.fs), assert(obj.fs == fs,'fs musi byt stejne');  else, obj.fs = fs; end
                    
                    %epochtime
                    if ~isempty(obj.epochtime)
                        assert(isequal(obj.epochtime,epochtime(1:2)),'epochtime musi byt stejne');
                    else
                        obj.epochtime = epochtime(1:2);
                    end
                    
                    %vyrazene epochy x kanaly
                    obj.RjEpochCh = cat(1,obj.RjEpochCh,RjEpochCh); %spojim pres kanaly, pocet epoch musi by stejny
                    
                    %jen kopie                                        
                    obj.orig(f).P = P; %psychopy data
                    obj.orig(f).epochData = epochData; %identita jednotlivych epoch
                    obj.orig(f).DatumCas = DatumCas; 
                    %TODO porovnat poradi podminek
                    
                    %todo H
                    
                    %budu mit spoustu chybejicich promennych
                    %nejake funkce mi urcite nebudou fungovat 
                    %budu doplnovat podle prvniho souboru v poradi?
                    
                    
                end
            end
        end
    end
    
end

