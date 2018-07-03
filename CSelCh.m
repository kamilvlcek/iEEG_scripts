classdef CSelCh < matlab.mixin.Copyable 
    %CSelCh Trida uchovavajici a poskutujici seznamy vybranych kanalu pro soubory CHIlbertMulti aj
    %   pred frekvence, kontrasty aj
    
    properties (Access = public)
        selCh;
        n; %pocet ulozenych dat
        filename;
    end
    
    methods (Access = public)
        function obj = CSelCh(filename)
            if exist('filename','var')                
                obj.Load(filename);            
            else
                obj.selCh = cell(1,5); %prvni sloupec je filename, druhy selCh, treti katstr, ctvrty freq , paty chnames
                obj.n = 0;
            end
        end
        function obj = SetSelCh(obj,selCh,filename,chnames, katstr,freq)            
            %ulozi vyber kanalu do teto tridy            
            if ~exist('chnames','var'), chnames = []; end
            if ~exist('katstr','var'), katstr = []; end
            if ~exist('freq','var'), freq = []; end
            if isa(selCh,'CiEEGData')
                %muzu predat jako parametr celou tridu CM aj, data z nim ziskam sam
                E = selCh;
                selCh = E.GetSelCh();
                filename = E.filename;
                if isprop(E,'label')
                    katstr = E.label; 
                else
                    katstr = []; 
                end
                if isprop(E, 'Hf') 
                    freq = [num2str(E.Hf(1)) '-' num2str(E.Hf(end))];
                else
                    freq = [];
                end
                chnames = {E.CH.H.channels(selCh).name}; %jmena vybranych kanaly
                obj.SetSelCh(selCh,filename,chnames,katstr,freq);
            else
                if obj.n>0
                    s = find(~cellfun(@isempty,strfind(obj.selCh(:,1),filename)),1); %najdu pouze prvni vyhovujici soubor
                else
                    s = [];
                end
                if isempty(s) %filename neexistuje, ulozim
                    obj.selCh(obj.n+1,:) = {filename,selCh,katstr,freq,chnames};
                    obj.n = obj.n +1;
                else %filename uz ulozen, radku s nim prepisu
                    obj.selCh(s,:) = {filename,selCh,katstr,freq,chnames}; 
                end
                disp([ num2str(numel(selCh)) ' selected channels saved']);                   
            end
        end
        function selCh = GetSelCh(obj,filename)
            %ziska vyber kanalu z teto tridy
            if isa(filename,'CiEEGData') %pokud predam jako parametr tridu, ulozi data primo do ni
                 E = filename;
                 s = find(~cellfun(@isempty,strfind(obj.selCh(:,1),E.filename)),1); %najdu pouze prvni vyhovujici soubor 
                 if ~isempty(s)
                     selCh = obj.selCh{s,2};                     
                     selChNames = obj.selCh{s,5};
                     for ich = 1:numel(selCh)
                        if ~isempty(obj.selCh{s,5})
                            assert(strcmp(selChNames{ich},E.CH.H.channels(selCh(ich)).name),['channel ' num2str(selCh(ich)) ' nesedi s popisem: ' selChNames{ich} ' vs.' E.CH.H.channels(selCh(ich)).name]);                                                        
                        end
                     end
                     E.plotRCh.selCh = selCh; %vlozim cisla kanalu do objektu, tam kam patri
                     if isempty(E.label) && ~isempty(obj.selCh{s,3})
                         E.label = obj.selCh{s,3}; %pokud je v objektu prazne label a tato trida ho obsahuje, vyplnim ho taky
                         disp(['nastaveno label: ' E.label]);
                     end
                     disp('ulozeno do vlastnosti tridy');
                 else
                     selCh = [];
                     disp('soubor nenalezen v obj.selCh');
                 end                     
            else
                s = find(~cellfun(@isempty,strfind(obj.selCh(:,1),filename)),1); %najdu pouze prvni vyhovujici soubor
                if ~isempty(s)
                    selCh = obj.selCh{s,2};
                else
                    disp('soubor nenalezen v obj.selCh');   
                    selCh = [];
                end
            end
        end
        function Save(obj,filename)
            if ~exist('filename','var') 
                filename = obj.filename;
            else
                obj.filename = filename;
            end
            selCh = obj.selCh; %#ok<NASGU,PROPLC>
            n = obj.n;      %#ok<NASGU,PROPLC>                  
            save(filename,'selCh','n','filename','-v7.3');  
            disp(['ulozeno do ' filename]); 
        end
        function obj = Load(obj,filename)            
            assert(exist(filename,'file')==2, 'soubor neexistuje');
            load(filename,'selCh','n','filename');              
            obj.n = n; %#ok<CPROPLC>
            obj.selCh = selCh; %#ok<CPROPLC>
            obj.filename = filename;
            disp(['nacten soubor ' filename]); 
        end
            
    end
    
end

