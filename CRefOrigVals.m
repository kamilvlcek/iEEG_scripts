classdef CRefOrigVals < matlab.mixin.Copyable 
    %CREFORIGVALS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        OrigVals = []; %matrix kategories x channels x 2 - maximalni hodnoty kanaly
        kats; %kategorie, ktere patri k OrigVals
        Eh; %odkaz na objekt CHilbertMulti
        setup;  %nastaveni podle testu
    end
    
    methods (Access = public)
        function obj = CRefOrigVals(E)
            obj.Eh = E;                
            obj.setup = eval(['setup_' E.PsyData.testname]); %nactu nastaveni
        end
        function GetData(obj)
            obj.OrigVals = zeros(numel(obj.Eh.Wp(obj.Eh.WpActive).kats),numel(obj.Eh.CH.H.channels),2);
            obj.kats = obj.Eh.Wp(obj.Eh.WpActive).kats;
            lastiFile = 0;
            for ch = 1:numel(obj.Eh.CH.H.channels) %pres vsechny bipolarni kanaly
                iFile = find(obj.Eh.els >= ch,1); %cislo souboru, ze ktereho nactu data pro aktualni kanal
                [ch1name,ch2name]= obj.extractChannelNames(obj.Eh.CH.H.channels(ch).name);
                if iFile ~= lastiFile
                    [pacientId,fnameOrig] = obj.extractPacient(obj.Eh.filenames{iFile});
                    E = pacient_load(pacientId,obj.Eh.PsyData.testname,fnameOrig);
                    valmax = zeros(E.channels,numel(obj.kats));
                    for k = 1:numel(obj.kats) %nactu si z tohohle souboru maxima pro vsechny kategorie
                        valmax(:,k) = E.ResponseTriggerTime(0.9,0.9,obj.kats(k));  
                    end
                    lastiFile = iFile; %posledni nactene cislo souboru
                end
                if ~isempty(E)
                    ch1 = E.CH.ChannelNameFind(ch1name);
                    ch2 = E.CH.ChannelNameFind(ch2name);
                    if ch1 < 1 || ch2 < 1
                        warning([chn2name 'or' chn2name 'not found in ' pacientId]);
                        continue;
                    end
                    for k = 1:numel(obj.kats) %index kategorie v aktualni statistice, zatim budu pracovat jen s jednoduchymi kontrasty, bez dvojic                        
                        obj.OrigVals(k,ch,1) = valmax(ch1,k);
                        obj.OrigVals(k,ch,2) = valmax(ch2,k);
                    end
                end
                
            end
        end
    end
    methods (Access = private)
        function [pacientId,fnameOrig] = extractPacient(obj,filename)
            %extrahuje ze jmena extraktu id pacienta a jmeno souboru s originalni referenci
            [pathstr,fname] = fileparts(filename);  
            path2= strrep(pathstr,obj.setup.basedir,'');
            pacientId = char(extractBefore(path2,'\')); %id pacienta, napriklad "p073 Pech VT6"
            mezery = strfind(fname,' ');
            fnamepart = extractBefore(fname,mezery(end)); %napriklad "PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2018-08"
            fnameOrig = char(strrep(fnamepart,'refBipo','refOrig'));
        end
    end
    methods (Static,Access = private)
        function [chn1,chn2]= extractChannelNames(biponame)
            
            %pro P82 (R-INS9-R-INS11) tohle nefunguje, a nevim jak to udelat pomoci regexp
%             [~,tok] = regexp(biponame,'\(([\w\d-]+)-([\w\d-]+)\)','match','tokens');
%             chn1 = tok{1}{1};
%             chn2 = tok{1}{2};
              n = char(extractBetween(biponame,'(',')'));
              pos = strfind(n,'-');
              if numel(pos) == 1 %A1-A2
                  chn1 = n(1:pos-1);
                  chn2 = n(pos+1:end);
              elseif numel(pos) == 3 %R-INS9-R-INS11
                  chn1 = n(1:pos(2)-1);
                  chn2 = n(pos(2)+1:end);
              else % teoreticky R-INS9-A2, pos = 
                  if pos(1) < diff(pos) %pokud je prvni jmeno delsi R-INS9-A2
                      chn1 = n(1:pos(2)-1);
                      chn2 = n(pos(2)+1:end);
                  else
                      chn1 = n(1:pos(2)-1);
                      chn2 = n(pos(2)+1:end);
                  end
              end
        end
    end
    
end

