classdef CRefOrigVals < matlab.mixin.Copyable 
    %CREFORIGVALS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        ValMax = []; %matrix kategories x channels x 2 - maximalni hodnoty kanaly
        TMax = []; %matrix kategories x channels x 2 - cas maxima
        kats; %kategorie, ktere patri k OrigVals
        Eh; %odkaz na objekt CHilbertMulti
        setup;  %nastaveni podle testu
        chnames = {}; %nazvy originalnich kanalu z kterych je to spocitane
        chnums = []; %cisla originalnich kanalu, kvuli kontrole
    end
    
    methods (Access = public)
        function obj = CRefOrigVals(E)
            obj.Eh = E;                
            obj.setup = eval(['setup_' E.PsyData.testname]); %nactu nastaveni
        end
        function GetData(obj)
            obj.ValMax = zeros(numel(obj.Eh.Wp(obj.Eh.WpActive).kats),numel(obj.Eh.CH.H.channels),2);
            obj.TMax = zeros(numel(obj.Eh.Wp(obj.Eh.WpActive).kats),numel(obj.Eh.CH.H.channels),2);
            obj.kats = obj.Eh.Wp(obj.Eh.WpActive).kats;
            obj.chnames = cell(numel(obj.Eh.CH.H.channels),6);
            obj.chnums = zeros(numel(obj.Eh.CH.H.channels),2);
            lastiFile = 0;
            for ch = 1:numel(obj.Eh.CH.H.channels) %pres vsechny bipolarni kanaly
                iFile = find(obj.Eh.els >= ch,1); %cislo souboru, ze ktereho nactu data pro aktualni kanal
                [ch1name,ch2name]= obj.extractChannelNames(obj.Eh.CH.H.channels(ch).name);                
                if iFile ~= lastiFile
                    [pacientId,fnameOrig] = obj.extractPacient(obj.Eh.filenames{iFile});
                    E = pacient_load(pacientId,obj.Eh.PsyData.testname,fnameOrig);
                    valmax = zeros(E.channels,numel(obj.kats));
                    tmax = zeros(E.channels,numel(obj.kats));
                    for k = 1:numel(obj.kats) %nactu si z tohohle souboru maxima pro vsechny kategorie
                        [valmax(:,k),tmax(:,k)] = E.ResponseTriggerTime(0.9,0.9,obj.kats(k));  
                    end
                    lastiFile = iFile; %posledni nactene cislo souboru
                end
                obj.chnames(ch,:) = {ch1name,ch2name,obj.Eh.CH.H.channels(ch).name,pacientId,fnameOrig,~isempty(E)};
                if ~isempty(E)
                    ch1 = E.CH.ChannelNameFind(ch1name);
                    ch2 = E.CH.ChannelNameFind(ch2name);
                    obj.chnums(ch,:) = [ch1,ch2];
                    if ch1 < 1 || ch2 < 1
                        warning([chn2name 'or' chn2name 'not found in ' pacientId]);
                        continue;
                    end
                    for k = 1:numel(obj.kats) %index kategorie v aktualni statistice, zatim budu pracovat jen s jednoduchymi kontrasty, bez dvojic                        
                        obj.ValMax(k,ch,1) = valmax(ch1,k);
                        obj.ValMax(k,ch,2) = valmax(ch2,k);
                        obj.TMax(k,ch,1) = tmax(ch1,k);
                        obj.TMax(k,ch,2) = tmax(ch2,k);
                    end
                end
                
            end
        end
        function PlotCh(obj,ch)
            figure('Name',['CRefOrigVals - ch' num2str(ch)]);
            hold on;
            baseColors = [0 1 0;  1 0 0; 0 0 1; 1 1 0; 1 0 1; 0 1 0];
            katnames = cell(1,numel(obj.kats));
            for k = 1:numel(obj.kats)
                katnames{k} = obj.Eh.PsyData.CategoryName(obj.kats(k));
                vals = squeeze(obj.ValMax(k,ch,:));
                times = squeeze(obj.TMax(k,ch,:));
                plot(times,vals,'-o','Color',baseColors(obj.kats(k),:));
                text(times(1),vals(1),obj.chnames{ch,1});
                text(times(2),vals(2),obj.chnames{ch,2});
            end
            xlim(obj.Eh.epochtime(1:2));
            legend(katnames);
            title(['channel ' num2str(ch) ', ' obj.chnames{ch,1} '-' obj.chnames{ch,2}]);                 
            hold off;
        end
        function E = PlotResponseCh(obj,ch)
            E = pacient_load(obj.chnames{ch,4},obj.Eh.PsyData.testname,obj.chnames{ch,5});
            if ~isempty(E)
                E.PlotResponseCh(obj.chnums(ch,1));
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

