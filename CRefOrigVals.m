classdef CRefOrigVals < matlab.mixin.Copyable 
    %CREFORIGVALS trida shromazdujici udaje o odpovedich v originalni referenci
    %   Detailed explanation goes here
    
    properties (Access = public)
        ValMax = []; %matrix kategories x channels x 2 - maximalni hodnoty kanaly
        TMax = []; %matrix kategories x channels x 2 - cas maxima
        kats; %kategorie, ktere patri k OrigVals
        Eh; %odkaz na objekt CHilbertMulti
        setup;  %nastaveni podle testu
        chnames = {}; %nazvy originalnich kanalu z kterych je to spocitane
        neuroLabels = {}; %jmena neurologyLabels 
        chnums = []; %cisla originalnich kanalu, kvuli kontrole
        PlotChH; %handle na obrazek
        PlotChCh; %aktualni zobrazeny kanal v PlotCh
        PlotRCh; %udaje o PlotResponseCh 
        chsignif; %channel significance relative to baseline
    end
    
    methods (Access = public)
        function obj = CRefOrigVals(E)
            if ~ strcmp(E.reference,'Bipolar')
                disp('CRefOrigVals: the original CiEEGData data should be bipolar');
            elseif isempty(E.PsyData)
                disp('CRefOrigVals: no PsyData');
            else
                obj.Eh = E;  
                obj.setup = eval(['setup_' E.PsyData.testname]); %nactu nastaveni
                obj.Load();
            end            
        end
        function GetData(obj)
            if isempty(obj.Eh)
                disp('no original CiEEGData data loaded');
                return; 
            end
            obj.ValMax = zeros(numel(obj.Eh.Wp(obj.Eh.WpActive).kats),numel(obj.Eh.CH.H.channels),2);
            obj.TMax = zeros(numel(obj.Eh.Wp(obj.Eh.WpActive).kats),numel(obj.Eh.CH.H.channels),2);
            obj.kats = obj.Eh.Wp(obj.Eh.WpActive).kats;
            obj.chnames = cell(numel(obj.Eh.CH.H.channels),6);
            obj.chnums = zeros(numel(obj.Eh.CH.H.channels),2);
            obj.neuroLabels = cell(numel(obj.Eh.CH.H.channels),2); %neurologyLabels pro ty dva originalni kanaly
            obj.chsignif = false(numel(obj.Eh.CH.H.channels),numel(obj.kats),3); % channels x kats x [bipolar and both unipolar signals]
            if isfield(obj.Eh.plotRCh,'selChSignum') && ~isempty(obj.Eh.plotRCh.selChSignum)
                signum = obj.Eh.plotRCh.selChSignum;  %jestli se maji pocitat jen vychylky nahoru, nebo dolu, nebo oboji
            else
                signum = [];
            end
            prumery = obj.Eh.IntervalyResp([],[],signum,0); %to get significance relative to baseline, kats in prumery are ordered according to current SetStatActive            
            lastiFile = 0;
            for ch = 1:numel(obj.Eh.CH.H.channels) %pres vsechny bipolarni kanaly
                iFile = find(obj.Eh.els >= ch,1); %cislo souboru, ze ktereho nactu data pro aktualni kanal
                [ch1name,ch2name]= obj.extractChannelNames(obj.Eh.CH.H.channels(ch).name);                
                if iFile ~= lastiFile %for the first channel of new patient
                    [pacientId,fnameOrig] = obj.extractPacient(obj.Eh.filenames{iFile});
                    E = pacient_load(pacientId,obj.Eh.PsyData.testname,fnameOrig);
                    valmax = zeros(E.channels,numel(obj.kats));
                    tmax = zeros(E.channels,numel(obj.kats));                    
                    for k = 1:numel(obj.kats) %nactu si z tohohle souboru maxima pro vsechny kategorie                       
                        [valmax(:,k),tmax(:,k)] = E.ResponseTriggerTime(0.9,0.9,obj.kats(k),[],signum);                          
                    end
                    lastiFile = iFile; %posledni nactene cislo souboru
                end
                obj.chnames(ch,:) = {ch1name,ch2name,obj.Eh.CH.H.channels(ch).name,pacientId,fnameOrig,~isempty(E)};
                if ~isempty(E) %for each channel
                    ch1 = E.CH.ChannelNameFind(ch1name);
                    ch2 = E.CH.ChannelNameFind(ch2name);
                    obj.chnums(ch,:) = [ch1,ch2];
                    if ch1 < 1 || ch2 < 1
                        warning([ch1name ' or ' ch2name ' not found in ' pacientId ' for channel ' num2str(ch) ]);
                        continue;
                    end
                    prumeryE = E.IntervalyResp([],[],signum,0);
                    for k = 1:numel(obj.kats) %index kategorie v aktualni statistice, zatim budu pracovat jen s jednoduchymi kontrasty, bez dvojic                        
                        obj.ValMax(k,ch,1) = valmax(ch1,k);
                        obj.ValMax(k,ch,2) = valmax(ch2,k);
                        obj.TMax(k,ch,1) = tmax(ch1,k);
                        obj.TMax(k,ch,2) = tmax(ch2,k);
                        obj.chsignif(ch,k,1)=prumery(ch,1,k)>0; %is the reponses in this kat and this ch signif relative to baseline
                        obj.chsignif(ch,k,2)=prumeryE(ch1,1,k)>0;
                        obj.chsignif(ch,k,3)=prumeryE(ch2,1,k)>0;
                    end
                    obj.neuroLabels(ch,:)={E.CH.H.channels(ch1).neurologyLabel, E.CH.H.channels(ch2).neurologyLabel};
                end
                
            end
        end
        function PlotCh(obj,ch)
            if isempty(obj.Eh)
                disp('no original CiEEGData data loaded');
                return; 
            end
            assert(~isempty(obj.ValMax),'CRefOrigVals: nejsou nactena data');
            assignhandle = false;
            if ~exist('ch','var'), ch = 1; end
            if isempty(obj.PlotChH) % Pokud je obj.PlotChH prazdne pole, ishghandle() vraci prazdne logicke pole, ktere se interpertuje jako true a nejde porovnat se skalarni logickou hodnotou isempty() (R2018a linux)
                obj.PlotChH = figure('Name',['CRefOrigVals - ch' num2str(ch)]);
                assignhandle = true;
            elseif ~ishghandle(obj.PlotChH)
                obj.PlotChH = figure('Name',['CRefOrigVals - ch' num2str(ch)]);
                assignhandle = true;
            else
                figure(obj.PlotChH);
                clf; %vymazu obrazek                 
            end
            hold on;
            baseColors = [0 1 0;  1 0 0; 0 0 1; 1 1 0; 1 0 1; 0 1 0];
            katnames = cell(1,numel(obj.kats));
            for k = 1:numel(obj.kats)
                katnames{k} = obj.Eh.PsyData.CategoryName(obj.kats(k));
                vals = squeeze(obj.ValMax(k,ch,:));
                times = squeeze(obj.TMax(k,ch,:));
                katcolor = baseColors(obj.kats(k)+iff(obj.kats(1)==0,1,0),:);
                plot(times,vals,'-o','Color',katcolor);                                
                text(times(1),vals(1),[obj.chnames{ch,1} ',' obj.neuroLabels{ch,1} ',' iff(obj.chsignif(ch,k,2),'*','')]);
                text(times(2),vals(2),[obj.chnames{ch,2}  ',' obj.neuroLabels{ch,2} ',' iff(obj.chsignif(ch,k,3),'*','')]);
            end
            xlim(obj.Eh.epochtime(1:2));
            legend(katnames);
            title(['channel ' num2str(ch) ', ' obj.chnames{ch,3}]);                 
            hold off;
            obj.PlotChCh = ch;
            if assignhandle
                methodhandle = @obj.hybejPlotCh;
                set(obj.PlotChH,'KeyPressFcn',methodhandle);
            end
        end
        function E = PlotResponseCh(obj,ch)
            if isempty(obj.Eh)
                disp('no original CiEEGData data loaded');
                return; 
            end
            plotfigure = true;
            if isfield(obj.PlotRCh,'fh') && ~isempty(obj.PlotRCh.fh) && ishghandle(obj.PlotRCh.fh) 
                if obj.PlotRCh.ch == ch  % Pokud je obj.PlotChH prazdne pole, ishghandle() vraci prazdne logicke pole, ktere se interpertuje jako true a nejde porovnat se skalarni logickou hodnotou isempty() (R2018a linux)
                    figure(obj.PlotRCh.fh);                                
                    plotfigure = false;
                end                
            end
            if plotfigure
                assert(~isempty(obj.ValMax),'CRefOrigVals: nejsou nactena data');
                nick = obj.chnames{ch,4};
                if  isfield(obj.PlotRCh,'nick') && strcmp(nick,obj.PlotRCh.nick) %pokud je nick pacienta stejny jako ulozeny
                    E = obj.PlotRCh.Eh;                    
                else
                    E = pacient_load(nick,obj.Eh.PsyData.testname,obj.chnames{ch,5});
                    obj.PlotRCh.Eh = E;
                    obj.PlotRCh.nick = nick;   
                    if isfield(obj.PlotRCh,'fh') && ~isempty(obj.PlotRCh.fh) && ishghandle(obj.PlotRCh.fh) 
                        fPosition = obj.PlotRCh.fh.Position; %ulozim si pozici predchoziho obrazku
                        close(obj.PlotRCh.fh);
                    end
                end
                if ~isempty(E)
                    E.PlotResponseCh(obj.chnums(ch,1));                    
                end
                obj.PlotRCh.ch = ch;
                obj.PlotRCh.fh = E.plotRCh.fh;         
                if exist('fPosition','var'), obj.PlotRCh.fh.Position = fPosition; end
            end           
        end
            
        function Save(obj)
            if isempty(obj.Eh)
                disp('no original CiEEGData data loaded');
                return; 
            end
            fname = saveFileName(obj.Eh.filename, 'CRefOrig');
            ValMax = obj.ValMax; %#ok<PROP,NASGU>
            TMax = obj.TMax; %#ok<PROP,NASGU>
            kats = obj.kats; %#ok<PROP,NASGU>
            setup = obj.setup; %#ok<PROP,NASGU>
            chnames = obj.chnames; %#ok<PROP,NASGU>
            chnums = obj.chnums; %#ok<PROP,NASGU>
            neuroLabels = obj.neuroLabels; %#ok<PROP,NASGU>
            chsignif = obj.chsignif; %#ok<PROP,NASGU>
            save(fname,'ValMax','TMax','kats','setup','chnames','chnums','neuroLabels','chsignif','-v7.3'); %do noveho souboru data z teto tridy
            disp(['saved to ' fname ]);
        end
        function Load(obj)
            if isempty(obj.Eh)
                disp('no original CiEEGData data loaded');
                return; 
            end
            fname = saveFileName(obj.Eh.filename, 'CRefOrig');
            if exist(fname,'file')
                V = load(fname);
                obj.ValMax = V.ValMax;  
                obj.TMax = V.TMax; 
                obj.kats = V.kats; 
                obj.setup = V.setup; 
                obj.chnames = V.chnames; 
                obj.chnums = V.chnums; 
                obj.neuroLabels = V.neuroLabels;
                if isfield(V,'chsignif')
                    obj.chsignif = V.chsignif;
                else
                    warning('the chsignif not found in the saved file');
                end
                disp(['loaded ' fname ]);
            else
                disp(['not found: ' fname ]);
            end
        end
        function ie = isEmpty(obj)
            %ISEMPTY returns true if the object does not contain any values
            ie = isempty(obj.ValMax);
        end
        function ExportXLS(obj)
            if isempty(obj.Eh)
                disp('no original CiEEGData data loaded');
                return; 
            end
            varsfirst = {'chnum' 'bipolarName','maxNeurologyLabel'};
            varskat = {'_valmax1' '_valmax2' '_valmaxdiff' '_signif' '_maxChName' '_maxNeurologyLabel' };
            cellOut = cell(size(obj.chnums,1), numel(varsfirst) + numel(varskat)*numel(obj.kats));
            variableNames = cell(1,numel(varsfirst)+5*numel(obj.kats));
            variableNames(1:numel(varsfirst)) = varsfirst;
            for k = 1:numel(obj.kats)   % projdu vsechny kategorie a nastavim headery pro prislusne sloupce
                katname = obj.Eh.PsyData.CategoryName(obj.kats(k));
                katname_vars = repmat({katname},1,numel(varskat));
                for kv = 1:numel(katname_vars)
                    katname_vars{kv} = [ katname_vars{kv} varskat{kv} ];
                end
                variableNames((k-1)*numel(varskat) + numel(varsfirst) + (1:numel(varskat))) = katname_vars;
            end
                        
            for ch = 1:size(obj.chnums,1)    % projdu vsechny kanaly
                lineOut = cell(1,1+5*numel(obj.kats));
                lineOut(1:2) = {ch , obj.chnames{ch,3}};
                neuroLabelsKats = cell(numel(obj.kats),2); %hodnoty pro kategorie zvlast
                for k = 1:numel(obj.kats)   % projdu vsechny kategorie (ruzne kategorie pro stejny kanal budou v tabulce pod sebou)
                    vals = squeeze(obj.ValMax(k,ch,:));
                    if vals(1) > vals(2)
                        maxId = 1;
                    else
                        maxId = 2;
                    end
                    maxChName = obj.chnames{ch,maxId};                   
                    maxChNeurologyLabel = obj.neuroLabels{ch,maxId};                   
                    neuroLabelsKats(k,:) = {maxChNeurologyLabel ,  abs(vals(1)-vals(2))} ; %for this kat, neuroLabel of largest respons and abs diff
                    lineOut((k-1)*numel(varskat)+ numel(varsfirst) + (1:numel(varskat))) = {vals(1), vals(2), abs(vals(1)-vals(2)), ...
                        num2str(double(squeeze(obj.chsignif(ch,k,:)))','%d_%d_%d'), maxChName, maxChNeurologyLabel};
                end
                if numel(unique(neuroLabelsKats(:,1)))==1 
                    %if for all kats (with signif bipolar response) the response is larger on the same origRef channel
                    lineOut(3) = neuroLabelsKats(1,1);
                elseif numel(unique(neuroLabelsKats(obj.chsignif(ch,:,1),1)))==1
                    %if for all kats (with signif bipolar response) the response is larger on the same origRef channel
                    lineOut(3) = unique(neuroLabelsKats(obj.chsignif(ch,:,1),1));
                else 
                    [~,im]= max(cell2mat(neuroLabelsKats(:,2))); %index neuroLabel s nejvyssim rozdilem mezi kanaly
                    lineOut(3) = {['? ' neuroLabelsKats{im,1}]};
                end   
                cellOut(ch, :) =  lineOut;
            end
            
            %export tabulky
            tablelog = cell2table(cellOut,'VariableNames', variableNames); 
            [~,mfilename,~] = fileparts(obj.Eh.hfilename); 
            mfilename = strrep(mfilename, ' ', '_');
            logfilename = ['OrigChannelResponse_' mfilename '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') ];  
            xlsfilename = fullfile('logs', [logfilename '.xls']);            
            writetable(tablelog, xlsfilename); %zapisu do xls tabulky            
            disp([ 'XLS table saved: ' xlsfilename]);
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
        function obj = hybejPlotCh(obj,~,eventDat)  
            %reaguje na udalosti v grafu PlotCh
            switch eventDat.Key
                case {'rightarrow','c'} %dalsi kanal
                    if obj.PlotChCh < size(obj.chnums,1) 
                        obj.PlotCh(obj.PlotChCh+1);
                    end
                case {'leftarrow','z'} %predchozi kanal    
                    if obj.PlotChCh > 1
                        obj.PlotCh(obj.PlotChCh-1);
                    end
                case 'pagedown' %skok o 10 kanalu dopred - zvysuju cislo kanalu
                    obj.PlotCh( min( [obj.PlotChCh + 10 , size(obj.chnums,1) ])); 
                case 'pageup' %skok 10 kanalu dozadu - snizuju cislo kanalu
                    obj.PlotCh( max( [obj.PlotChCh - 10 , 1]));                    
                case 'home' %skok na prvni kanal
                    obj.PlotCh( 1);                    
                case 'end' %skok na posledni kanal
                    obj.PlotCh( size(obj.chnums,1)); 
                case 'e' %zobrazi kanal v originalnim CHilbertMulti
                    obj.Eh.PlotResponseCh(obj.PlotChCh);
                    figure(obj.PlotChH);
                case 'o' %zobrazi kanal v originalin referenci
                    obj.PlotResponseCh(obj.PlotChCh);
                    figure(obj.PlotChH);
            end
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

