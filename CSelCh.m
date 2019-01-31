classdef CSelCh < matlab.mixin.Copyable 
    %CSelCh Trida uchovavajici a poskutujici seznamy vybranych kanalu pro soubory CHIlbertMulti aj
    %   pred frekvence, kontrasty aj
    
    properties (Access = public)
        selCh;
        n; %pocet ulozenych dat
        filename;
        reference;
    end
    %#ok<*PROP>
    %#ok<*PROPLC>
    methods (Access = public)
        function obj = CSelCh(filename)
            if exist('filename','var')                
                obj.Load(filename);            
            else
                obj.selCh = cell(1,9); %sloupce filename, selCh,3:katstr=label, freq , chnames,6: chval, chlabels, 8: datestr, MNI
                obj.n = 0;
            end
        end
        function obj = SetSelCh(obj,selCh,filename,chnames, katstr,freq,chvals,neurologyLabels,MNI,reference)            
            %ulozi vyber kanalu do teto tridy            
            if ~exist('chnames','var'), chnames = []; end
            if ~exist('katstr','var'), katstr = []; end
            if ~exist('freq','var'), freq = []; end
            if isa(selCh,'CiEEGData')
                %muzu predat jako parametr celou tridu CM aj, data z nim ziskam sam
                E = selCh;
                selCh = E.GetSelCh(); %vraci E.plotRCh.selCh
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
                iCh = 1:size(selCh,1); %index kanalu, ktere ukladam 
                chnames = {E.CH.H.channels(iCh).name}'; %jmena vybranych kanaly - 2018-09-07 - ted jsou SelCH indexy, takze ukladam jmena vsech kanalu
                [katname,interval,signum] = CHilbertMulti.GetLabelInfo(E.label);
                [prumery, ~,~,~,katsnames,neurologyLabels] = E.IntervalyResp(interval,[],signum,1);    %chci vykresli obrazek        
                chvals = prumery(:,1,contains(katsnames,katname));
                MNI_xyz = cat(2,{E.CH.H.channels(iCh).MNI_x}',{E.CH.H.channels(iCh).MNI_y}',{E.CH.H.channels(iCh).MNI_z}'); %cell(ch,3)
                obj.SetSelCh(selCh,filename,chnames,katstr,freq,chvals,neurologyLabels,MNI_xyz,E.reference);
            else
                if obj.n>0
                    s = find(~cellfun(@isempty,strfind(obj.selCh(:,1),filename)),1); %najdu pouze prvni vyhovujici soubor
                else
                    s = [];
                end                 
                if isempty(s) %filename neexistuje, ulozim
                    obj.selCh(obj.n+1,:) = {filename,selCh,katstr,freq,chnames,chvals,neurologyLabels,datestr(now),MNI};
                    obj.n = obj.n +1;
                else %filename uz ulozen, radku s nim prepisu
                    obj.selCh(s,:) = {filename,selCh,katstr,freq,chnames,chvals,neurologyLabels,datestr(now),MNI}; 
                end
                if ~isempty(obj.reference)
                    assert(strcmp(obj.reference,reference),'soubory musi mit stejnou referenci');
                else
                    obj.reference = reference;
                end
                disp([ num2str(numel(find(any(selCh,2)))) ' selected channels saved']);                   
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
                     for ich = 1:size(selCh,1)
                        if ~isempty(obj.selCh{s,5})
                            if ~strcmp(selChNames{ich},E.CH.H.channels(ich).name)
                                warning(['channel ' num2str(selCh(ich)) ' nesedi s popisem: ' selChNames{ich} ' vs.' E.CH.H.channels(selCh(ich)).name]);     
                                selCh = [];
                                break;
                            end
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
            selCh = obj.selCh; %#ok<NASGU>
            n = obj.n;      %#ok<NASGU>  
            ref = obj.reference; %#ok<NASGU>  
            save(filename,'selCh','n','filename','ref','-v7.3');  
            disp(['ulozeno do ' filename]); 
        end
        function obj = Load(obj,filename)            
            assert(exist(filename,'file')==2, 'soubor neexistuje');
            load(filename,'selCh','n','filename','ref');              
            obj.n = n; %#ok<CPROPLC>
            obj.selCh = selCh; %#ok<CPROPLC>
            obj.filename = filename;
            if exist('ref','var') obj.reference = ref; else obj.reference = 'nan'; end
            disp(['nacten soubor ' filename]); 
        end
        function obj = SortByLabel(obj)
            %seradi soubory v selCh podle labels
            labels = obj.selCh(:,3);
            [~,il] = sort(labels);
            obj.selCh(1:end,:) = obj.selCh(il,:);
        end  
        function CM = LoadCM(obj,n)
            %nacte a vraci CM soubor, podle vybrane radky 
            CM = CHilbertMulti(obj.selCh{n,1});
        end
        function CM = BrainPlot(obj,n)
            %nacte CM soubor a zobrazi vcetne BrainPlotu
            CM = CHilbertMulti(obj.selCh{n,1});             
            CM.PlotResponseCh();
            CM.CH.ChannelPlot2D(1,CM.plotRCh,@CM.PlotResponseCh,obj.selCh{n,3}); %vlozim i label jako posledni parametr
            disp([ 'loaded: ' obj.selCh{n,3} ]); %vypisu co jsem vlastne nahral
        end
        function GetInfo(obj,n)
            fprintf('line %i with %i channels: %s \n', n, numel(obj.selCh{n,5}), obj.selCh{n,3} ); %vypisu co jsem vlastne nahral
        end
        function [ChNames,ChVals,ChLabels,ChNum,Intervals] = GetTables(obj,selCh,chlabels,notchnlabels,sort,normalize)
            %TODO - k selCh to bude chtit jeste parametr katstr, abych mel data napr jen z SceneXFace nebo SceneXObject
            if ~exist('selCh','var') || isempty(selCh), selCh = [1 0 0 0 0 0]; end
            selCh(numel(selCh)+1:6)= zeros(1,6-numel(selCh));            
            if ~exist('chlabels','var') || isempty(chlabels), chlabels = {}; end
            if ~exist('notchnlabels','var') || isempty(notchnlabels), notchnlabels = {}; end
            if ~exist('normalize','var') || isempty(normalize) , normalize = 1; end
            if ~exist('sort','var') || isempty(sort), sort = 1; end
            ChNames = obj.selCh{1,5};  %seznam jmen kanalu pres vsechny soubory                     
            [invervalstr,katname]=obj.GetIntervalKat(1);
            Intervals = {invervalstr}; %seznam vsech casovych intervalu v poradu podle razeni selCh
            KatNames = {katname};
            for n = 2:obj.n  
                ChNames = union(ChNames, obj.selCh{n,5});  %ziskam serazeny seznam vsech kanalu ve vsech souborech                
                [invervalstr,katname]=obj.GetIntervalKat(n);
                Intervals = union(Intervals,invervalstr);
                KatNames = union(KatNames,katname);
            end
            ChVals = zeros(numel(ChNames),numel(Intervals)); %tam budou hodnoty kanalu
            ChNum = zeros(numel(ChNames),numel(Intervals));  %pocty prekryvajicich se kontrastu v tomto kanalu a intervalu
            ChLabels = cell(numel(ChNames),1); %labely kanalu 
            ChSelections = cell(numel(ChNames),1); %oznaceni kanalu v selCh(:,2)
            for n = 1:obj.n %pro vsechny radky = importovane CM soubory               
               iI = contains(Intervals,obj.GetIntervalKat(n)); %index v poli Intervals
               for ch = 1:size(obj.selCh{n,2},1) %pro vsechny kanaly v kazdem CM souboru
                   if any(obj.selCh{n,2}(ch,logical(selCh))) %pokud je selected alespon jeden z selCh
                       iU = contains(ChNames,obj.selCh{n,5}(ch)); 
                       if  ChVals(iU,iI) == 0
                           ChVals(iU,iI) = obj.selCh{n,6}(ch);  %prvni hodnota pro tento kanal                           
                       else
                           %pokud bude vice hodnot nez dve, bude tohle spatne
                           ChVals(iU,iI) = mean([obj.selCh{n,6}(ch) ChVals(iU,iI)]);  %prumer pro rozdil vuci jedne a druhe kategorii
                       end
                       ChNum(iU,iI) = ChNum(iU,iI) + 1;
                       if isempty(ChLabels{iU})
                           ChLabels{iU} = obj.selCh{n,7}{ch}; %staci jednou pro kazdy kanal
                       end
                       ChSelections{iU} = union(ChSelections{iU},obj.GetSelectedMarks(n,ch),'stable'); %TODO - najde to prekryv v odpovedi v jakemkoliv case
                   end
               end
            end
            
            [ChVals,ChNames,ChNum,ChLabels,ChSelections,notLstr] = obj.ChValsFilter(ChVals,ChNames,ChNum,ChLabels,ChSelections,chlabels,notchnlabels);
            if sort == 1
                [ChVals,ChNames,ChNum,ChLabels,ChSelections] = obj.ChValsSort(ChVals,ChNames,ChNum,ChLabels,ChSelections);
            end
            if normalize == 1
                ChVals = obj.ChValsNormalize(ChVals); 
            end
            
            fh = figure('Name',['CSelCh GetTables ' cell2str(chlabels)]);
            ax1 = axes('Position',[0 0 1 1],'Visible','off'); 
            ax2 = axes('Position',[0.18 0.18 0.75 0.75]); % https://nl.mathworks.com/help/matlab/creating_plots/placing-text-outside-the-axes.html
            imagesc(ax2,ChVals);            
            title([cell2str(chlabels) notLstr ' - ' num2str(selCh,'%.0f-%.0f-%.0f-%.0f-%.0f-%.0f')]);
            set(gca, 'YTick', 1:numel(ChNames), 'YTickLabel', ChNames) % 10 ticks 
            set(gca, 'XTick', 1:numel(Intervals), 'XTickLabel', Intervals) % 10 ticks 
%             if numel(chlabels) > 0
                for iL = 1:numel(ChLabels)
                    text(0.6,iL,ChLabels{iL},'Color','white');
                    text(7,iL,cell2str(ChSelections{iL}),'Color','white');
                end
%             end
            colorbar;
            axes(ax1); % sets ax1 to current axes
            text(0.1,0.1, cell2str(KatNames) ); %pozici jsem si vyzkousel empiricky
            
        end
        function ExtractBrainPlotData(obj)
            BPD = struct; 
            %{ 
            BPD.intervals - matrix
            BPD.katstr - cell katstr
            BPD.reference
            BPD.Hf 50:5:150
            BPD.selCh cell(intervals,kat)
            BPD.NAMES cell(intervals,kat)
            !!! BPD.MNI
            BPD.VALS
            BPD.NLABELS 
            xx BPD.EPI
            %}
            
        end
    end
     %  --------- privatni metody ----------------------
    methods (Access = private)
        function [intervalstr,katname] = GetIntervalKat(obj,n)
           [katname,interval,~] = CHilbertMulti.GetLabelInfo(obj.selCh{n,3}); %inverval prvniho souboru
           intervalstr = num2str(interval,'%.1f-%.1f');
        end    
        function SelMarks = GetSelectedMarks(obj,n,ch)
            %vraci znacky S F O pro jeden kanal podle selCh indexu 1-6 (=fghjkl)
            %podle PPA testu - f-znamena odpoved na prvni podminka kontrastu a ostatni na treti podminku neuvedenou v kontrastu
           [katname,~,~] = CHilbertMulti.GetLabelInfo(obj.selCh{n,3});
           C = strsplit(katname,'X');  %katname je napriklad FaceXScene
           Marks = {C{1}(1),C{2}(1)}; %prvni znaky z Face a Scene
           if Marks{1} == 'S'
              Marks{2} = regexprep(Marks{2},{'F','O','X'},{'X','F','O'}); %nejak to funguje, asi to zkousi trikrat dokolecka, dokud ma co nahrazovat, takze proto X
           elseif Marks{1} == 'O' 
               Marks{2} = regexprep(Marks{2},{'F','S','X'},{'X','F','S'});
           elseif Marks{1} == 'F' 
               Marks{2} = regexprep(Marks{2},{'S','O','X'},{'X','S','O'});
           end            
           SelMarks = Marks( [obj.selCh{n,2}(ch,1)>0 , any(obj.selCh{n,2}(ch,2:6)>0) ]); %vraci napr {F S}
        end
    end
    methods (Access = private, Static)
        function [ChVals,ChNames,ChNum,ChLabels,ChSelections,notLstr]=ChValsFilter(ChVals,ChNames,ChNum,ChLabels,ChSelections,chlabels,notchnlabels)
             %odfiltruju kanaly bez oznaceni pres vsechny intervaly
            iChVals = max(ChVals,[],2)==0; 
            ChVals(iChVals,:) = [];
            ChNames(iChVals) = [];
            ChNum(iChVals,:) = [];
            ChLabels(iChVals) = [];
            ChSelections(iChVals) = [];           
            
            if numel(chlabels) > 0
                iL = contains(ChLabels,chlabels);
                if numel(notchnlabels) > 0
                    iLx = contains(ChLabels,notchnlabels);
                    iL = iL & ~iLx;
                    notLstr = [ ' not:' cell2str(notchnlabels)];
                else
                    notLstr = '';
                end
                ChVals = ChVals(iL,:);
                ChNames = ChNames(iL);
                ChNum = ChNum(iL,:);
                ChLabels = ChLabels(iL);
                ChSelections = ChSelections(iL); 
            else                
                notLstr = '';
            end
        end
        function [ChVals,ChNames,ChNum,ChLabels,ChSelections] = ChValsSort(ChVals,ChNames,ChNum,ChLabels,ChSelections)
            %seradim kanaly podle casu a velkosti odpovedi            
            [Max,iMax] = max(ChVals,[],2); %Nejvyssi odpoved v kazdem radku  + jeji index v radku
            [iMax2,iiMax] = sort(iMax); %serazene indexy max odpovedi (prvni v radku, druhe v radku ....) a jejich indexy 
            MaxVals = [iMax2, Max(iiMax), iiMax]; %serazene indexy max odpovedi + ty maximalni odpovedi serazene podle jejich indexu v radku
            MaxVals = sortrows(MaxVals,[1 -2]); %ziskam serazene iiMax v tretim sloupci
            ChVals = ChVals(MaxVals(:,3),:);
            ChNames = ChNames(MaxVals(:,3));
            ChLabels = ChLabels(MaxVals(:,3));
            ChSelections = ChSelections(MaxVals(:,3));
            ChNum = ChNum(MaxVals(:,3),:);
        end
        function ChVals = ChValsNormalize(ChVals)
            for ch = 1:size(ChVals,1)
                iCh = ChVals(ch,:)>0;
                ChVals(ch,iCh) = ChVals(ch,iCh)./max(ChVals(ch,iCh));
            end
        end
    end
end

