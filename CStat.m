classdef CStat < handle
    %CSTAT Class for custom statistical functions
    %   Kamil Vlcek, FGU AVCR, since 2016 06
    
    properties   (Access = public) 
        plotAUC % data k auc plotu 
        %{ 
        Eh  - odkaz na celou instanci CiEEGData
        setup -  nastaveni, vznika v AUTPlotIni a nemeni se
        kategories - E.Wp.kats (cisla kategorii ve statistice)
        time - [max(E.baseline(2),E.epochtime(1)) E.epochtime(2)]
        sig - signifikance auc krivek - kanaly x kombinace kategorii, plni se v ROCAnalysis
        aucdata - struct array (channels) se sloupci AUC a AVG
        selCh - kopie E.plotRCh.selCh
        %}
        plotAUC_m; %data k aucplotM plotu
    end
    
    methods (Access = public)
        function obj = CStat(~)
            obj.AUCReset();
            obj.AUCPlotIni();
        end
        function obj = AUCPlotIni(obj)
            %nastavi defaultni hodnoty pro graf AUCPlot            
            obj.plotAUC.corelplot = 0; %defaulte hezci plot, ale nemozny do corelu
            hue = 0.8;
            obj.plotAUC.setup.colorskat = {[0 0 0],[0 1 0],[1 0 0],[0 0 1]; [hue hue hue],[hue 1 hue],[1 hue hue],[hue hue 1]}; % prvni radka - prumery, druha radka errorbars = svetlejsi            
            obj.plotAUC.setup.colorkomb = [nan 2 1; 2 nan 4; 1 4 nan]; %index barvy kombinace kategorii
            obj.plotAUC.setup.legendkomb = [nan 1 2 ; 1 nan 3 ; 2 3 nan ]; % do ktereho pole legendy se ma ukladat kombinace kategorii                                 
            obj.plotAUC.katplot = ones(1,max(max(obj.plotAUC.setup.legendkomb))); %vic kategorii nikdy nebude - ktera kombinace kategorii se maji kreslit
        end
        function obj = AUCReset(obj)
            obj.plotAUC.aucdata = struct('AUC',{},'AVG',{}); %empty struct array            
            obj.plotAUC.sig = [];
        end
        
        function [obj] = AUCPlot(obj,ch,E, time,kategories)
            %vykresli AUC krivku pro vybrany kanal. 
            %pouzije data z plotAUC, pokud je potreba zavola funkci ROCAnalysis
            %time je vyhodnocovany cas ve vterinach, ale vzdy se bere z E
            %kategories jsou E.Wp.kats (cisla kategorii ve statistice), ale parametr se nikdy v kodu nepouziva, bere se z E            
            
            %overim jestli se nezmenily parametry krivek - zatim mi staci kategorie, kvuli zmene statistiky
            if isfield(obj.plotAUC,'kategories') 
                if exist('kategories','var') && ~isempty(kategories) && ~isequal(obj.plotAUC.kategories,kategories)
                    obj.AUCReset();               
                elseif exist('E','var') && ~isequal(obj.plotAUC.kategories,E.Wp(E.WpActive).kats)
                    obj.AUCReset();               
                end               
            end
                        
            if ch > numel(obj.plotAUC.aucdata) || isempty( obj.plotAUC.aucdata(ch).AUC)
                 %assert(~isempty(E),['pro kanal' num2str(ch) 'nejsou spocitana data']);
                 if ~exist('time','var'), time = []; end
                 if ~exist('kategories','var'), kategories = []; end
                 [AUCall,AVGall] = obj.ROCAnalysis(E,ch,time,kategories);                 
                 PsyData = E.PsyData;
                 kategories = obj.plotAUC.kategories;
                 time = obj.plotAUC.time;                 
            else
                if ~isfield(obj.plotAUC,'Eh') %pokud jsem nacetl CS.plotAUC ze souboru, nejsou v nem pole Eh a PsyData
                    obj.plotAUC.Eh = E;
                    obj.plotAUC.PsyData = E.PsyData;
                end
                AUCall = obj.plotAUC.aucdata(ch).AUC;
                AVGall = obj.plotAUC.aucdata(ch).AVG;
                if ~exist('time','var') || isempty(time), time = obj.plotAUC.time; end
                if ~exist('kategories','var') || isempty(kategories), kategories = obj.plotAUC.kategories; end
                if ~exist('E','var') || isempty(E)
                    PsyData = obj.plotAUC.PsyData; 
                else
                    PsyData = E.PsyData;
                    obj.plotAUC.selCh = E.plotRCh.selCh; %vyber kanalu fghjkl * pocet kanalu - updatuju pri kazdem vykresleni auc
                end                
            end
            obj.plotAUC.ch = ch;
            
            %barvy jako v PlotResponseCh            
            katnames = PsyData.CategoryName(kategories,[]);
            legenda = cell(1,numel(kategories));
            obj.plotAUC.katnames = cell(1,numel(kategories));
            
            if isfield(obj.plotAUC,'fh') && (verLessThan('matlab','9.0') || isvalid(obj.plotAUC.fh)) %isvalid je od verze 2016
                figure(obj.plotAUC.fh); %pouziju uz vytvoreny graf
                clf(obj.plotAUC.fh); %graf vycistim
            else               
               figurename = ['AUC waveform for ch ' num2str(ch) ];               
               obj.plotAUC.fh = figure('Name',figurename); %ulozim si pojmenovani kategorii 
            end
                        
            for k = 1:numel(kategories)-1
            for l = k+1:numel(kategories)
                sig = iff(obj.plotAUC.sig(ch,obj.plotAUC.setup.legendkomb(k,l)),'*',''); %hvezdicku pro signifikatni AUC krivku
                obj.plotAUC.katnames{obj.plotAUC.setup.legendkomb(k,l)} = [katnames{l} ' X ' katnames{k}];
                if obj.plotAUC.katplot(obj.plotAUC.setup.legendkomb(k,l)) > 0 %pokud se tahle kombinace kategorii ma kreslit                    
                    legenda{obj.plotAUC.setup.legendkomb(k,l)} = [ obj.plotAUC.katnames{obj.plotAUC.setup.legendkomb(k,l)}  ' ' sig];                

                    AUC = AUCall{k,l};
                    MP = AVGall{k,l};                
                    MN = AVGall{l,k};

                    %kod podle PlotResponseCh, aby stejne barvy pro PPA test
                    colorkat_kl = obj.plotAUC.setup.colorkomb(kategories(k), kategories(l));
                    color_kl = cell2mat(obj.plotAUC.setup.colorskat(:,colorkat_kl));                    

                    %PRVNI plot je AUC krivka
                    subplot(2,1,1); 
                    X = linspace(time(1),time(2),size(AUC,1));
                    hold on;
                    if obj.plotAUC.corelplot ==1
                        ciplot(AUC(:,2), AUC(:,3), X,color_kl(2,:)); %nepruhledny, ale jde okopirovat do corelu
                    else
                        plotband(X, AUC(:,1), AUC(:,3) - AUC(:,1), color_kl(2,:)); %nejlepsi, je pruhledny, ale nejde kopirovat do corelu
                    end
                    plot(X,AUC(:,1),'.-','Color',color_kl(1,:));
                    line([X(1) X(end)],[.5 .5]);
                    title(['AUC for channel ' num2str(ch) '/' num2str(obj.plotAUC.channels)]);
                    ylim([0.3 1]);

                    %do DRUHEHO plot vykreslim rozdil power obou kategorii
                    subplot(2,1,2); 
                    title('power');
                    hold on;

                    %prvni a druha kategorie - jejich rozdil                
                    plot(X,MP-MN,'o-','Color',color_kl(1,:)); %barvy podle PlotResponseCh 
                end
            end
            end
            legenda = legenda(~cellfun('isempty',legenda)); %ymazu prazdne polozky, ktere se nevykresluji
            if ~isempty(legenda), legend(legenda); end
            set(obj.plotAUC.fh,'KeyPressFcn',@obj.hybejAUCPlot); 
        end
        function [AUCall,AVGall,obj] = ROCAnalysis(obj,E, channels,time,kategories)
            %vypocita ROC data a ulozi do  obj.plotAUC.aucdata
            %time je cas, ve kterem chci spocitat ROC ve vterinach
            %kategories je seznam cisel kategorii, pro jejichz kombinace se ma vypocitat a vykreslit AUC krivka
            %ch je seznam kanalu k vyhodnoceni, pokud prazne, spocita pro vsechny kanaly                      
                       
            if ~exist('channels','var') || isempty(channels), channels = 1:E.channels; end
            if ~exist('time','var') || isempty(time), time = [max(E.baseline(2),E.epochtime(1)) E.epochtime(2)]; end
            if ~exist('kategorie','var') || isempty(kategorie), kategories = E.Wp(E.WpActive).kats; end
            if ~isfield(obj.plotAUC,'sig') || isempty(obj.plotAUC.sig)
                obj.plotAUC.sig = zeros(E.channels,max(max(obj.plotAUC.setup.legendkomb))); %signifikance auc krivek - kanaly x kombinace kategorii
            end
            sample = round((time - E.epochtime(1))*E.fs); %cisla vzorku v E.d od-do
            katnames = E.PsyData.CategoryName(kategories,[]); %do lengendy grafu
            
            fprintf('computing AUC data for channels (of %u): ',numel(channels));
            
            obj.plotAUC.PsyData = E.PsyData;
            obj.plotAUC.time = time;
            obj.plotAUC.kategories = kategories;            
            obj.plotAUC.channels = E.channels;  % pocet kanalu
            obj.plotAUC.selCh = E.plotRCh.selCh; %vyber kanalu fghjkl * pocet kanalu
            obj.plotAUC.selChNames = E.plotRCh.selChNames; %vyber kanalu fghjkl * pocet kanalu
            obj.plotAUC.Eh = E; %handle na celou strukturu, abych mohl volat ROCAnalysis z AUCPlotM
            
            for ch = channels %jde po sloupcich
                fprintf('%u,',ch);
                %musim vyradit spatne epochy            
                AUCall = cell(numel(kategories)); %tam bud davat AUC data pro vsechny kombinace kategorii
                AVGall = cell(numel(kategories)); %tam budou prumery rozdilu mezi kategoriemi
                
                for k = 1:numel(kategories)-1 %index nizsi kategorie
                for l = k+1:numel(kategories) %index vyssi kategorie - jsou serazene ve statistice tak, ze dulezitejsi maji vyssi cisla, napr 2=Face x 3=Object x 1=Scene = [2 3 1]

                    [~,~,~,iEpP] = E.CategoryData(kategories(l),[],[],ch); %ziskam parametr iEp se seznamem vsech validnich epoch pro tento kanal
                    [~,~,~,iEpN] = E.CategoryData(kategories(k),[],[],ch); %.... chci mit tu nejdulezitejsi kategorii (l=vyssi cislo) jako prvni, aby se od ni odecitaly ostatni                                                        

                    %prvni a druha kategorie - prumer power
                    dataP = squeeze(E.d(sample(1):sample(2),ch,iEpP));
                    dataN = squeeze(E.d(sample(1):sample(2),ch,iEpN));
                    MP = mean(dataP,2); %prumer pres epochy - zustavaji samples jako rozmer
                    MN = mean(dataN,2);                   

                    if numel(sample) == 1 %chci udela ROC krivku jen pro jeden bod v case
                        AUC = CStat.ROCKrivka(E.epochData(iEpP | iEpN,:),squeeze(E.d(sample,ch,iEpP | iEpN)),{katnames{l},katnames{k}},1); %udela i graf
                    else %udelam ROC krivku pro ze vsech auc
                        AUC = zeros(diff(sample)+1,3); %samples x auc hodnota + confidence intervals
                        for s = sample(1):sample(2)
                            auc = CStat.ROCKrivka(E.epochData(iEpP | iEpN,:),squeeze(E.d(s,ch,iEpP | iEpN)),{katnames{l},katnames{k}},0); %zadny graf nedela
                            ci =  CStat.AUCconfI(auc,[sum(iEpP) sum(iEpN)],0.05);
                            AUC(s-sample(1)+1,1:3) = [auc ci];
                        end                                                           
                    end
                    AUCall{k,l} = AUC; %samples x [auc, dolni ci, horni ci]                   
                    obj.plotAUC.sig(ch,obj.plotAUC.setup.legendkomb(k,l)) = any( all(AUC(:,2:3)>0.5,2) | all(AUC(:,2:3)<0.5,2) ); %krivka AUC je signifikancni, pokud oba CI prekroci 0.5 jednom nebo druhym smerem
                    AVGall{k,l} = MP; % pozitivni data - prvni vyssi kategorie - l
                    AVGall{l,k} = MN; % negativni data - druha nizsi kategorie - k
                end
                end
                obj.plotAUC.aucdata(ch).AUC = AUCall;
                obj.plotAUC.aucdata(ch).AVG = AVGall;
            end
            fprintf('... done \n');
        end
        function obj = AUCPlotM(obj,channels,chSelection,selch)
            %vykresli vic AUC krivek pres sebe z vybranych kanalu
            %channels je seznam kanalu k vykresleni
            %chSelection je pojmenovani vyberu kanalu, napr [Scene]
            %selch je cislo kanalu z chsort(channels), ktery je zobrazen vyrazne
            assert(sum(obj.plotAUC.katplot)==1,'one kontrast for plot is not selected'); %nechci kreslit vic kontrastu jako ScenexObject a ScenexFace dohromady
            
            params = {'channels','chSelection','selch'}; %zkusim hromadne zpracovani parametru 
            for p = 1:numel(params)            
                if ~exist(params{p},'var')  || eval(['isempty(' params{p} ')'])
                    eval([params{p} '=' 'obj.plotAUC_m.' params{p} ';']); %touhle velmi nedoporucovanou metodou
                else
                    eval(['obj.plotAUC_m.' params{p} '='  params{p} ';']  );            
                end
            end            
                            
            chantodo = channels(channels>numel(obj.plotAUC.aucdata)); %spocitam si nespocitane kanaly predem najednou. Zatim neresim spocitani i tech chybejicich kanalu
            if ~isempty(chantodo)
                obj.ROCAnalysis(obj.plotAUC.Eh,chantodo,obj.plotAUC.time,obj.plotAUC.kategories); 
            end
            sig = obj.plotAUC.sig(channels,logical(obj.plotAUC.katplot)); %index kanalu se signifik auc krivkami ke kresleni
            
            legenda = cell(1,sum(sig));                        
            ColorSet = distinguishable_colors(sum(sig)); %ruzne barvy pro ruzne kanaly
            ploth = zeros(1,sum(sig)); %handly na ploty, kvuli selektivni legende            
                        
            if isfield(obj.plotAUC_m,'fh') && (verLessThan('matlab','9.0') || isvalid(obj.plotAUC_m.fh)) %isvalid je od verze 2016
               figure(obj.plotAUC_m.fh); %pouziju uz vytvoreny graf
               clf(obj.plotAUC_m.fh); %graf vycistim
               figurenew = 0;  %obnovuju uz drive vytvoreny graf=figure            
               selchReal = obj.plotAUC_m.chsort(selch); %je vybrany kanal podle velikosti AUC krivky               
            else                                  
               obj.plotAUC_m.fh = figure('Name','AUC for multiple channels');
               figurenew = 1; %vytvoril jsem novy graf - okno
               selchReal = 0; %neni zadny vybrany kanal
            end
            
            if exist('chSelection','var'), ChSelText = [' chnls: ' cell2str(obj.plotAUC.selChNames{chSelection}) ]; else, ChSelText = ''; end
            figuretitle= ['AUCPlotM kontrast: ' obj.plotAUC.katnames{find(obj.plotAUC.katplot)}  ChSelText   ]; %#ok<FNDSB>
            if figurenew, disp(figuretitle); end            
            ileg = 1; %specialni index na signif kanaly - legendu a barvy            
            for ch = 1:numel(channels)                
                if channels(ch) > numel(obj.plotAUC.aucdata) || isempty( obj.plotAUC.aucdata(channels(ch)).AUC)                     
                    obj.ROCAnalysis(obj.plotAUC.Eh,channels(ch),obj.plotAUC.time,obj.plotAUC.kategories);                                        
                end                
                for k = 1:numel(obj.plotAUC.kategories)-1
                for l = k+1:numel(obj.plotAUC.kategories)
                    if obj.plotAUC.katplot(obj.plotAUC.setup.legendkomb(k,l)) > 0 ... %pokud se tahle kombinace kategorii ma kreslit
                    && ~isempty( obj.plotAUC.aucdata(channels(ch)).AUC) %a AUC data existuji
                        
                        AUC = obj.plotAUC.aucdata(channels(ch)).AUC{k,l};
                        X = linspace(obj.plotAUC.time(1),obj.plotAUC.time(2),size(AUC,1));                                             
                        if obj.plotAUC.sig(channels(ch),obj.plotAUC.setup.legendkomb(k,l)) && sig(ch) % nechci mit true u dodatecne pocitanych kanalu (ROCAnalysis kousek vyse) , || figure new tu nemuzu dat, protoze pokud neni ch prvni, nejsou vsechny sig naplnene od zacatku
                            legenda{ileg} = ['ch' num2str(channels(ch)) ' ' obj.plotAUC.Eh.CH.H.channels(channels(ch)).name];
                            style = iff(any(AUC(:,2)>0.75) || any(AUC(:,3)<0.25),'o-','.-');
                            LineWidth = iff(ch==selchReal,3,1); %vybrany kanal je nakresleny tluste
                            ploth(ileg) = plot(X,AUC(:,1),style,'Color',ColorSet(ileg,:),'LineWidth',LineWidth); %kazdy kanal jinou barvou, pokud budu kreslit vic krivek pro jeden kanal, bude to asi zmatek
                            if ch==selchReal, selchH = ploth(ileg); end                                                        
                            ileg = ileg + 1;
                        else  %nesignif  auc krivka
                            style = iff(ch==selchReal,'.-','.:'); %ybrany kanal neni teckovane, ostani ano
                            ph = plot(X,AUC(:,1),style,'Color',[.5 .5 .5]); %seda barva
                            if ch==selchReal, selchH = ph; end
                        end
                        hold on;                                                
                    end
                end
                end
            end
            line([X(1) X(end)],[.5 .5]);             
            if selch>0 %kdyz poprve graf vykreslim, neni zadny vybrany kanal
                uistack(selchH, 'top');  %vybrany kanal dam na popredi
                txt = sprintf('ch: %i(%i), %s: %s, %s',channels(selchReal),selch, obj.plotAUC.Eh.CH.H.channels(channels(selchReal)).name, ...
                    obj.plotAUC.Eh.CH.H.channels(channels(selchReal)).neurologyLabel,obj.plotAUC.Eh.CH.H.channels(channels(selchReal)).ass_brainAtlas );
                text(.05,.1,txt);
            end
            legenda = legenda(~cellfun('isempty',legenda)); %vymazu prazdne polozky, ktere se nevykresluji
            legend(ploth,legenda);
            ylim([0 1]);
            
            title(figuretitle);           
            if figurenew %jen pokud kreslim graf poprve
                disp(['significant channels: ' num2str(channels(logical(sig)),'%i ')]);
                disp(['not significant channels: ' num2str(channels(~logical(sig)),'%i ')]);
                obj.AUCPlotMMax(channels);    %spocitam si razeni kanalu podle maxima AUC krivek            
            end
            set(obj.plotAUC_m.fh,'KeyPressFcn',@obj.hybejAUCPlotM); 
        end
        function [chsort] = AUCPlotMMax(obj,channels)
            %chci ziskat poradi kanalu podle maxima AUC krivky
            chmax = zeros(numel(channels),1);
            katplotind = find(obj.plotAUC.setup.legendkomb == find(obj.plotAUC.katplot));
            [l,k] = ind2sub(size(obj.plotAUC.setup.legendkomb),katplotind(1)); %ziskam kategorie, mezi kterymi jsou vykreslene kontrasty
            for ch =1:numel(channels)
                aucdata = obj.plotAUC.aucdata(channels(ch)).AUC{k,l}; %kvuli ladeni
                [~,imax] = max(abs(aucdata(:,1)-.5)); %maximum z auc krivky, nikoli ci 
                chmax(ch) = aucdata(imax,1)-.5;
            end
            [~,chsort] = sort(chmax, 'descend');
            obj.plotAUC_m.chsort = chsort;
            obj.plotAUC_m.chmax = chmax; %ulozim i hodnoty, asi nemuzu ukladat pro vsechny kanaly, protoze pak by neplatilo chsort
        end
        function AUCPlotBrain(obj,selch,vals)
            %volam funkci na vykresleni 3D obrazku mozku ChannelPlot
            %vals - muzu dodat hodnoty na vykresleni, defaultne jsou pouzite maxima AUC krivek. 
            if ~exist('vals','var')
                vals = abs(obj.plotAUC_m.chmax)+.5; %chmax jsou hodnoty -.5 az .5. Chci zobrazovat negativni rozliseni jako pozitivni
            end 
            obj.plotAUC.Eh.CH.ChannelPlot([],0,vals,... %param chnvals
                obj.plotAUC_m.channels,... %chnsel jsou cisla kanalu, pokud chci jen jejich vyber
                obj.plotAUC_m.chsort(selch)); %selch je jedno zvyraznene cislo kanalu - index v poli chnsel
        end       
      
        function AUC2XLS(obj)
            %vypise seznam kanalu z grafu AUCPlotM do xls souboru 
            %vola se pomoci x z grafu
            channels = obj.plotAUC_m.channels;   %XXX: predpokladam, ze obj.plotAUC_m.channels uz obsahuje vsechny spocitane kanaly            
            cellout = cell(numel(channels),15); % z toho bude vystupni xls tabulka s prehledem vysledku            
            for ch = 1:numel(channels)  %XXX: iterace pres kanaly, predpokladam, ze je jen jedna platna kombinace {k,l} nize!            
                channelHeader = obj.plotAUC.Eh.CH.H.channels(channels(ch));                
                for k = 1:numel(obj.plotAUC.kategories)-1
                for l = k+1:numel(obj.plotAUC.kategories)
                    if obj.plotAUC.katplot(obj.plotAUC.setup.legendkomb(k,l)) > 0 ... %pokud se tahle kombinace kategorii ma kreslit
                    && ~isempty( obj.plotAUC.aucdata(channels(ch)).AUC) %a AUC data existuji
                        
                        AUC = obj.plotAUC.aucdata(channels(ch)).AUC{k,l};
                        X = linspace(obj.plotAUC.time(1),obj.plotAUC.time(2),size(AUC,1));
                        
                        [amax, idx, idxHalf] = cMax(AUC(:,1));
                        tmax = X(idx);
                        thalf = X(idxHalf);                        
                        ci_u = AUC(idx, 3);
                        ci_l = AUC(idx, 2);                     
                        sig = obj.plotAUC.sig(channels(ch), obj.plotAUC.setup.legendkomb(k,l));
                        
                        cellout(ch, :) = {  channels(ch),   channelHeader.name, channelHeader.neurologyLabel, channelHeader.MNI_x, ...
                                channelHeader.MNI_y, channelHeader.MNI_z, channelHeader.seizureOnset, channelHeader.interictalOften, ...
                                mat2str(channelHeader.rejected), tmax, thalf, amax, ci_u, ci_l, sig };                                          
                    end
                end
                end                
            end
            
            tablelog = cell2table(cellout, ...
                'VariableNames', {'channel' 'name'  'neurologyLabel'  'MNI_x'  'MNI_y'  'MNI_z'  'seizureOnset'  'interictalOften'  ...
                    'rejected'  'tmax'  'thalf'  'aucmax'  'ci_u'  'ci_l'  'significance'   
                });
            obj.plotAUC_m.xlsvals = cell2mat(cellout(:,10:12)); %ulozim hodnoty tmax, thalf a aucmax
            %TODO: Identifikace nazvu souboru? 
            kat = strrep([obj.plotAUC.katnames{find(obj.plotAUC.katplot)}], ' ', '_'); %#ok<FNDSB>
            chnls = regexprep(cell2str(obj.plotAUC.selChNames{obj.plotAUC_m.chSelection}), {' ','[',']'}, {'_','(',')'}); %writetable cant use [] in filenames
            [~,mfilename,~] = fileparts(obj.plotAUC.Eh.hfilename);
            mfilename = strrep(mfilename, ' ', '_');
            logfilename = ['AUCPlotM_' kat '_chnls_' chnls '_' mfilename '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') ];  
            xlsfilename = fullfile('logs', [logfilename '.xls']);            
            writetable(tablelog, xlsfilename); %zapisu do xls tabulky            
            disp([ 'xls tables saved: ' xlsfilename]);
        end
    end
    methods (Static,Access = public)        
        function W = Wilcox2D(A,B,print,fdr,msg,RjEpChA,RjEpChB)
            %Wilcox2D(A,B,print,fdr,msg,RjEpChA,RjEpChB)
            %srovna dve 3D matice proti sobe, ohledne hodnot v poslednim rozmeru
            %A musi mit oba prvni rozmery > rozmery B, 
            %B muze mit jeden nebo oba prvni rozmer = 1 - pak se porovnava se vsemi hodnotami v A
            %pokud fdr=1 nebo prazne, provadi fdr korekci
            %pokud print = 0 nebo prazne, netiskne nic
            if ~exist('print','var'), print = 0; end
            if ~exist('fdr','var') || isempty(fdr), fdr = 1; end %min striktni je default           
            if ~exist('msg','var'), msg = ''; end
            if ~exist('RjEpChA','var'), RjEpChA = false(size(A,2),size(A,3)); end
            if ~exist('RjEpChB','var'), RjEpChB = false(size(B,2),size(B,3)); end
            W = zeros(size(A,1),size(A,2));
           
            if print, fprintf(['Wilcox Test 2D - ' msg ': ']); end
            for j = 1:size(A,1) % napr cas
                if print && mod(j,50)==0, fprintf('%d ', j); end %tisknu jen cele padesatky
                for k = 1:size(A,2) %napr kanaly                  
                   aa = squeeze (A(j,k,~RjEpChA(k,:))); %jen nevyrazene epochy 
                   bb = squeeze (B( min(j,size(B,1)) , min(k,size(B,2)) , ~RjEpChB(k,:) )); %jen nevyrazene epochy
                   if numel(aa) >= 2 && numel(bb) >= 2 
                      W(j,k) = ranksum(aa,bb); % Statistics and Machine Learning Toolbox
                   else
                      W(j,k) = 1; %pokud jen jedna hodnota, nelze delat statistika
                   end
                end
            end
            if fdr
                if fdr == 2, method = 'dep'; else method = 'pdep'; end
                [~, ~, adj_p]=fdr_bh(W,0.05,method,'no'); %dep je striktnejsi nez pdep
                W = adj_p; %prepisu puvodni hodnoty korigovanymi podle FDR
            end
            if print, fprintf('%d .. done\n',j); end
        end
        function p_corrected = PermStat(A,B,print,msg,RjEpChA,RjEpChB) %#ok<INUSD>
            %P = PermStat(A,B,print,msg,RjEpChA,RjEpChB)            
            % A a B jsou cas x channels x epochs
            % RjEpChA a B jsou channels x epochs
            % funkce od Radka Bortela pro permutacni statistiku - 29.1.2018
            if ~exist('print','var'), print = 0; end
            if ~exist('msg','var'), msg = ''; end
            if ~exist('RjEpChA','var'), RjEpChA = false(size(A,2),size(A,3)); end %#ok<NASGU>
            if ~exist('RjEpChB','var'), RjEpChB = false(size(B,2),size(B,3)); end %#ok<NASGU>                          
            %vyrazovani epoch pro jednotlive kanaly se zatim nepouziva, funkce to neumoznuje
            if (size(A,3)==0 || size(B,3)==0 ) 
                warning('Permutation Test 2D - zadne epochy');                
            end
          
            if print, fprintf(['Permutation Test 2D - ' msg ': ']); end
            timer = tic; %zadnu merit cas
            p_corrected=CmpPerm(A,B,30000,2000); %pocty permutaci, cisla doporucena od Radka
            cas = toc(timer); %ukoncim mereni casu a vypisu
            if print, fprintf(' .. done in %.1f min \n',cas/60); end            
        end
        function [W] = Klouzaveokno(W,oknosirka, funkce,dimension)
            % oknosirka je v poctu bodu, funkce muze byt min, max, mean
            % 27.4.2015 - vynato z wilcoxmap, 21.6.2016 presunuto do CStat      
                        
            if oknosirka >= 1 
                if dimension == 1 %pokud chci pocitat klouzave okno v prvnim rozmeru, transponuju na zacatku i na konci
                    W = W';
                end                
                pulsirka = ceil(oknosirka/2); %ktery sloupec se povazuje za pulku okna - tam se hodnota ulozi
                W2 = zeros(size(W,1),size(W,2)); %musim udelat kopii, jinak si prepisuju hodnoty ze kterych pak pocitam
                for sloupec = 1:size(W,2) 
                    iW = max([1 sloupec-pulsirka+1]) : min([size(W,2) sloupec-pulsirka+oknosirka]); 
                    switch funkce
                        case 'min'
                            W2(:,sloupec)=min(W(:,iW),[],2);
                        case 'max'
                            W2(:,sloupec)=max(W(:,iW),[],2);
                        case 'mean'
                            W2(:,sloupec)=mean(W(:,iW),2);
                        otherwise
                            W2(:,sloupec)=W(:,sloupec); %pokud jina funkce, nemenim matici
                    end 
                end
                if dimension == 1 
                    W = W2';
                else
                    W = W2;
                end
            end
        end
        
        function N = round(n,dig)
            %zakrouhuje na dany pocet desetinnych mist
            N = round(n * 10^dig) / 10 ^ dig;
        end
        function [frequencies,fft_d_abs]=Fourier(dd,fs,method)
            %[frequencies,fft_d_abs]=Fourier(dd,fs) 
            % spocita Fourierovu fransformaci dat z jednoho kanalu 
            %predpoklada data s casem v prvnim rozmeru, ostatni rozmery nejsou nebo 1
            % 4.5.2017 kopie z CiEEGdata.Fourier
            if ~exist('method','var'), method = 'pwelch'; end 
            switch method   
                
                case 'fft'
            %1. varianta podle MikeXCohen
            frequencies       = linspace(0,fs/2,length(dd)/2+1); %maximalni frekvence je fs/2, ale frekvencni rozliseni je N/2+1
            fft_d    = fft(dd)/length(dd); %fast fourier transform
           
            %vezmu jen tolik frekvenci, kolik je v frequencies - realnych frekvenci. ostatni jsou imaginarni frekvence
            fft_d_abs = abs(fft_d(1:length(frequencies)))*2; %dvema nasobim kvuli tem imaginarnim frekvencim. Viz MikeCohenP.            
            %prvni frekvence je DC
            fft_d_abs = 10*log10(fft_d_abs);
            
                case 'periodogram'
            %2. varianta z Matlab Answers https://uk.mathworks.com/matlabcentral/answers/114942-how-to-calculate-and-plot-power-spectral-density-of-a-given-signal
            %ale je to v podstate stejne jako fft
            [pxx,frequencies] = periodogram(dd,[],length(dd),fs);  %tohle je jen abs(fft())^2/n
            fft_d_abs = 10*log10(pxx);
            
                case 'pwelch'
            %3. pouziti pwelch, viz https://stackoverflow.com/questions/27079289/on-the-use-and-understanding-of-pwelch-in-matlab            
            M = 8*fs; %round(length(dd)/20); 
            [pxx, frequencies] = pwelch(dd,M,round(M/2),fs*20,fs); %2.-4. parametr od Radka Bortela 29.06.2017
            fft_d_abs = 10*log10(pxx);
            end
        end
            
        function [filter_result] = FIR(freq,dd,fs,firtype)
            %[filter_result] = FIR(freq,dd,fs) 
            % provede FIR filter podle Cohen Ch 14.  
            % 4.5.2017   
            if ~exist('firtype','var'), firtype = 'fir1'; end
            assert(strcmp(firtype,'fir1') || strcmp(firtype,'firls'), ['neznamy typ kernelu ' firtype]);
            
            nyquist =fs/2;                      
                        
            if strcmp(firtype,'fir1')
                %FIR1 - jednodussy kernel doporuceny v eeglab
                if numel(freq)==1 %highpass
                    filter_order       = round(3*(fs/freq(1))); %trojnasobek dolni frekvence - delka filter kernelu 
                    filterweights = fir1(filter_order,freq./nyquist, 'high'); %vytvorim filter kernel 
                elseif freq(1)==0 %lowpass
                    filter_order       = round(3*(fs/freq(2))); %trojnasobek dolni frekvence - delka filter kernelu 
                    filterweights = fir1(filter_order,freq(2)./nyquist); %vytvorim filter kernel 
                else
                    filter_order       = round(3*(fs/freq(1))); %trojnasobek dolni frekvence - delka filter kernelu 
                    filterweights = fir1(filter_order,freq./nyquist); %vytvorim filter kernel
                end
                %kdyz u fir1 dam jen jednu frekvenci, automaticky pocita lowpass - viz eeglab:eegfilt
            else
                %FIRLS
                %freqspread    = (freq(2)-freq(1))/2; % Hz +/- the center frequency
                %center_freq   = freq(1) + freqspread;   
                assert(numel(freq)>1 && freq(1)>0,'firls muze byt jen bandpass');
                filter_order       = round(3*(fs/freq(1))); %trojnasobek dolni frekvence - delka filter kernelu 
                transwid      = .10; % transition zone withth
                ffrequencies  = [ 0 (1-transwid)*(freq(1)) (freq(1)) (freq(2)) (1+transwid)*(freq(2)) nyquist ]/nyquist;
                ffrequencies( ffrequencies>1 ) = 1; %pokud jsem nastavil jako pasmo=nyquist, zlobilo by to jinak
                ffrequencies( ffrequencies<0 ) = 0; %pokud jsem nastavil jako pasmo=0, zlobilo by to jinak
                idealresponse = [ 0 0 1 1 0 0 ];
                filterweights = firls(filter_order,ffrequencies,idealresponse); %vytvorim filter kernel                          
            end
            filter_result = filtfilt(filterweights,1,dd);
            
        end
        
        function [AUC] = ROCKrivka(epochData,dd,katnames,kresli)
            %11/2018 - spocita a vykresli roc krivku
            %epochData - cellarray z CiEEGData
            %dd - eegdata pro jeden channel a jeden sample
            %sample - cislo sample
            %katnames - jmena dvou kategorii, positive a negative
            if ~exist('kresli','var'), kresli = 0; end
            labels = epochData(:,1);
            scores = dd;
            posclass = katnames{1}; %napr Scene
            if numel(katnames) > 1
                negclass = katnames{2}; %napr Face
            else
                negclass = 'all';
            end
            
            %labels = cell array/int array nazvu kategorii - co radek to event, napr Scene, Object, Face
            %scores = array se scores pro jednotlive - co radek to event - power - prumer nebo urcity cas
            %posclass = string - Positive class label - nazev vyhodnocovane kategorie
            %options - 'NegClass','versicolor' - negative class, default je all
            
            [X,Y,~,AUC] = perfcurve(labels,scores,posclass,'NegClass',negclass);            
            if kresli
                figure('Name',['ROC curve for ' posclass 'x' negclass ]);
                subplot(1,2,1);
                plot(X,Y);     
                subplot(1,2,2);
                posdata = dd(not(cellfun('isempty',strfind(epochData(:,1),posclass))));
                negdata = dd(not(cellfun('isempty',strfind(epochData(:,1),negclass))));            
                plot(rand(size(posdata))*0.5+1,posdata,'or');
                hold on;
                plot(rand(size(negdata))*0.5+2,negdata,'ob');
            end
            
        end
        
        function ci = AUCconfI(auc,n,p)
            % http://www.real-statistics.com/descriptive-statistics/roc-curve-classification-table/auc-confidence-interval/
            %n jsou velikosti dvou vzorku
            %p = 0.05 napr
            
            q0 = auc*(1-auc);
            q1 = auc/(2-auc)-auc^2;
            q2 = 2*auc^2/(1+auc)-auc^2;
            se = sqrt( (q0 + (n(1)-1)*q1 + (n(2)-1)*q2) / (n(1)*n(2)) );
            zcrit = norminv(1-p/2); %two tailed z critical value from p value
            ci = [auc - se*zcrit , auc + se*zcrit];
        end
    end
    methods (Access = private)
        function obj = hybejAUCPlot(obj,~,eventDat)
            switch eventDat.Key
                case {'u','i','o','p'} %jina pismena nez f-l, aby se to nepletlo                       
                    ik = find('uiop'==eventDat.Key); %index 1-4
                    obj.plotAUC.katplot( ik ) = 1 - obj.plotAUC.katplot(ik);
                    obj.AUCPlot(obj.plotAUC.ch);                      
                case {'c'}
                    obj.plotAUC.corelplot = 1 - obj.plotAUC.corelplot;
                    obj.AUCPlot(obj.plotAUC.ch);
                case {'f','g','h','j','k','l'}                    
                    channels = find(obj.plotAUC.selCh(:,'fghjkl'==eventDat.Key))'; %cisla musi byt v radce
                    obj.AUCPlotM(channels,find('fghjkl'==eventDat.Key),0); %#ok<FNDSB> %povinne ted uvadim predvybrany kanal
            end            
        end    
        function obj = hybejAUCPlotM(obj,~,eventDat)
            kresli = 1;
            switch eventDat.Key                
                case 'uparrow' % dalsi kanal v poradi
                    selch = max(1,obj.plotAUC_m.selch-1);                     
                case 'downarrow' %predchozi kanal v poradi
                    selch = min(numel(obj.plotAUC_m.channels),obj.plotAUC_m.selch+1);                    
                case 'home' % dalsi kanal v poradi   
                    selch = 1;                    
                case 'end' % dalsi kanal v poradi 
                    selch = numel(obj.plotAUC_m.channels);                    
                case 'pageup' % dalsi kanal v poradi
                    selch = max(1,obj.plotAUC_m.selch-10);                    
                case 'pagedown' %predchozi kanal v poradi
                    selch = min(numel(obj.plotAUC_m.channels),obj.plotAUC_m.selch+10);   
                case 'return' %zobrazeni mozku
                    selch = max(1,obj.plotAUC_m.selch); %cislo jednoho vybraneho kanalu: musim neco priradit - puvodni kanal, ale ne 0
                    obj.AUCPlotBrain(selch);
                case 'x' %export kanalu do xls tabulky
                    obj.AUC2XLS(); 
                    kresli = 0;
                otherwise     
                    kresli = 0;
            end
            if kresli 
                obj.AUCPlotM([],[],selch); %prekresli sumarni plot AUC krivek
                obj.plotAUC.Eh.PlotResponseCh(obj.plotAUC_m.channels(obj.plotAUC_m.chsort(selch))); %prekresli graf PlotResponseCh
                figure(obj.plotAUC_m.fh);
            end
        end
        
    end
    
end

