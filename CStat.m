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
        function c = ColorKomb(obj,kat1,kat2) 
            %nova funkce kvuli jinym kategoriim nez u PPA plotu. U Menrot jsou treba kategorie od 0
            if min([kat1 kat2])==0
                kategorie = [kat1 kat2]+1;
                kat1 = kategorie(1:numel(kat1));
                kat2 = kategorie(numel(kat1)+1:end);
            end
            c = obj.plotAUC.setup.colorkomb(kat1(1), kat2(1)); %pokud jsou kat1 a kat2 arrays, vratim jen prvni hodnotu
        end
        function [obj] = AUCPlot(obj,ch,E, time,kategories)
            %vykresli AUC krivku pro vybrany kanal. 
            %pouzije data z plotAUC, pokud je potreba zavola funkci ROCAnalysis
            %time je vyhodnocovany cas ve vterinach, ale vzdy se bere z E
            %kategories jsou E.Wp.kats (cisla kategorii ve statistice), ale parametr se nikdy v kodu nepouziva, bere se z E            
            assert(numel(ch)==1, 'only one channels should be specified');
            
            if ~exist('E','var') || isempty(E)
                E = obj.plotAUC.Eh;
                ch = obj.plotAUC.Eh.CH.sortorder(ch);%cislo kanalu musi odpovidat aktualnimu filtru a razeni
            else
                ch = E.CH.sortorder(ch);                
            end
            if ~isfield(obj.plotAUC,'reversekats') || isempty(obj.plotAUC.reversekats), obj.plotAUC.reversekats = 0; end %if the categories should be in original or reversed order
            
            if exist('time','var') && numel(time) == 1 %chci jen vykreslit jednu ROC krivku projeden bod
                obj.ROCAnalysis(E,ch,time,kategories); 
                return;
            end
            %overim jestli se nezmenily parametry krivek - zatim mi staci kategorie, kvuli zmene statistiky
            if isfield(obj.plotAUC,'kategories') 
                if exist('kategories','var') && ~isempty(kategories) && ~isequal(obj.plotAUC.kategories,kategories) %if there is argument kategories and there are different categories in stored ROC results
                    obj.AUCReset();               
                elseif exist('E','var') && ~isequal(obj.plotAUC.kategories,E.Wp(E.WpActive).kats) %if the kategories in active stat are different from saved categories
                    obj.AUCReset();               
                end               
            end
            if ch > numel(obj.plotAUC.aucdata) || obj.AUCisEmpty(ch)
%                 || isempty( obj.plotAUC.aucdata(ch).AUC) ... %if no data from this channel are computed                    
%                  || isempty(iff(obj.plotAUC.reversekats,obj.plotAUC.aucdata(ch).AUC{2,1},obj.plotAUC.aucdata(ch).AUC{1,2})) %or if they are computed in other direction (SxO vs OxS)
                 %the ROC data for this channel should be computed
                 %assert(~isempty(E),['pro kanal' num2str(ch) 'nejsou spocitana data']);
                 if ~exist('time','var'), time = []; end
                 if ~exist('kategories','var'), kategories = []; end
                 [AUCall,AVGall] = obj.ROCAnalysis(E,ch,time,kategories);  %computes the ROC data, and stores them to obj.plotAUC.aucdata(ch)             
                 PsyData = E.PsyData;
                 kategories = obj.plotAUC.kategories;
                 time = obj.plotAUC.time;                 
            else %the ROC data are computed already
                if ~isfield(obj.plotAUC,'Eh') %pokud jsem nacetl CS.plotAUC ze souboru, nejsou v nem pole Eh a PsyData
                    obj.plotAUC.Eh = E;
                    obj.plotAUC.PsyData = E.PsyData;
                end 
                AUCall = obj.plotAUC.aucdata(ch).AUC; % cellarray: kats x kats - load the previously computed AUC data
                AVGall = obj.plotAUC.aucdata(ch).AVG; % cellarray: kats x kats - difference between kats
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
            
            for kk = 1:numel(kategories)-1
            for ll = kk+1:numel(kategories)
                if obj.plotAUC.reversekats, [l, k] = deal(kk,ll); else,  [k, l] = deal(kk,ll); end %swap both variables content
                sig = iff(obj.plotAUC.sig(ch,obj.plotAUC.setup.legendkomb(k,l)),'*',''); %hvezdicku pro signifikatni AUC krivku
                obj.plotAUC.katnames{obj.plotAUC.setup.legendkomb(k,l)} = [katnames{l} ' X ' katnames{k}];
                if obj.plotAUC.katplot(obj.plotAUC.setup.legendkomb(k,l)) > 0 %pokud se tahle kombinace kategorii ma kreslit                    
                    legenda{obj.plotAUC.setup.legendkomb(k,l)} = [ obj.plotAUC.katnames{obj.plotAUC.setup.legendkomb(k,l)}  ' ' sig];                

                    AUC = AUCall{k,l};
                    MP = AVGall{k,l};  %samples x 1 (epochs mean)              
                    MN = AVGall{l,k};

                    %kod podle PlotResponseCh, aby stejne barvy pro PPA test
                    colorkat_kl = obj.ColorKomb(cellval(kategories,k), cellval(kategories,l));
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
                    contrastKeys = 'uiop';
                    text(0.1,0.9,contrastKeys(logical(obj.plotAUC.katplot)));
                    if obj.plotAUC.reversekats, text(0.7,0.9,'r'), end
                    
                    %do DRUHEHO plotu vykreslim rozdil power obou kategorii
                    subplot(2,1,2); 
                    title('power difference');
                    hold on;

                    %prvni a druha kategorie - jejich rozdil                
                    plot(X,MP-MN,'o-','Color',color_kl(1,:)); %barvy podle PlotResponseCh 
                end
            end
            end
            legenda = legenda(~cellfun('isempty',legenda)); %vymazu prazdne polozky, ktere se nevykresluji
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
            if ~exist('kategories','var') || isempty(kategories), kategories = E.Wp(E.WpActive).kats; end
            if ~isfield(obj.plotAUC,'sig') || isempty(obj.plotAUC.sig)
                obj.plotAUC.sig = zeros(E.channels,max(max(obj.plotAUC.setup.legendkomb))); %signifikance auc krivek - kanaly x kombinace kategorii
            end
            sample = round((time - E.epochtime(1))*E.fs); %cisla vzorku v E.d od-do
            katnames = E.PsyData.CategoryName(kategories,[]); %do lengendy grafu
            
            fprintf('computing AUC data for channels (of %u): ',numel(channels));
            if(numel(sample)>1)
                obj.plotAUC.PsyData = E.PsyData;
                obj.plotAUC.time = time;
                obj.plotAUC.kategories = kategories;            
                obj.plotAUC.channels = E.channels;  % pocet kanalu
                obj.plotAUC.selCh = E.plotRCh.selCh; %vyber kanalu fghjkl * pocet kanalu
                obj.plotAUC.selChNames = E.plotRCh.selChNames; %vyber kanalu fghjkl * pocet kanalu
                obj.plotAUC.Eh = E; %handle na celou strukturu, abych mohl volat ROCAnalysis z AUCPlotM
            end
            for ch = channels %jde po sloupcich
                fprintf('%u,',ch);
                %musim vyradit spatne epochy            
                if ~isfield(obj.plotAUC.aucdata,'AUC') || isempty(obj.plotAUC.aucdata)  || length(obj.plotAUC.aucdata) < ch || isempty(obj.plotAUC.aucdata(ch).AUC)  %are there any computed AUC values for this channel?
                    AUCall = cell(numel(kategories)); %tam bud davat AUC data pro vsechny kombinace kategorii
                    AVGall = cell(numel(kategories)); %tam budou prumery rozdilu mezi kategoriemi
                else
                    AUCall = obj.plotAUC.aucdata(ch).AUC;
                    AVGall = obj.plotAUC.aucdata(ch).AVG; 
                end
                
                for kk = 1:numel(kategories)-1 %index nizsi kategorie
                for ll = kk+1:numel(kategories) %index vyssi kategorie - jsou serazene ve statistice tak, ze dulezitejsi maji vyssi cisla, napr 2=Face x 3=Object x 1=Scene = [2 3 1]
                    if obj.plotAUC.reversekats, [l, k] = deal(kk,ll); else,  [k, l] = deal(kk,ll); end %swap both variables content
                    
                    [~,~,~,iEpP] = E.CategoryData(cellval(kategories,l),[],[],ch); %ziskam parametr iEp se seznamem vsech validnich epoch pro tento kanal
                    [~,~,~,iEpN] = E.CategoryData(cellval(kategories,k),[],[],ch); %.... chci mit tu nejdulezitejsi kategorii (l=vyssi cislo) jako prvni, aby se od ni odecitaly ostatni                                                        

                    %prvni a druha kategorie - prumer power
                    if(numel(sample)>1) %from time=sample range - normal case
                        dataP = squeeze(E.d(sample(1):sample(2),ch,iEpP));
                        dataN = squeeze(E.d(sample(1):sample(2),ch,iEpN));
                    else %chci mit zachovanou moznot vyplotit ROC krivku pro jeden sample
                        dataP = squeeze(E.d(sample(1),ch,iEpP));
                        dataN = squeeze(E.d(sample(1),ch,iEpN));
                    end
                    MP = mean(dataP,2); %prumer pres epochy - zustavaji samples jako prvni rozmer
                    MN = mean(dataN,2);                   

                    if numel(sample) == 1 %chci udela ROC krivku jen pro jeden bod v case
                        auc = CStat.ROCKrivka(E.epochData(iEpP | iEpN,:),squeeze(E.d(sample,ch,iEpP | iEpN)),{katnames{l},katnames{k}},1); %udela i graf
                        ci =  CStat.AUCconfI(auc,[sum(iEpP) sum(iEpN)],0.05);
                        %AUC = [auc ci];
                        disp(['AUC = ' num2str(auc) '+-' num2str(auc-ci(1)) '=[' num2str(ci(1)) ';' num2str(ci(2)) ']']);
                    else %udelam ROC krivku for all time samples in the range sample(1) - (2)
                        AUC = zeros(diff(sample)+1,3); %samples x auc hodnota + confidence intervals
                        for s = sample(1):sample(2)
                            auc = CStat.ROCKrivka(E.epochData(iEpP | iEpN,:),squeeze(E.d(s,ch,iEpP | iEpN)),{katnames{l},katnames{k}},0); %zadny graf nedela
                            ci =  CStat.AUCconfI(auc,[sum(iEpP) sum(iEpN)],0.05);
                            AUC(s-sample(1)+1,1:3) = [auc ci];
                        end                                                           
                        AUCall{k,l} = AUC; %samples x [auc, dolni ci, horni ci]                                           
                        AVGall{k,l} = MP; % pozitivni data - samples x 1 (epochs mean) = of first higher category l
                        AVGall{l,k} = MN; % negativni data = of second lower category k  - lower means earlier in kategories vector
                        %the significance of AUC curve - if both confidence intervals are over or below 0.5, in the same direction
                        obj.plotAUC.sig(ch,obj.plotAUC.setup.legendkomb(k,l)) = any( all(AUC(:,2:3)>0.5,2) | all(AUC(:,2:3)<0.5,2) ); %krivka AUC je signifikancni, pokud oba CI prekroci 0.5 jednom nebo druhym smerem                    
                    end                    
                end
                end
                if(numel(sample)>1)
                    obj.plotAUC.aucdata(ch).AUC = AUCall;
                    obj.plotAUC.aucdata(ch).AVG = AVGall;
                end
            end
            fprintf('... done \n');
        end
        function obj = AUCPlotM(obj,channels,chSelection,selch)
            %vykresli vic AUC krivek pres sebe z vybranych kanalu
            %channels je seznam kanalu k vykresleni
            %chSelection je pojmenovani vyberu kanalu, napr [Scene], index v obj.plotAUC.selChNames, 0 pro all channels
            %selch je cislo kanalu z channels(chsort) = chnums, ktery je zobrazen vyrazne            
            assert(sum(obj.plotAUC.katplot)==1,'one kontrast for plot is not selected'); %nechci kreslit vic kontrastu jako ScenexObject a ScenexFace dohromady
            
            if ~isempty(channels) %has to be above the arguments initialization below (params=)
                channels = intersect(obj.plotAUC.Eh.CH.sortorder,channels); %chci jen kanaly, ktere odpovidaji filtru podle sortorder
            end
            
            % POZOR - parameters inicialization %
            params = {'channels','chSelection','selch'}; %zkusim hromadne zpracovani parametru 
            for p = 1:numel(params)            
                if ~exist(params{p},'var')  || eval(['isempty(' params{p} ')']) 
                    eval([params{p} '=' 'obj.plotAUC_m.' params{p} ';']); %touhle velmi nedoporucovanou metodou
                else
                    eval(['obj.plotAUC_m.' params{p} '='  params{p} ';']  );            
                end
            end            
            if ~isfield(obj.plotAUC_m,'plotmean'), obj.plotAUC_m.plotmean = 0; end
            if isfield(obj.plotAUC_m,'fh') && (verLessThan('matlab','9.0') || isvalid(obj.plotAUC_m.fh)) %isvalid je od verze 2016
               figure(obj.plotAUC_m.fh); %pouziju uz vytvoreny graf
               clf(obj.plotAUC_m.fh); %graf vycistim
               figurenew = 0;  %obnovuju uz drive vytvoreny graf=figure                           
            else                                  
               obj.plotAUC_m.fh = figure('Name','AUC for multiple channels');
               figurenew = 1; %vytvoril jsem novy graf - okno               
            end            
            
            chantodo = obj.AUCisEmpty(channels); %channels where the AUC data are missing
%             chantodo = find(cellfun(@isempty,{obj.plotAUC.aucdata.AUC})); %najdu kanaly, ktere jsou nespocitane v aucdata
%             chantodo = find(cellfun(@isempty,{obj.plotAUC.aucdata.AUC{1,2}}));
%             chantodo = intersect(chantodo,channels); %z tech nespocitanych chci jen ty, ktere se maji ted vykreslit
            %chantodo = channels(channels>numel(obj.plotAUC.aucdata)); %spocitam si nespocitane kanaly predem najednou. Zatim neresim spocitani i tech chybejicich kanalu
            if ~isempty(chantodo)
                obj.ROCAnalysis(obj.plotAUC.Eh,chantodo,obj.plotAUC.time,obj.plotAUC.kategories); 
            end                     
            if isfield(obj.plotAUC_m,'chsort') &&  figurenew == 0 %pokud jsou kanaly serazene jinak a neni to novy obrazek
                chnums = channels(obj.plotAUC_m.chsort); %only reorder the channels
            else
                chnums = channels;            %channels are absolute ch numbers
            end
%             sig = obj.plotAUC.sig(chnums,logical(obj.plotAUC.katplot)); % ch x 1 double - index kanalu se signifik auc krivkami ke kresleni
                    
            if exist('chSelection','var') && ~isempty(obj.plotAUC.selChNames)
                if ~isempty(chSelection)  && chSelection > 0
                    ChSelText = [' chnls: ' cell2str(obj.plotAUC.selChNames{chSelection}) ]; %list of channels to the figure title
                else
                    ChSelText = ' chnls: all'; % vsechny kanaly zobrazuju, nefiltruju je podle selCh
                end
            else
                ChSelText = ''; 
            end
            figuretitle= ['AUCPlotM kontrast: ' obj.plotAUC.katnames{find(obj.plotAUC.katplot)}  ChSelText   ]; %#ok<FNDSB>
            if figurenew, disp(figuretitle); end            
            ileg = 1; %specialni index na signif kanaly - legendu a barvy            
            [AUCData,AUCSig]=obj.AUCPlotMGetData(chnums); %AUCData: time x 3 x channels, AUCSig: channels x 1
            T = linspace(obj.plotAUC.time(1),obj.plotAUC.time(2),size(AUCData,1));  %1 x time in sec                                           
            if obj.plotAUC_m.plotmean
                %mean of the significant channels
                %TODO merge with obj.AUCPlotMean
                labels = {obj.plotAUC.Eh.CH.brainlabels(chnums).label};
                [ulabels,barvy]=obj.plotAUC.Eh.CH.channelPlot.GetBrainlabelsSaved();
                legendStr = cell(numel(ulabels),1);
                ploth = zeros(1,numel(ulabels)); 
                colorsErrorBars = obj.plotAUC.Eh.PL.GetColorsErrorBars(barvy); 
                for j = 1:numel(ulabels)
                    ich = contains(labels,ulabels{j});
                    if sum(ich)>0
                        M = mean(AUCData(:,1,ich),3);
                        E = std(AUCData(:,1,ich),[],3)/sqrt(sum(ich)); %std err of mean                        
                        ciplot(M+E, M-E, T,colorsErrorBars(j,:) ); %colorkatk(2,:)
                        %plotband(T, M, E, colorsErrorBars(j,:)); %nejlepsi, je pruhledny, ale nejde kopirovat do corelu
                        hold on;
                        ploth(j) = plot(T,M,'LineWidth',2,'Color',barvy(j,:));
                        legendStr{j}=[ ulabels{j} ' x' num2str(sum(ich))];
                    end
                end
            else
                legendStr = cell(1,sum(AUCSig));                        
                ColorSet = distinguishable_colors(sum(AUCSig)); %ruzne barvy pro ruzne kanaly
                ploth = nan(1,sum(AUCSig)); %handly na ploty, kvuli selektivni legende 
                for ch = 1:size(AUCData,3)                   
                    if AUCSig(ch) % nechci mit true u dodatecne pocitanych kanalu (ROCAnalysis kousek vyse) , || figure new tu nemuzu dat, protoze pokud neni ch prvni, nejsou vsechny sig naplnene od zacatku
                        legendStr{ileg} = ['ch' num2str(chnums(ch)) ' ' obj.plotAUC.Eh.CH.H.channels(chnums(ch)).name];
                        style = iff(any(AUCData(:,2,ch)>0.75) || any(AUCData(:,3,ch)<0.25),'o-','.-'); %channels with high values (>0.75) plotted with circles
                        LineWidth = iff(ch==selch,3,1); %vybrany kanal je nakresleny tluste
                        ploth(ileg) = plot(T,AUCData(:,1,ch),style,'Color',ColorSet(ileg,:),'LineWidth',LineWidth); %kazdy kanal jinou barvou, pokud budu kreslit vic krivek pro jeden kanal, bude to asi zmatek
                        if ch==selch, selchH = ploth(ileg); end                                                        
                        ileg = ileg + 1;
                    else  %nesignif  auc krivka
                        style = iff(ch==selch,'.-','.:'); %ybrany kanal neni teckovane, ostani ano
                        ph = plot(T,AUCData(:,1,ch),style,'Color',[.5 .5 .5]); %seda barva
                        if ch==selch, selchH = ph; end
                    end
                    hold on;                                                
                end

                if selch>0 % if there is some selected channel - kdyz poprve graf vykreslim, neni zadny vybrany kanal
                    uistack(selchH, 'top');  %vybrany kanal dam na popredi
                    chnum = channels(obj.plotAUC_m.chsort(selch));
                    txt = sprintf('ch: %i(%i), %s: %s, %s',chnum,selch, obj.plotAUC.Eh.CH.H.channels(chnum).name, ...
                        obj.plotAUC.Eh.CH.H.channels(chnum).neurologyLabel,obj.plotAUC.Eh.CH.H.channels(chnum).ass_brainAtlas );
                    text(.05,.1,txt);
                    if isfield(obj.plotAUC_m,'xlsvals')
                        if size(obj.plotAUC_m.xlsvals,1) ==  numel(obj.plotAUC_m.channels)
                            txt = sprintf('aucmax: %.3f, tmax %.3f',obj.plotAUC_m.xlsvals(selch,1),obj.plotAUC_m.xlsvals(selch,2) );                        
                        else
                            txt = 'xlsvals have different no of values. Run AUC2XLS to recompute to see aucmax and tmax value.';
                        end
                        text(.05,.05,txt);
                    end
                end
                
                if figurenew %jen pokud kreslim graf poprve
                    disp(['significant channels: ' num2str(channels(logical(AUCSig)),'%i ')]);
                    disp(['not significant channels: ' num2str(channels(~logical(AUCSig)),'%i ')]);
                    obj.AUCPlotMMax(channels);    %spocitam si razeni kanalu podle maxima AUC krivek            
                end
            end
            ylim([0.4 1]);
            line([T(1) T(end)],[.5 .5]);             
%             legendStr = legendStr(~cellfun('isempty',legendStr)); %vymazu prazdne polozky, ktere se nevykresluji
%             ploth = ploth(~isnan(ploth)); %vymazu prazdne polozky, ktere se nevykresluji
            emptyh = ploth==0;
            ploth(emptyh)=[];
            legendStr(emptyh)=[];
            if ~isempty(legendStr) 
                legend(ploth,legendStr,'Location','best'); 
            end
            title(figuretitle);           
            set(obj.plotAUC_m.fh,'KeyPressFcn',@obj.hybejAUCPlotM); 
        end
        function [AUCData,AUCSig]=AUCPlotMGetData(obj,chnums)          
            %gets auc data to be plotted - individually or averaged
            AUCData=zeros(size(obj.plotAUC.aucdata(chnums(1)).AUC{1,2},1),3,numel(chnums));
            AUCSig = false(numel(chnums,1));
            for ch = 1:numel(chnums)                     
                chnum = chnums(ch); %absolute channel number - only shortcut               
                if chnum > numel(obj.plotAUC.aucdata) || isempty( obj.plotAUC.aucdata(chnum).AUC)                     
                    obj.ROCAnalysis(obj.plotAUC.Eh,chnum,obj.plotAUC.time,obj.plotAUC.kategories); %if there are still missing values for this channel, compute them now
                end                
                for kk = 1:numel(obj.plotAUC.kategories)-1
                for ll = kk+1:numel(obj.plotAUC.kategories)
                    if obj.plotAUC.reversekats, [l, k] = deal(kk,ll); else,  [k, l] = deal(kk,ll); end %swap both variables content
                    if obj.plotAUC.katplot(obj.plotAUC.setup.legendkomb(k,l)) > 0 ... %pokud se tahle kombinace kategorii ma kreslit
                        && ~isempty( obj.plotAUC.aucdata(chnum).AUC) %a AUC data existuji
                        
                        AUC = obj.plotAUC.aucdata(chnum).AUC{k,l}; %time x 3 - auc curve and confidence intervals
                        AUCData(:,:,ch)=AUC; %collect the data                        
                        AUCSig(ch) = obj.plotAUC.sig(chnum,obj.plotAUC.setup.legendkomb(k,l)); % nechci mit true u dodatecne pocitanych kanalu (ROCAnalysis kousek vyse) , || figure new tu nemuzu dat, protoze pokud neni ch prvni, nejsou vsechny sig naplnene od zacatku
                                    %sig = obj.plotAUC.sig(chnums,logical(obj.plotAUC.katplot));
                    end
                end
                end
            end
        end
        function [chsort] = AUCPlotMMax(obj,channels)
            %chci ziskat poradi kanalu podle maxima AUC krivky
            chmax = zeros(numel(channels),1);
            katplotind = find(obj.plotAUC.setup.legendkomb == find(obj.plotAUC.katplot));
            [ll,kk] = ind2sub(size(obj.plotAUC.setup.legendkomb),katplotind(1)); %ziskam kategorie, mezi kterymi jsou vykreslene kontrasty           
            if obj.plotAUC.reversekats, [l, k] = deal(kk,ll); else,  [k, l] = deal(kk,ll); end %swap both variables content
            for ch =1:numel(channels)
                aucdata = obj.plotAUC.aucdata(channels(ch)).AUC{k,l}; %kvuli ladeni
                [~,imax] = max(abs(aucdata(:,1)-.5)); %maximum z auc krivky, nikoli ci 
                chmax(ch) = aucdata(imax,1)-.5;
            end
            [~,chsort] = sort(chmax, 'descend');
            obj.plotAUC_m.chsort = chsort;
            obj.plotAUC_m.chmax = chmax; %ulozim i hodnoty, asi nemuzu ukladat pro vsechny kanaly, protoze pak by neplatilo chsort
        end
        function AUCPlotBrain(obj,selch,vals,rangeZ)
            %volam funkci na vykresleni 3D obrazku mozku ChannelPlot
            %selch je cislo kanalu podle poradi (podle velikosti chmax)
            %vals - muzu dodat hodnoty na vykresleni, defaultne jsou pouzite maxima AUC krivek. 
            if ~exist('selch','var') || isempty(selch)
                selch = 1; %index in obj.plotAUC_m.chsort = channels sorted by AUC value - first is max                
            end
            if ~exist('vals','var') || isempty(vals)
                sig = logical(obj.plotAUC.sig(obj.plotAUC_m.channels,logical(obj.plotAUC.katplot)));  %jestli jsou AUCkrivky signifikangni
                vals = obj.plotAUC_m.chmax+.5; %chmax jsou hodnoty -.5 az .5. chciz rozsah 0-1                         
                selchval = vals(obj.plotAUC_m.chsort(selch)); %zjistim hodnotu, kterou chci v mozku oznacit
                vals = vals(sig); %vezem jen signif kanaly
                channels =  obj.plotAUC_m.channels(sig);                
            else
                channels =   obj.plotAUC_m.channels;
                selchval = vals(selch);
            end           
            if ~exist('rangeZ','var') || isempty(rangeZ) 
               if isfield(obj.plotAUC.Eh.CH.channelPlot.plotCh3D,'rangeZ')
                   rangeZ=obj.plotAUC.Eh.CH.channelPlot.plotCh3D.rangeZ;  %copy rangeZ from existing channelplot 3D
               else
                   rangeZ=[0 1]; 
               end
            end
            
            selchs = find(vals==selchval);  %najdu znovu index hodnoty ze signif kanalu
            obj.plotAUC.Eh.CH.ChannelPlotProxy(vals,... %param chnvals
                channels,... %chnsel jsou cisla kanalu, pokud chci jen jejich vyber
                iff(~isempty(selchs),selchs,0),[],...%selch je jedno zvyraznene cislo kanalu - index v poli chnsel
                'AUCPlotBrain', rangeZ); 
            set(obj.plotAUC.Eh.CH.channelPlot.plotCh3D.fh, 'WindowButtonDownFcn', {@obj.hybejPlot3Dclick, selch});
        end
        function AUC2XLS(obj, val_fraction, int_fraction)
            %vypise seznam kanalu z grafu AUCPlotM do xls souboru 
            %vola se pomoci stlaceni x z grafu AUCPlotM
            
            %pokud neni specifikovan parametr 'fraction', zobrazi se dialogove okno pro zadani procent z maxima            
            if nargin == 1
                prompt = {'Value trigger percentage:', 'Integral trigger percentage:'};
                default = {'90', '50'};
                percent = inputdlg(prompt, 'XLS Export', [1 30], default);
                if isempty(percent)
                    disp('XLS export cancelled');
                    return;
                end
                val_fraction = str2double(percent{1})/100; % procenta auxmax, u kterych se ma zjistit cas
                int_fraction = str2double(percent{2})/100; % procenta plochy pod krivkou, u kterych se ma zjistit cas
            end
            channels = obj.plotAUC_m.channels;   %XXX: predpokladam, ze obj.plotAUC_m.channels uz obsahuje vsechny spocitane kanaly
            cellout = cell(numel(channels),16); % z toho bude vystupni xls tabulka s prehledem vysledku
            for ch = 1:numel(channels)  %XXX: iterace pres kanaly, predpokladam, ze je jen jedna platna kombinace {k,l} nize!
                channelHeader = obj.plotAUC.Eh.CH.H.channels(channels(ch));                
                for kk = 1:numel(obj.plotAUC.kategories)-1
                for ll = kk+1:numel(obj.plotAUC.kategories)
                    if obj.plotAUC.reversekats, [l, k] = deal(kk,ll); else,  [k, l] = deal(kk,ll); end %swap both variables content
                    
                    if obj.plotAUC.katplot(obj.plotAUC.setup.legendkomb(k,l)) > 0 ... %pokud se tahle kombinace kategorii ma kreslit
                    && ~isempty( obj.plotAUC.aucdata(channels(ch)).AUC) %a AUC data existuji
                        
                        AUC = obj.plotAUC.aucdata(channels(ch)).AUC{k,l}; %columns - AUC, lower CI, upper CI
                        X = linspace(obj.plotAUC.time(1),obj.plotAUC.time(2),size(AUC,1));
                        %TODO use only time, when categories are significantly different
                        [amax, idx, idxFrac] = cMax(AUC(:,1), val_fraction, 0.5); %Nalezne maximum / minimum krivky curve, i jeho podilu (napr poloviny) a jeho parametry
                        tmax = X(idx);
                        thalf = X(idxFrac); %cas poloviny maxima, nebo jineho podilu                       
                        ci_u = AUC(idx, 3); %confidence intervals - upper
                        ci_l = AUC(idx, 2);                     
                        sig = obj.plotAUC.sig(channels(ch), obj.plotAUC.setup.legendkomb(k,l));
                        
                        tint = cIntegrate(X, AUC(:,1), int_fraction, 2, 0.5); % integrace s posunem minima krivky do nuly
                        
                        cellout(ch, :) = {  channels(ch),   channelHeader.name, channelHeader.neurologyLabel, channelHeader.MNI_x, ...
                                channelHeader.MNI_y, channelHeader.MNI_z, channelHeader.seizureOnset, channelHeader.interictalOften, ...
                                mat2str(channelHeader.rejected),  amax, ci_l, ci_u, sig , tmax, thalf, tint  };                                          
                    end
                end
                end                
            end
            
            tablelog = cell2table(cellout, ...
                'VariableNames', {'channel' 'name'  'neurologyLabel'  'MNI_x'  'MNI_y'  'MNI_z'  'seizureOnset'  'interictalOften'  ...
                    'rejected' 'aucmax' 'ci_l' 'ci_u'   'significance'  'tmax'  ['t' percent{1}]   ['tint' percent{2}]
                });
            obj.plotAUC_m.xlsvals = cell2mat(cellout(:,[10 14 15])); %ulozim hodnoty aucmax, tmax, thalf 
            %TODO: Identifikace nazvu souboru? 
            kat = strrep([obj.plotAUC.katnames{find(obj.plotAUC.katplot)}], ' ', '_'); %#ok<FNDSB>
            if obj.plotAUC_m.chSelection==0
                selChNames = 'all';
            else
                selChNames = obj.plotAUC.selChNames{obj.plotAUC_m.chSelection}; %iff nejde pouzit, protoze se vyhodnotit oba argumenty vzdy
            end            
            chnls = regexprep(cell2str(selChNames), {' ','[',']'}, {'_','(',')'}); %writetable cant use [] in filenames
            [~,mfilename,~] = fileparts(obj.plotAUC.Eh.hfilename);
            mfilename = strrep(mfilename, ' ', '_');
            logfilename = ['AUCPlotM_' kat '_chnls_' chnls '_' mfilename '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') ];  
            xlsfilename = fullfile('logs', [logfilename '.xls']);            
            writetable(tablelog, xlsfilename); %zapisu do xls tabulky            
            disp([ 'XLS table saved: ' xlsfilename]);
        end
        function Scatter(obj,names)
            %vykresli scatterplot z hodnot valmax a tmax, exportovanych drive do xls tabulky
            assert(isfield(obj.plotAUC_m,'xlsvals') && ~isempty(obj.plotAUC_m.xlsvals),'no xls data');
            if ~exist('names','var'), names = 1; end
            figure('Name','AUCPlotM Scatter');
            scatter(obj.plotAUC_m.xlsvals(:,2),obj.plotAUC_m.xlsvals(:,1),'filled');
            xlim(obj.plotAUC.Eh.epochtime(1:2));
            ylim([0 1]);
            if(names)
                for ch = 1:numel(obj.plotAUC_m.channels)
                    text(obj.plotAUC_m.xlsvals(ch,2)+0.02,obj.plotAUC_m.xlsvals(ch,1)-0.02,num2str(obj.plotAUC_m.channels(ch)));
                end
            end
        end
        function AUCPlotMean(obj)
            %vykresli graf s prumerem auc krivek, pro vybrane kategorie
            ikat = find(obj.plotAUC.katplot); % index dvojice kategorii v legendkomb
            [l,k] = find(obj.plotAUC.setup.legendkomb==ikat);
            channels = obj.plotAUC_m.channels;
            aucch = obj.plotAUC.aucdata(channels(1)).AUC{k(1),l(1)};  %docasna hodnota pro prvni kanal v poradi
            X = linspace(obj.plotAUC.time(1),obj.plotAUC.time(2),size(aucch,1));
            figure('Name','AUCPlotMean');
            AUC = zeros(numel(channels),size(aucch,1)); %hodnoty pro vsechy kanaly, ktere budu prumerovat
            for ch = 1:numel(channels)                
                aucch  = obj.plotAUC.aucdata(channels(ch)).AUC{k(1),l(1)};                 
                plot(X,aucch(:,1)); 
                AUC(ch,:) = aucch(:,1);
                hold on;
            end
            M = mean(AUC,1);
            E = std(AUC,[],1)/sqrt(size(AUC,1));            
            ciplot(M+E, M-E, X,[1 0.8 0.8]);
            plot(X,M,'-','LineWidth',3,'Color',[1 0 0]);
            ylim([0 1]);
            line([X(1) X(end)],[.5 .5]);             
        end
    end
    methods (Static,Access = public)        
        function W = Wilcox2D(A,B,print,fdr,msg,RjEpChA,RjEpChB,paired)
            %Wilcox2D(A,B,print,fdr,msg,RjEpChA,RjEpChB)
            %srovna dve 3D matice proti sobe, ohledne hodnot v poslednim rozmeru
            %A musi mit oba prvni rozmery > rozmery B, 
            %B muze mit jeden nebo oba prvni rozmer = 1 - pak se porovnava se vsemi hodnotami v A
            %pokud fdr=1,2 nebo prazne, provadi fdr korekci
            %pokud print = 0 nebo prazne, netiskne nic
            if ~exist('print','var'), print = 0; end
            if ~exist('fdr','var') || isempty(fdr), fdr = 1; end %min striktni je default           
            if ~exist('msg','var') || isempty(msg), msg = ''; end
            if size(A,3)==1 %if there are only two dimension of A
                A = reshape(A,size(A,1),1,[]); %assume that samples are in the second dimension and make it the third dimesions
            end            
            if ~exist('RjEpChA','var') || isempty(RjEpChA), RjEpChA = false(size(A,2),size(A,3)); end
            if ~exist('RjEpChB','var') || isempty(RjEpChB), RjEpChB = false(size(B,2),size(B,3)); end
            if ~exist('paired','var')
                paired = iff(numel(B)==1,1,0);  %automatically paired, when B is single value
            end %pokud se ma pouzit parovy neparametricky test, defaulte ne            
            W = zeros(size(A,1),size(A,2));
            
            if print, fprintf(['Wilcox Test 2D - ' msg ' (' num2str(size(A)) '): ']); end
            for j = 1:size(A,1) % napr cas
                if print && mod(j,50)==0, fprintf('%d ', j); end %tisknu jen cele padesatky
                for k = 1:size(A,2) %napr kanaly   
                   if paired %pri parovem testu musim porovnavat stejny kanal, takze musi vyradit epochy parove
                       RjEpChA_k = RjEpChA(k,:) | RjEpChB(k,:); %binarni OR
                       RjEpChB_k = iff(numel(RjEpChB)>1,RjEpChA(k,:) | RjEpChB(k,:),0);
                   else
                       RjEpChA_k = RjEpChA(k,:);
                       RjEpChB_k = RjEpChB(k,:);
                   end    
                   aa = squeeze (A(j,k,~RjEpChA_k)); %jen nevyrazene epochy 
                   bb = squeeze (B( min(j,size(B,1)) , min(k,size(B,2)) , ~RjEpChB_k )); %jen nevyrazene epochy
                   if numel(aa) >= 2 
                      if paired %the signed rank with a single value bb is one-sample test against this value
                        W(j,k) = signrank(aa,bb); %  Wilcoxon signed rank test  paired, two-sided , Statistics and Machine Learning Toolbox  
                      elseif numel(bb) >= 2  %randsum need two sample vectorss
                        W(j,k) = ranksum(aa,bb); % Wilcoxon rank sum test, non-paired, Statistics and Machine Learning Toolbox
                      else
                          W(j,k) = 1; %if only one bb for nonpaired,, no stats can be computed
                      end
                   else
                      W(j,k) = 1; %if only one aa value, no stats can be computed
                   end
                end
            end
            if fdr > 0
                if fdr == 2, method = 'dep'; else method = 'pdep'; end %#ok<SEPEX>
                [~, ~, adj_p]=fdr_bh(W,0.05,method,'no'); %dep je striktnejsi nez pdep 
                %optionally log the uncorredted and corectec values
%                 logcell = cell(size(W,1)*2+2,size(W,2)+5);
%                 logcell{1,1}=method;  logcell{1,2}=msg; logcell{1,5}=datestr(now, 'yyyy-mm-dd_HH-MM-SS'); 
%                 logcell(2:size(W,1)+1,1:size(W,2)) = num2cell(W);
%                 logcell(size(W,1)+2,1) = {'adj_p'};
%                 logcell(size(W,1)+3:end,1:size(W,2)) = num2cell(adj_p);
%                 if(size(logcell,2)>256), logcell = logcell(:,1:256); end %xls can export more columns
%                 xlswrite(['logs\Wilcox2D_' msg '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.xls'],logcell); %zapisu do xls tabulky
                W = adj_p; %prepisu puvodni hodnoty korigovanymi podle FDR
            else
                if print, fprintf('no fdr ...'); end
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
                    W = W'; %vysledne rozmery channel x time(samples)
                end                
                pulsirka = ceil(oknosirka/2); %ktery sloupec se povazuje za pulku okna - tam se hodnota ulozi
                W2 = zeros(size(W,1),size(W,2)); %musim udelat kopii, jinak si prepisuju hodnoty ze kterych pak pocitam
                for sloupec = 1:size(W,2)  %pro vsechny samply
                    iW = max([1 sloupec-pulsirka+1]) : min([size(W,2) sloupec-pulsirka+oknosirka]); %1234 pro sloupec 1, 234567 pro sloupec 4 atd
                    switch funkce
                        case 'min'
                            W2(:,sloupec)=min(W(:,iW),[],2);
                        case 'max'
                            W2(:,sloupec)=max(W(:,iW),[],2); % pro vsechny kanaly najednou
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
        function [AUC,X,Y] = ROCKrivka(epochData,dd,katnames,kresli)
            %11/2018 - spocita a vykresli roc krivku
            %epochData - cellarray z CiEEGData
            %dd - eegdata pro jeden channel a jeden sample
            %sample - cislo sample
            %katnames - jmena dvou kategorii, positive a negative
            if ~exist('kresli','var'), kresli = 0; end
            %labels = epochData(:,1);
            labels = CStat.epochData2Labels(epochData,katnames);
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
        function [labels]= epochData2Labels(epochData,katnames)
            %pokud je v katnames vice kategorii, musim k tomu uzpusobit i katnames 
            labels = epochData(:,1);
            if iscell(katnames) && strcmp(katnames{1}(1:2),'{[') %pokud se jedna o dvojici kategorii
                for k = 1:numel(katnames)
                    for l = 1:numel(labels)
                        if strfind(katnames{k},labels{l})
                           labels{l} = katnames{k};
                        end
                    end
                end
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
        function [timeB,timeK]=StatDiffStart(channels,Wp,kategories,plevel,maintain)  
            %vraci casy zacatku signifikantnich rozdilu vuci baseline a kategorii vuci sobe
            %nezohlednuje smer rozdilu, signum, jako ktere se pouziva treba v CiEEGData.SelChannelStat
            if ~exist('kategories','var'), kategories = Wp.kats; end
            if ~exist('plevel','var'), plevel = 0.05; end
            if ~exist('maintain','var'), maintain = 1; end %for how many samples should the significance stay to be detected. 1=1/fs
                %Bastin 2013jneurosci used maintain=8. But we use already the sliding window in CiEEGStat.WilcoxCat. 
            sigvector = true(1,maintain);
            timeB = NaN(numel(channels),numel(Wp.kats)); % %casy rozdilu vuci baseline
            timeK = NaN(numel(channels),numel(Wp.kats),numel(Wp.kats)); % casy rozdilu mezi kat
            Tr = linspace(Wp.baseline(2),Wp.epochtime(2),size(Wp.D2,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            for k = 1:numel(kategories) %pro vsechny zadane kategorie
                %index where this kategory number is in the list of categories for this stat
                %ik = find(Wp.kats==kategories(k)); %this does not work for cell arrays
                for ikats = 1:numel(Wp.kats)
                    if cellval(Wp.kats,ikats)==cellval(kategories,k)
                        ik = ikats;
                        break;
                    end
                end                
                iWp = Wp.WpKatBaseline{ik,1}(:,channels)  <= plevel; 
                for ch = 1:numel(channels)
                    iWpfirst = strfind(iWp(:,ch)',sigvector); %index of occurences of sigvector - 8 significant differences in sequence
                    if ~isempty(iWpfirst)
                        timeB(ch,k) = Tr(iWpfirst(1)); %time of start of significance
                    end
                end
                for l = k+1:numel(kategories)
                    %index where this kategory number is in the list of categories for this stat
                    %il = find(Wp.kats==kategories(l)); %this does not work for cell arrays
                    for ikats = 1:numel(Wp.kats)
                        if cellval(Wp.kats,ikats)==cellval(kategories,l)
                            il = ikats;
                            break;
                        end
                    end                    
                    if ~isempty(Wp.WpKat{ik,il}) %pro jistotu, nevim v jakem poradi prijdou kategorie
                        WpKat = Wp.WpKat{ik,il}; 
                    else 
                        WpKat = Wp.WpKat{il,ik}; %jestli neni obsazena ta prvni kombinace, tahle bude
                    end
                    iWp = WpKat(:,channels)  <= plevel;  
                    for ch = 1:numel(channels)
                        iWpfirst = strfind(iWp(:,ch)',sigvector); %index of occurences of sigvector - 8 significant differences in sequence
                        if ~isempty(iWpfirst)
                            timeK(ch,k,l) = Tr(iWpfirst(1)); %time of start of significance
                            timeK(ch,l,k) = Tr(iWpfirst(1)); %time of start of significance - currently WpKat is not directional
                        end
                    end
                    
                end
            end
        end
        function [p,t]=CorrelStat(r,n)
            %returns the significance of correlation coeeficient, r can be a vector or matrix
            %according to MXCohen par 34.7 to be used to test significance of ISCP difference - 2021/07
            %z = 0.5*log((1+r)/(1-r)); %fischer z transform of r - not needed
            t = r.*sqrt((n-2)./(1-r.^2)); %t-statistics
            p = 1-tcdf(abs(t),n-2); %two tailed one-sample t-test            
        end
    end
    methods (Access = private)
        function obj = hybejAUCPlot(obj,~,eventDat)
            switch eventDat.Key
                case {'u','i','o','p'} %jina pismena nez f-l, aby se to nepletlo                       
                    ik = find('uiop'==eventDat.Key); %index 1-4
                    if ik <= numel(obj.plotAUC.katplot) %pokud je dost kontrastu mezi kategoriemi k zobrazeni
                        obj.plotAUC.katplot( ik ) = 1 - obj.plotAUC.katplot(ik);
                    end
                    %obnovim jednu krivku AUC
                    obj.AUCPlot(find(obj.plotAUC.Eh.CH.sortorder==obj.plotAUC.ch));  %#ok<FNDSB>
                case {'c'} %zmeni barvy pasu stderr na nepruhledne, aby se daly kopirovat do corelDRAW / a zpet
                    obj.plotAUC.corelplot = 1 - obj.plotAUC.corelplot;
                    obj.AUCPlot(find(obj.plotAUC.Eh.CH.sortorder==obj.plotAUC.ch));%#ok<FNDSB>
                case {'f','g','h','j','k','l'}                    
                    channels = find(obj.plotAUC.selCh(:,'fghjkl'==eventDat.Key))'; %indexy kanalu se znackou f-l
                    %vytvorim multiple AUC graf:
                    obj.AUCPlotM(channels,find('fghjkl'==eventDat.Key),0); %#ok<FNDSB> %povinne ted uvadim predvybrany kanal
                case {'a'} %chci zobrazit krivky ze vsech kanalu
                    channels = 1:size(obj.plotAUC.selCh,1);
                    obj.AUCPlotM(channels,0,0); %#ok<FNDSB> % 0 znamena vsechy kanaly, povinne ted uvadim predvybrany kanal
                case 'r' %reverse order of categories in the ROC analysis
                    obj.plotAUC.reversekats = 1 - obj.plotAUC.reversekats;
                    obj.AUCPlot(find(obj.plotAUC.Eh.CH.sortorder==obj.plotAUC.ch),obj.plotAUC.Eh);%#ok<FNDSB>
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
                case 'm' %plot mean instead of individual channels
                    obj.plotAUC_m.plotmean = 1  - obj.plotAUC_m.plotmean;                    
                case 'd' %just redraws
                    kresli = 1;
                    selch = obj.plotAUC_m.selch;
                otherwise     
                    kresli = 0;
            end
            if kresli 
                obj.AUCPlotM([],[],selch); %prekresli sumarni plot AUC krivek
                if selch > 0
                    ch = find(obj.plotAUC.Eh.CH.sortorder==obj.plotAUC_m.channels(obj.plotAUC_m.chsort(selch))); %index vybraneho kanalu v CH.sortorder
                    obj.plotAUC.Eh.PlotResponseCh(ch); %#ok<FNDSB> %prekresli graf PlotResponseCh
                end
                figure(obj.plotAUC_m.fh);
            end
        end        
        function hybejPlot3Dclick(obj, ~, ~, ~)
          mousept = get(gca,'currentPoint');
          p1 = mousept(1,:); p2 = mousept(2,:); % souradnice kliknuti v grafu - predni a zadni bod
          displayedChannels = obj.plotAUC.Eh.CH.H.channels(obj.plotAUC_m.channels); % zobrazene kanaly
          coordinates = [displayedChannels.MNI_x; displayedChannels.MNI_y; displayedChannels.MNI_z];    % souradnice zobrazenych kanalu
          closestChannel = findClosestPoint(p1, p2, coordinates, 2);    % najdu kanal nejblize mistu kliknuti
          if closestChannel  % pokud jsem nejaky nasel:
             ch = obj.plotAUC.Eh.CH.channelPlot.plotCh3D.chnsel(closestChannel);
             if isfield(obj.plotAUC.Eh.CH.channelPlot.plotCh3D, 'fh') && isvalid(obj.plotAUC.Eh.CH.channelPlot.plotCh3D.fh)
                obj.plotAUC.Eh.CH.channelPlot.highlightChannel(ch);
             end
             obj.AUCPlotM([],[],find(obj.plotAUC_m.chsort == closestChannel, 1)); % oznacim vybrany kanal v AUC plotu
          else  % pokud se zadny kanal nenasel (kliknuti mimo)
             if isfield(obj.plotAUC.Eh.CH.channelPlot.plotCh3D, 'fh') && isvalid(obj.plotAUC.Eh.CH.channelPlot.plotCh3D.fh)
                obj.plotAUC.Eh.CH.channelPlot.highlightChannel(0);
             end
             %TODO: Zrusit zvyrazneni v AUC plotu
          end
        end    
        function chtodoNum = AUCisEmpty(obj,channels)
            if ~exist('channels','var') || isempty(channels), channels = 1:numel(obj.plotAUC.aucdata); end
            chtodo = false(numel(channels),1);
            for ch = 1:numel(channels) %ch is just index in channels, no real channel number
                kk = 1; ll = 2;
                if obj.plotAUC.reversekats, [l, k] = deal(kk,ll); else,  [k, l] = deal(kk,ll); end %swap both variables content
                if length(obj.plotAUC.aucdata)<channels(ch) || isempty(obj.plotAUC.aucdata(channels(ch)).AUC) || isempty(obj.plotAUC.aucdata(channels(ch)).AUC{k,l})
                    chtodo(ch) = true; 
                end              
            end
            if numel(channels) == 1
                chtodoNum = chtodo(1); %if only one channel, returs boolean
            else
                chtodoNum = channels(chtodo); %otherwise returns real channel numbers 
            end
        end
        
    end
    
end

