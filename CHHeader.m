classdef CHHeader < matlab.mixin.Copyable %je mozne kopirovat pomoci E.copy();
    %CHHEADER Trida na praci s headerem od Jirky Hammera
    %Kamil Vlcek, FGU AVCR, since 2016 04
    
    properties (Access = public)
        H; %data od Jirky Hammera
        %E; % unikatni jmena elektrod, napr RG, I, GC ...
        chgroups; %cell array, kazda bunka jsou cisla kanalu ve skupine
        els;
        RjCh;
        BrainAtlas_zkratky; %tabulka od Martina Tomaska
        filterMatrix; %kopie filterMatrix, vytvari se pri zmene reference
        sortorder; %index serazenych kanalu
        sortedby; %podle ceho jsou kanaly serazeny
        plotCh2D; %udaje o 2D grafu kanalu ChannelPlot2D, hlavne handle
        plotCh3D; %udaje o 3D grafu kanalu ChannelPlot, hlavne handle
        reference; %aby trida vedela, jestli je bipolarni nebo ne        
        classname; %trida v ktere je Header vytvoren
    end
    %#ok<*PROPLC>
    
    events
        FilterChanged
    end
    
    methods (Access = public)
        function obj = CHHeader(H,filename,reference,classname)
            %konstruktor
            %filename se zatim na nic nepouziva
            obj.H = H;     
            if exist('filename','var') && ~isempty(filename)
               obj.H.filename = filename;
            end
              %nefunguje protoze name je napriklad RG12 a ne jen RG, takze neni unikatni pro kazdou elektrodu
%             E = cell(size(H.electrodes,2),1); %#ok<*PROP>
%             for iE= 1:size(H.electrodes,2)
%                 E{iE}=H.electrodes(iE).name;
%             end
%             [U,idx]=unique(cell2mat(E)); 
%             obj.E = cell(numel(U),1);
%             for iE = 1:numel(idx)
%                 obj.E{iE}= E{idx};
%             end
             obj = obj.SelChannels(); 
             obj.sortorder = 1:numel(obj.H.channels); %defaultni sort order
             obj.sortedby = '';
             if exist('reference','var'), obj.reference = reference; else, obj.reference = []; end
             if exist('classname','var'), obj.classname = classname; else, obj.classname = []; end
        end
        
        function [obj, chgroups, els] = ChannelGroups(obj,chnsel,subjname,forcechgroups)
            %vraci skupiny kanalu (cisla vsech channels na elekrode) + cisla nejvyssiho kanalu v na kazde elektrode v poli els
            % subjname - jestli jsou skupiny podle subjname (napr p173, u CHibertMulti) nebo elektrod (napr A)  
            if ~exist('subjname','var') || isempty(chnsel), subjname = iff( strcmp(obj.classname,'CHilbertMulti')); end %defaultne podle classname, pokud neni CHilbertMulti, pouzivaji se jmena elektrod
            if ~exist('chnsel','var') || isempty(chnsel) %pokud nenam zadny vyber kanalu, pouziva se pri load dat
            if ~exist('forcechgroups','var') || isempty(forcechgroups), forcechgroups = 0; end %chci prepocitat chgroups 
                if isempty(obj.chgroups) || forcechgroups               
                    chgroups = obj.getChannelGroups(subjname); %pouziju vlastni funkci, getChannelGroups_kisarg je hrozne pomale
                    els = zeros(1,numel(chgroups)); 
                    for j = 1:numel(chgroups)
                        els(j)=max(chgroups{j});
                    end
                    els = sort(els);
                    obj.chgroups = chgroups; 
                    obj.els = els;
                else
                    chgroups = obj.chgroups; 
                    els = obj.els;
                end
            else %pokud mam definovan vyber kanalu v chsel. Pouzivam z grafu ChannelPlot
                if obj.plotCh3D.allpoints %pokud chci zobrazovat i ostatni kanal jako tecky
                   chgroups = {chnsel,setdiff(obj.H.selCh_H,chnsel)}; %do druhe skupiny dam ostatni kanaly
                elseif subjname
                   chgroups = obj.getChannelGroups(subjname,chnsel);  
                else
                   chgroups = {chnsel}; %pokud mam vyber kanalu, zatim to nechci resit a vsechny v jedne skupine - bez ohledu na elektrody, jako cellarray
                end
                els = obj.els;
            end
        end
        function [obj,els2plot ] = ElsForPlot(obj)
            %vrati cisla nejvyssiho kanalu pro zobrazeni - kdyz je nejaka elektroda moc dlouha, tak ji rozdeli
            els2plot = zeros(1,numel(obj.els));
            e0 = 0; %cislo elektrody z minuleho cyklu
            iEls =  1; %index v els2plot;
            kontaktumax = 15;
            for e = 1:numel(obj.els)
                while obj.els(e)-e0 > kontaktumax %pokud je nejaka elektroda moc dlouha, rozdelim ji do zobrazeni po padesati kouscich
                    els2plot(iEls) = e0 + kontaktumax;
                    e0 = e0+kontaktumax;
                    iEls = iEls+1;
                    if e==numel(obj.els)
                        els2plot = [els2plot, obj.els(e)]; %#ok<AGROW> %11.7.2019 - pokud jsem uz u posledni elektrody, musi jeste pridat zvyvajici kontakty 
                    end
                end
                if obj.els(e)-e0 > 3 %pokud je elektroda prilis kratka, nebudu ji zobrazovat samostatne
                    els2plot(iEls) = obj.els(e);
                    e0= obj.els(e);
                    iEls = iEls+1;
                end
            end 
            els2plot(els2plot==0)=[];
        end
        function obj = RejectChannels(obj,RjCh)
            %ulozi cisla vyrazenych kanalu - kvuli pocitani bipolarni reference 
            obj.RjCh = RjCh;            
        end
        function chs = GetChOK(obj)
            %vraci cisla kanalu, ktera jsou SEEG a nejsou vyrazena
            chs = setdiff(obj.H.selCh_H,obj.RjCh);
        end
        function chs = ChannelsN(obj)
            %vraci celkovy pocet kanalu
            chs = length(obj.H.channels);
        end
        function [ch,name,MNI,brainAtlas] = ChannelNameFind(obj,name)
            % najde kanal podle jmena kontaktu, napr 'R5', name musi byt na zacatku jmena, 
            % ale nemusi byt cele
            % vraci jeho cislo, cele jmeno, mni souradnice a pojmenovani podle brain atlasu 
            MNI = [];
            brainAtlas = '';
            for ch = 1:size(obj.H.channels,2)
                if 1==strfind(obj.H.channels(1,ch).name,name)
                    if isfield(obj.H.channels,'MNI_x') %pokud mni souradnice existuji
                        MNI = [obj.H.channels(ch).MNI_x obj.H.channels(ch).MNI_y obj.H.channels(ch).MNI_z];
                    end
                    if isfield(obj.H.channels,'ass_brainAtlas')
                        brainAtlas = obj.H.channels(ch).ass_brainAtlas;
                    end
                    name = obj.H.channels(1,ch).name;
                    return;
                end
            end
            ch = 0;
        end
        function [XYZ,obj] = ChannelPlot(obj,pohled,labels,chnvals,chnsel,selch,roi,popis,rangeZ)
            %zobrazi 3D obrazek elektrod v MNI prostoru. Obrazek ma rozmery podle rozmeru mozku
            %pohled muze urcti smer pohledu s-sagital,c-coronal,h-horizontal
            %chnsel jsou cisla kanalu, pokud chci jen jejich vyber - musi byt stejny pocet jako chvals (hodnoty k vykresleni)
            %selch je jedno zvyraznene cislo kanalu - index v poli chnsel
            %roi je zvyraznena krychlova oblast [ x y z edge]
            %popis je text k zobrazeni na obrazku
%             if ~exist('pohled','var') || isempty(pohled), pohled = ''; end            
            
            params = {'pohled','labels','chnvals','chnsel','selch','roi','popis','rangeZ'}; %zkusim hromadne zpracovani parametru touhle nedoporucovanou metodou
            iSEEG = contains({obj.H.channels.signalType},'SEEG'); %index kanalu s EEG signalem
            for p=1:numel(params) %parametry, ktere se ukladaji do obj.plotCh3D
                if ~exist(params{p},'var') || eval(['isempty(' params{p} ')']) %pokud neni vstupni promenna nebo je prazdna
                    if isfield(obj.plotCh3D,params{p}) %pokud ale existuje ulozena hodnota
                        eval([ params{p} ' = obj.plotCh3D.' params{p} ';']); %tak ji pouziju
                    else 
                        switch params{p}
                            case 'pohled'
                                pohled =  ''; %defaultni pohled nedefinovany = horizontal
                            case 'labels'
                                labels = 0;  %defaultne se nezobrazuji neurology labels, ale jmena kanalu
                            case 'chnvals'
                                chnvals = zeros(1, numel(obj.H.channels(iSEEG))); %default same nuly
                            case 'chnsel'
                                chnsel = 1:numel(obj.H.channels(iSEEG)) ; %default vsechny kanaly
                            case 'selch'
                                selch = []; %default vsechny kanaly
                            case 'roi'
                                roi = []; %default zadne
                            case 'popis'
                                popis = ''; %default zadny text
                            case 'rangeZ'
                                rangeZ = [];                            
                        end
                        eval(['obj.plotCh3D.' params{p} ' = ' params{p} ';']); %nastavim ulozenou hodnotu na default
                    end
                else
                    eval(['obj.plotCh3D.' params{p} ' = ' params{p} ';']); %podle vstupni promenne zmeni ulozenou hodnotu
                end
            end           
            if ~isfield(obj.plotCh3D,'allpoints'), obj.plotCh3D.allpoints = 0; end
            if ~isfield(obj.plotCh3D,'zoom'), obj.plotCh3D.zoom = 0; end            
            if ~isfield(obj.plotCh3D,'reorder'), obj.plotCh3D.reorder = 0; end   %defaultne se neprerazuji kanaly podle velikosti
            if ~isfield(obj.plotCh3D,'lines'), obj.plotCh3D.lines = 0; end   %defaultne se nespojuji pacienti spojnicemi            
            assert(numel(chnvals) == numel(chnsel), 'unequal size of chnvals and chnsel');
            nblocks = numel(chnvals); %pocet barev bude odpovidat poctu kanalu
            cmap = parula(nblocks+1); %+1 protoze hodnoty se budou zaokrouhlovat nahoru nebo dolu
            reverse = 0; %jestli obratit barevnou skalu a velikosti
            if isempty(rangeZ)
                rangeZ = [min(chnvals) max(chnvals)];                 
            elseif rangeZ(1) > rangeZ(2) %pokud dam minmax v obrazenem poradi, barvy i velikosti taky v obracenem poradi
                reverse = 1;
                rangeZ = flip(rangeZ);            
                cmap = flip(cmap,1);
            end
            
            chnvalsN = chnvals - rangeZ(1); %odectu mimimum
            chnvalsN = chnvalsN/diff(rangeZ); % normalization  - podelim maximem - hodnoty jsou [0;1]          
            chnvalsN(isnan(chnvalsN)) = 0; % in case of all zeros, nan nahradim 0
            chnvalsN(chnvalsN<0) = 0; chnvalsN(chnvalsN>1) = 1; %omezim rozsah na [0;1];            
            clrs = cmap(round(nblocks*chnvalsN)+1, :); % color values, prevedu na rozsah 1-nblocks a priradim barvy
            sizes = 20+200*iff(reverse,1-chnvalsN,chnvalsN); %velikosti kulicek 
            %if reverse, sizes = flip(sizes); end
            if isfield(obj.H.channels,'MNI_x')
                if isfield(obj.plotCh3D,'fh') && ishandle(obj.plotCh3D.fh)
                    figure(obj.plotCh3D.fh); %pouziju uz vytvoreny graf
                    clf(obj.plotCh3D.fh); %graf vycistim                     
                else
                    obj.plotCh3D.fh = figure('Name','ChannelPlot 3D in MNI');                     
                    obj.plotCh3D.isColormapReversed = 0;
                end          
                               
                [obj,chgroups] = obj.ChannelGroups(chnsel,obj.plotCh3D.lines); %rozdeli kanaly po elektrodach do skupin. 
                 %Pokud chnsel, jsou vsecny v jedne skupine. Ale pokud obj.plotCh3D.allpoints, ve druhe skupine jsou ostatni kanaly
                
                %objekt se dobre uklada i pri poradi return values XYZ,obj
                XYZ = struct('X',0,'Y',0,'Z',0);
                for chg = 1:size(chgroups,2) 
                    group = chgroups{chg};                     
                    X = [obj.H.channels(group).MNI_x];
                    Y = [obj.H.channels(group).MNI_y];
                    Z = [obj.H.channels(group).MNI_z];                                             
                    
                    linestyle = iff(numel(chgroups)>1 && ~obj.plotCh3D.allpoints,'-','.'); %cara bude jina pokud je pouzite chnsel
                    plot3(X,Y,Z,linestyle,'LineWidth',2);
                    if chg==1, hold on; end                         
                    if labels == 1 %neurology labels
                        names = {obj.H.channels(group).neurologyLabel};
                    elseif labels == 2 %channel names
                        names = {obj.H.channels(group).name};
                    elseif labels == 3 %pacient names
                        names = {obj.H.channels(group).name};
                        names = cellstr(extractBefore(names,' ')); %vsechno pred mezerou - pro CHilbertMulti
                    end
                    iZ = mod(1:numel(Z), 2); iZ(iZ == 0) = -1;                    
                    if ~(chg>1 && obj.plotCh3D.allpoints) %prvni skupiny do barevnych kulicek vzdy; 
                        %druhou skupinu chci jen pokud zobrazuju vsechny (chnsel je prazdne) nebo pokud nejsou v druhe skupine ostatni kanaly
                        XYZ(chg) = struct('X',X,'Y',Y,'Z',Z); %export pro scatter3 nize, ktery zobrazi ruzne velke a barevne kulicky
                        if labels>0
                            text(X+abs(iZ)*0.5,Y,Z+iZ*0.5,names,'FontSize', 7);
                        end
                    end
                end
                % Plot with different colors and sizes based on chnvals
                if isempty(chnvals)  %indexy vsech kanalu, nemam zadne hodnoty k vykresleni
                    isizes = 1:obj.H.channels;
                elseif obj.plotCh3D.allpoints  %zobrazuju pozice vsech kanalu jako tecek (dve skupiny kanalu v chgroups - barevne kulicky + tecky)
                    isizes = find(chnsel==[chgroups{1}]); %indexy v poli chnsel pro pouziti v poli sizes, find pracuje i hromadne
                else %indexy vsech kanalu ve vsech skupinach
                    isizes = find(chnsel==[chgroups{:}]); %indexy v poli chnsel pro pouziti v poli sizes, find pracuje i hromadne
                end   
                X = [XYZ.X]; Y = [XYZ.Y]; Z = [XYZ.Z]; %souradnice pres vsechny pole struct XYZ
                if obj.plotCh3D.reorder %pokud chci seradi body podle velikosti, tak aby v prislusnem pohledu byly nejvetsi v popredi
                    switch obj.plotCh3D.pohled
                        case 'h'
                            Z = sortBlikeA(sizes,Z); %nejvetsi hodnoty na nejvyssich souradnicich Z
                        case 'c'
                            Y = sortBlikeA(-sizes,Y); %nejvetsi hodnoty na nejnizsich souradnicich y
                        case 's'
                            X = sortBlikeA(sizes,X);
                    end
                    annotation('textbox', [.6 0.15 .2 .1], 'String', 'REORDERED', 'EdgeColor', 'none');
                end
                scatter3(X,Y,Z,sizes(isizes),clrs(isizes,:),'filled'); %ruzne velke a barevne krouzky vsech kanalu najednou
               
                if ~isempty(selch) && selch>0
                    scatter3(XYZ.X(selch),XYZ.Y(selch),XYZ.Z(selch),max(sizes),[0 0 0]);
                end
                
                if ~isempty(roi) && numel(roi)>=4 %[x y z edge]
                    for r = 1:size(roi,1)                    
                        plotcube([roi(r,4) roi(r,4) roi(r,4)], roi(r,1:3),0,[0 0 0]); %ROI jako kostka
                        text(roi(r,1)+roi(r,4)/2,  roi(r,2)+roi(r,4)/2 , roi(r,3)+roi(r,4)/2, num2str(r),'FontSize', 10, 'Color',[1 0 0]);
                    end
                end
                
                xlabel('MNI X'); %levoprava souradnice
                ylabel('MNI Y'); %predozadni souradnice
                zlabel('MNI Z'); %hornodolni
                
                switch pohled
                    case ''
                        if isfield(obj.plotCh3D,'view')
                            view(obj.plotCh3D.view); 
                        else
                            view([0 0 1]); %shora - pokud neni zadny ulozeny
                        end
                    case 's' %sagital = levoprava
                        view([1 0 0]); %zleva
                    case 'c' %coronal = predozadni
                        view([0 -1 0]); %zepredu
                    case 'h' %horizontal = hornodolni   
                        view([0 0 1]); %shora
                end
                if obj.plotCh3D.zoom == 0
                    axis([-75 75 -120 80 -75 85]); %zhruba velikost mozku        
                else 
                    axis([ min([XYZ.X]) max([XYZ.X]) min([XYZ.Y]) max([XYZ.Y]) min([XYZ.Z]) max([XYZ.Z]) ] );
                end
                text(-70,0,0,'LEVA');        
                text(70,0,0,'PRAVA');   
                text(0,65,0,'VPREDU');        
                text(0,-115,0,'VZADU');                                 
               
                if isfield(obj.plotCh3D,'boundary') && obj.plotCh3D.boundary
                    obj.Plot3DBoundary();
                else
                    load('GMSurfaceMesh.mat'); %seda hmota v MNI
                    scatter3(GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3),'.','MarkerEdgeAlpha',.2);
                end
                
                if(max(chnvals)>0)
                    colorbar;
                    if (reverse && ~obj.plotCh3D.isColormapReversed) || (~reverse && obj.plotCh3D.isColormapReversed)
                        oldcmap = colormap;
                        colormap( flipud(oldcmap) ); %prevratim colomapu, jinak se zobrazuje defaultni poradi, bez ohledu na moje prehozeni                  
                        obj.plotCh3D.isColormapReversed = 1 - obj.plotCh3D.isColormapReversed;                     
                    end
                    caxis(rangeZ); 
                end %barevna skala, jen pokud jsou ruzne hodnoty kanalu
                if obj.plotCh3D.zoom < 2, axis equal;  end %maximalni zoom je bez stejnych os
                title(popis);
                if isfield(obj.plotCh3D,'background') && obj.plotCh3D.background==0
                    set(gca,'color','none'); %zadne bile pozadi, pak ani v corelu
                end
                %rozhybani obrazku            
                set(obj.plotCh3D.fh,'KeyPressFcn',@obj.hybejPlot3D);
            else
                disp('No MNI data');
            end
        end             
        function ChannelPlot2D(obj,chsel,plotRCh,plotChH,label)
            %vstupni promenne
            %plotRCh - cela struktura plotRCh z CiEEGData
            if ~exist('chsel','var') %promenna na jeden cerveny kanal
                if isfield(obj.plotCh2D,'chsel')
                    chsel = obj.plotCh2D.chsel;
                else
                    chsel = 1;  
                    obj.plotCh2D.chsel = 1;
                end
            else
                obj.plotCh2D.chsel = chsel;
            end
            chselo = obj.sortorder(chsel); %jeden kanal, ktery je zobrazeny v PlotResponseCh
                                           %chsel je index v sortorder, chselo je skutecne cislo kanalu
            if ~exist('plotRCh','var') 
                if isfield(obj.plotCh2D,'selCh') %promenna na vic oznacenych kanalu f-l, podle obj.PlotRCh.SelCh
                    selCh = obj.plotCh2D.selCh;
                else
                    selCh = []; 
                    obj.plotCh2D.selCh = [];
                end
                if isfield(obj.plotCh2D,'selChNames') %pojmenovani vyberu kanalu pomoci f-l
                    selChNames = obj.plotCh2D.selChNames;
                else
                    selChNames = cell(1,6); 
                    obj.plotCh2D.selChNames = selChNames;
                end
            else
                if isfield(plotRCh,'selCh')
                    obj.plotCh2D.selCh = plotRCh.selCh;
                    selCh = plotRCh.selCh; 
                else
                    selCh = []; 
                end
                if isfield(plotRCh,'selChNames')
                    obj.plotCh2D.selChNames = plotRCh.selChNames;
                    selChNames = plotRCh.selChNames;
                else
                    selChNames = cell(1,6);
                    obj.plotCh2D.selChNames = selChNames;                    
                end
            end
            if exist('plotChH','var')  %handle na funkci z CiEEGData @obj.PlotResponseCh
                obj.plotCh2D.plotChH = plotChH;
            end
            if ~isfield(obj.plotCh2D,'marks')  %handle na funkci z CiEEGData @obj.PlotResponseCh
                obj.plotCh2D.marks = [1 1 1 1 1 1]; %ktere znacky fghjjkl se maji zobrazovat
            end
            if ~exist('label','var') %promenna z CM oznacujici nejaky label celeho souboru 
                if isfield(obj.plotCh2D,'label')
                    label = obj.plotCh2D.label;
                else
                    label = ''; 
                    obj.plotCh2D.label = '';
                end
            else
                obj.plotCh2D.label = label;
            end
            if ~isfield(obj.plotCh2D,'chseltop'), obj.plotCh2D.chseltop = 0; end %jestli se ma vybrany kanal zobrazovat na popredi ostatnych  - zlute kolecko
            if ~isfield(obj.plotCh2D,'names'), obj.plotCh2D.names = 1; end %jestli se maji vypisovat jmena kanalu
            if ~isfield(obj.plotCh2D,'lines'), obj.plotCh2D.lines=1; end %defaltne se kresli cary mezi kanaly jedne elektrody
            if ~isfield(obj.plotCh2D,'transparent'), obj.plotCh2D.transparent=0; end %defaltne se kresli body nepruhledne
            if ~isfield(obj.plotCh2D,'chshow'), obj.plotCh2D.chshow = 1:numel(obj.H.channels); end %defaltne se kresli body nepruhledne
            if ~isfield(obj.plotCh2D,'ch_displayed'), obj.plotCh2D.ch_displayed=obj.plotCh2D.chshow; end %defaltne jsou zobrazeny vsechny vybrane kanaly (podle filtru)
            if ~isfield(obj.plotCh2D,'chshowstr'), obj.plotCh2D.chshowstr = ''; end   %defaultne bez filtrovani
            %------------------------- vytvoreni figure -----------------------------------
            x = [obj.H.channels(:).MNI_x];
            y = [obj.H.channels(:).MNI_y];
            z = [obj.H.channels(:).MNI_z];            
            load('GMSurfaceMesh.mat'); %seda hmota v MNI
            if isfield(obj.plotCh2D,'boundary') && obj.plotCh2D.boundary && ~isfield(obj.plotCh2D,'BrainBoundaryXY') %trva docela dlouho nez se to spocita
                obj.plotCh2D.BrainBoundaryXY = boundary(GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2)); %vnejsi hranice mozku
                obj.plotCh2D.BrainBoundaryYZ = boundary(GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3));
            end
            
            size_ch = 12; %velikosti krouzko oznacujicich kanaly
            size_selCh = 50; %7;
            x_text = -100;
            if isfield(obj.plotCh2D,'fh') && ishandle(obj.plotCh2D.fh)
                figure(obj.plotCh2D.fh); %pouziju uz vytvoreny graf
                clf(obj.plotCh2D.fh); %graf vycistim
            else
                obj.plotCh2D.fh = figure('Name',['ChannelPlot2D - ' label]);                     
            end            
                   
            if isfield(obj.plotCh2D,'chshow') && numel(obj.plotCh2D.chshow) < numel(obj.H.channels) %pokud chci zobrazovat jen cast kanalu podle chshow
                els = obj.plotCh2D.chshow;
                els0 = obj.plotCh2D.chshow; %nebudu resit zacatky a konec elektrod
                chshow = obj.plotCh2D.chshow;
            else
                els = obj.els; %konce elektrod/pacientu u CM
                els0 = [1 obj.els(1:end-1)+1]; %zacatky elektrod/pacientu u CM
                chshow = 1:numel(obj.H.channels); 
            end
            
            subplot(1,2,1);
            % ----------------- axialni plot   ---------------------------           
            if isfield(obj.plotCh2D,'boundary') && obj.plotCh2D.boundary
                %defaultne budu vykreslovat scatter, ale kvuli kopirovani se bude hodit i jen boundary
                plot(GMSurfaceMesh.node(obj.plotCh2D.BrainBoundaryXY,1),GMSurfaceMesh.node(obj.plotCh2D.BrainBoundaryXY,2));
            else
                scatter(GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2),'.','MarkerEdgeAlpha',.1); %seda hmota normalizovaneho mozku
            end           
            hold on;             
            
            for ie = 1:numel(els) 
                plotstyle = iff(obj.plotCh2D.lines,'-o','o'); %,'ok' pro cernobile
                if obj.plotCh2D.lines >= 0 %-1 znamena, ze se nemaji zobrazovat neoznacene kanaly pomoci fghjkl, cili tady se nekresli nic
                    plot(x(els0(ie):els(ie)),y(els0(ie):els(ie)),plotstyle); %plot kontaktu jedne elektrody
                end
                if obj.plotCh2D.names > 0
                for ch = els0(ie):els(ie)
                        if obj.plotCh2D.lines >= 0 || isempty(selCh) || any(selCh(ch,logical(obj.plotCh2D.marks))) %pokud je kanal oznacen jednou ze zobrazenych znacek
                            if obj.plotCh2D.names == 2
                                th = text(x(ch),y(ch),obj.H.channels(ch).name); %jmeno kanalu
                            else
                                th = text(x(ch),y(ch),num2str(ch)); %cislo kazdeho kanalu
                            end
                            th.FontSize = 8;
                        end
                end
                end                              
            end
            if ~isempty(chsel) %pokud je vybrany nejaky kanal                
                h_selection = plot(x(chselo),y(chselo),'o','MarkerSize',size_ch,'MarkerEdgeColor','y','MarkerFaceColor','y'); 
                chstr = iff(isempty(obj.sortedby),num2str(chsel), [ num2str(obj.sortorder(chsel)) '(' obj.sortedby  num2str(chsel) ')' ]);
                title( [ 'channel ' chstr ]);
                
            end
            if ~isempty(selCh) %vybery kanalu fghjkl
                barvy = 'gbrcmk';
                ch_displayed = cell(1,6);
                for m = size(selCh,2):-1:1 %jednu znacku za druhou
                   if  obj.plotCh2D.marks(m) %pokud se ma znacka zobrazovat
                       ch = find(selCh(:,m)); %seznam cisel vybranych kanalu pro danou znacku
                       ch = intersect(chshow,ch); 
                       %plot(x(ch),y(ch),'o','MarkerSize',size_selCh,'MarkerEdgeColor',barvy(m),'MarkerFaceColor',barvy(m));
                       sh = scatter(x(ch),y(ch),size_selCh,barvy(m),'filled');
                       if obj.plotCh2D.transparent, alpha(sh,.5); end %volitelne pridani pruhlednosti
                       ch_displayed{m} = ch';
                   end
                end
                if obj.plotCh2D.lines < 0
                    obj.plotCh2D.ch_displayed = unique(cell2mat(ch_displayed))'; %ktere kanaly jsou opravu zobrazeny, do radku
                else
                    obj.plotCh2D.ch_displayed = chshow;
                end
            end
            if obj.plotCh2D.chseltop, uistack(h_selection, 'top'); end
            text(-70,70,'LEVA');
            text(55,70,'PRAVA');  
            axis equal;  
            if isfield(obj.plotCh2D,'grid') && obj.plotCh2D.grid==1
                grid on;
            end
                
            set(gca, 'XTick',-70:10:70); %xticks(-70:10:70); %xtics jsou az od 2016b
            set(gca, 'YTick',-100:10:70); %yticks(-100:10:70); %ytics jsou az od 2016b
            xlabel('MNI X'); %levoprava souradnice
            ylabel('MNI Y'); %predozadni souradnice
            if isfield(obj.plotCh2D,'background') && obj.plotCh2D.background==0
                set(gca,'color','none'); %zadne bile pozadi, pak ani v corelu
            end
           
            subplot(1,2,2);
            % ----------------- sagitalni plot ---------------------------
            if isfield(obj.plotCh2D,'boundary') && obj.plotCh2D.boundary
                plot(GMSurfaceMesh.node(obj.plotCh2D.BrainBoundaryYZ,2),GMSurfaceMesh.node(obj.plotCh2D.BrainBoundaryYZ,3));
            else
                scatter(GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3),'.','MarkerEdgeAlpha',.1);   %seda hmota normalizovaneho mozku
            end
            hold on;     
            for ie = 1:numel(els)  
                plotstyle = iff(obj.plotCh2D.lines,'-o','o');
                if obj.plotCh2D.lines >= 0 %-1 znamena, ze se nemaji zobrazovat neoznacene kanaly pomoci fghjkl, cili tady se nekresli nic
                    plot(y(els0(ie):els(ie)),z(els0(ie):els(ie)),plotstyle); %plot kontaktu jedne elektrody
                end
                if obj.plotCh2D.names
                for ch = els0(ie):els(ie)
                    if obj.plotCh2D.lines >= 0 || isempty(selCh) || any(selCh(ch,logical(obj.plotCh2D.marks))) %pokud je kanal oznacen jednou ze zobrazenych znacek
                        if obj.plotCh2D.names == 2
                            th = text(y(ch),z(ch),obj.H.channels(ch).name); %jmeno kanalu
                        else
                            th = text(y(ch),z(ch),num2str(ch)); %cislo kazdeho kanalu
                        end                    
                        th.FontSize = 8;
                    end
                end                
                end
            end  
            if ~isempty(chsel) %pokud je vybrany nejaky kanal
                h_selection = plot(y(chselo),z(chselo),'o','MarkerSize',size_ch,'MarkerEdgeColor','y','MarkerFaceColor','y'); 
                
                text(x_text,110,[ obj.H.channels(1,chselo).name]);
                text(x_text,100,[ obj.H.channels(1,chselo).neurologyLabel ',' obj.H.channels(1,chselo).ass_brainAtlas]);
                if  isfield(obj.H.channels,'MNI_x') %vypisu MNI souradnice
                    text(x_text,90,[ 'MNI: ' num2str(round(obj.H.channels(1,chselo).MNI_x)) ', ' num2str(round(obj.H.channels(1,chselo).MNI_y )) ', ' num2str(round(obj.H.channels(1,chselo).MNI_z))]);
                else
                    text(x_text,90,'no MNI');
                end                
            end
            if ~isempty(selCh) %hromadne vybrane kanaly, zobrazne cernym koleckem                
                barvy = 'gbrcmk';
                klavesy = 'fghjkl'; %abych mohl vypsat primo nazvy klaves vedle hvezdicky podle selCh
                for m = size(selCh,2):-1:1 %jednu znacku za druhou - naposled ty prvni aby byly nahore
                    if  obj.plotCh2D.marks(m) %pokud se ma znacka zobrazovat
                       ch = find(selCh(:,m)); %seznam cisel vybranych kanalu pro danou znacku
                       ch = intersect(chshow,ch); 
                       if ~isempty(ch) %pokud jsou takove nejake vybrane kanaly
                           %plot(y(ch),z(ch),'o','MarkerSize',size_selCh,'MarkerEdgeColor',barvy(m),'MarkerFaceColor',barvy(m));
                           sh = scatter(y(ch),z(ch),size_selCh,barvy(m),'filled');
                           if obj.plotCh2D.transparent, alpha(sh,.5); end %volitelne pridani pruhlednosti
                           
                           th = text(x_text+m*10,-90,klavesy(m), 'FontSize', 15,'Color',barvy(m)); %legenda k barvam kanalu dole pod mozkem
                           th.BackgroundColor = [.6 .6 .6];
                           if ~isempty(selChNames) && ~isempty(selChNames{m})
                             text(x_text+70,-60-m*7,cell2str(selChNames{m}), 'FontSize', 9,'Color',barvy(m)); %popisy znacek f-l                           
                           end
                       end
                    end
                end
                if any(selCh(chselo,:),2)==1 %pokud je aktualni kanal jeden z vybranych                
                    klavesy = 'fghjkl'; %abych mohl vypsat primo nazvy klaves vedle hvezdicky podle selCh
                    text(x_text,80,['*' klavesy(logical(selCh(chselo,:)))], 'FontSize', 12,'Color','red');
                end
                if ~isempty(label)
                    text(x_text,-75,strrep(label,'_','\_'), 'FontSize', 10,'Color','blue' );
                end
                if isfield(obj.plotCh2D,'chshowstr') && ~isempty(obj.plotCh2D.chshowstr)
                    text(0,-90,['chshow:' obj.plotCh2D.chshowstr] ,'Color','red');
                end
            end
            if obj.plotCh2D.chseltop, uistack(h_selection, 'top'); end
            axis equal;
            if isfield(obj.plotCh2D,'grid') && obj.plotCh2D.grid==1
                grid on;
            end
            set(gca, 'YTick',-80:10:80); %yticks(-80:10:80);
            set(gca, 'XTick',-100:10:70); %xticks(-100:10:70);
            xlabel('MNI Y'); %predozadni souradnice
            ylabel('MNI Z'); %hornodolni
            if isfield(obj.plotCh2D,'background') && obj.plotCh2D.background==0
                set(gca,'color','none'); %zadne bile pozadi, pak ani v corelu
            end
            
            %rozhybani obrazku            
            set(obj.plotCh2D.fh,'KeyPressFcn',@obj.hybejPlot2D); 
            set(obj.plotCh2D.fh, 'WindowButtonDownFcn', @obj.hybejPlot2Dclick);
        end
%         function obj = SaveAUCPlotHandle(obj,fh) 
%             obj.plotCh2D.plotAUCH = fh; %ulozim handle na CStat.AUCPlot funkci,abych ji mohl volat z grafu ChannelPlot2D
%         end
        function tag= PacientTag(obj)
            %vraci tag pacienta, napriklad p73
            if isfield(obj.H,'patientTag'), tag = obj.H.patientTag; else, tag=obj.H.subjName; end
        end
        function [MNI_coors]= GetMNI(obj,channels)   
            %vraci koordinaty MNI pro Jirkovy skripty na SEEG-vizualizaci
            if ~exist('channels','var'), channels = 1:size(obj.H.channels,2); end
            MNI_coors = repmat(struct('MNI_x',0,'MNI_y',0,'MNI_z',0),1,numel(channels));
            
            for ch = 1:numel(channels)
                if strcmpi(obj.H.channels(ch).signalType,'SEEG')
                     MNI_coors(ch).MNI_x = obj.H.channels(ch).MNI_x;
                     MNI_coors(ch).MNI_y = obj.H.channels(ch).MNI_y;
                     MNI_coors(ch).MNI_z = obj.H.channels(ch).MNI_z;
                end
            end
        end
        function [names,neurologyLabels]=GetChNames(obj,channels)
            %vrati jmena vsech kanalu jako cell array
            if ~exist('channels','var'), channels = 1:size(obj.H.channels,2); end
            names = cell(numel(channels),1);
            neurologyLabels = cell(numel(channels),1);
            for ch = 1: numel(channels)
                if strcmpi(obj.H.channels(ch).signalType,'SEEG')
                    names{ch}=obj.H.channels(ch).name;
                    neurologyLabels{ch}=obj.H.channels(ch).neurologyLabel;
                end
            end
        end
        function [epiInfo]=GetChEpiInfo(obj,channels)
            %vrati info epilepsii v kanalech - 1vs0vsNaN - seizureOnset a interictalOften dohromady
            if ~exist('channels','var'), channels = 1:size(obj.H.channels,2); end
            if isfield(obj.H.channels,'seizureOnset')
                epiInfo = double([obj.H.channels(channels).seizureOnset]' |  [obj.H.channels(channels).interictalOften]'); %vrati 1 pokud je jedno nebo druhe 1
            else
                epiInfo = nan(numel(channels),1); %neznam epiinfo
                warning(['no epiinfo in the header: ' obj.H.subjName]);
            end
        end
            
        function tch = GetTriggerCh(obj)
            %ziska cislo trigerovaciho kanalu, pokud existuje
            tch = find(strcmp({obj.H.channels.signalType}, 'triggerCh')==1);         
        end
        function [tch,obj] = SetTriggerCh(obj,tch,trigger)
            %nastavi, ze kanal je trigger
            if ~exist('trigger','var'), trigger = 1; end 
            if trigger == 1  %defaultne nastavi, ze kanal je trigger
                obj.H.channels(tch).signalType = 'triggerCh';
                obj.H.triggerCH = tch;
            else
                obj.H.channels(tch).signalType = 'SEEG'; %nebo nastavi ze to naopak neni trigger
                obj.H.triggerCH = [];
            end
            tch = find(strcmp({obj.H.channels.signalType}, 'triggerCh')==1);
        end
        function [obj]= SetFilterMatrix(obj,filterMatrix)
            %ulozi filterMatrix pri zmene reference pro pozdejsi pouziti pri RjEpochCh
            obj.filterMatrix = filterMatrix;
        end
        function [mozek] = GetBrainNames(obj) %#ok<MANU>
            % najde popisy mist v mozku podle tabulky od Martina    
            % TODO, co kdyz jmeno obsahuje pomlcku nebo zavorku?
            BA = load('BrainAtlas_zkratky.mat'); %nactu tabulku od Martina Tomaska, 
            obj.BrainAtlas_zkratky = BA.BrainAtlas_zkratky; %ulozim do obj.BrainAtlas_zkratky
            % sloupce plnynazev, structure,region, zkratka, alt zkratka 2-4            
            mozek = cell(numel(obj.sortorder),6);            
            for ich=1:numel(obj.sortorder)  %kanaly podle aktualniho filtru              
                ch = obj.sortorder(ich);
                mozek{ich,1} = ch;
                mozek{ich,2} = obj.H.channels(ch).name; %jmeno kanalu
                mozek{ich,3} = obj.H.channels(ch).neurologyLabel; %neurologyLabel
                [C,matches] = strsplit(obj.H.channels(ch).neurologyLabel,{'(',')','-','/'},'CollapseDelimiters',false);                
                structure = cell(numel(C),2);
                for j = 1:numel(C) %pro kazdou strukturu v neurologylabel
                    [C{j},structure(j,:)] = obj.brainlabel(C{j});  %predelam zkratku na plny nazev                  
                end
                mozek{ich,4} = strjoin(C,matches);  %puvodni label zase slozim, pri pouziti plnych jmen struktur
                %structure(cellfun(@isempty,C),:) = []; %vymazu ty casti structure, pro ktere nebyl zadny match
                mozek{ich,5} = CHHeader.joinNoDuplicates(structure(:,1),matches); 
                mozek{ich,6} = CHHeader.joinNoDuplicates(structure(:,2),matches);                 
            end
        end
        function [seizureOnset,interIctal]=GetSeizures(obj)
            %vrati indexy kanalu, ktere maji seizureOnset=1 nebo interictalOften=1
            if isfield(obj.H.channels,'seizureOnset')
               seizureOnset = find(cellfun(@any,{obj.H.channels.seizureOnset}));
               interIctal = find(cellfun(@any,{obj.H.channels.interictalOften}));
            else
               seizureOnset = [];
               interIctal = [];
            end    
        end
        function obj = ChangeReference(obj,ref)
            assert(any(ref=='heb'),'neznama reference, mozne hodnoty: h e b');
            H = obj.H;  %kopie headeru
            switch ref %jaky typ reference chci
                case 'h'  %headbox          
                    filterSettings.name = 'car'; % options: 'car','bip','nan'
                    filterSettings.chGroups = 'perHeadbox';        % options: 'perHeadbox' (=global CAR), OR 'perElectrode' (=local CAR per el. shank)
                case 'e'  %elektroda
                    filterSettings.name = 'car'; % options: 'car','bip','nan'
                    filterSettings.chGroups = 'perElectrode';        % options: 'perHeadbox' (=global CAR), OR 'perElectrode' (=local CAR per el. shank)
                case 'b'  %bipolarni
                    filterSettings.name = 'bip'; % options: 'car','bip','nan'           
            end
            filterMatrix = createSpatialFilter_kisarg(H, numel(H.selCh_H), filterSettings,obj.RjCh); %ve filterMatrix uz nejsou rejectovane kanaly                        
                %assert(size(rawData,2) == size(filterMatrix,1));            
            if ref=='b' %u bipolarni reference se mi meni pocet kanalu
                H.channels = struct; 
                selCh_H = zeros(1,size(filterMatrix,2)); 
                for ch = 1:size(filterMatrix,2) % v tehle matrix jsou radky stare kanaly a sloupce nove kanaly - takze cyklus pres nove kanaly
                    oldch = find(filterMatrix(:,ch)==1); %cislo stareho kanalu ve sloupci s novym kanalem ch
                    fnames = fieldnames(obj.H.channels(oldch)); %jmena poli struktury channels
                    for f = 1:numel(fnames) %postupne zkopiruju vsechny pole struktury, najednou nevim jak to udelat
                        fn = fnames{f};
                        H.channels(ch).(fn) = obj.H.channels(oldch).(fn); 
                    end
                    H.channels(ch).name = ['(' H.channels(ch).name '-' obj.H.channels(filterMatrix(:,ch)==-1).name ')']; %pojmenuju kanal jako rozdil
                    H.channels(ch).neurologyLabel = ['(' H.channels(ch).neurologyLabel '-' obj.H.channels(filterMatrix(:,ch)==-1).neurologyLabel ')'];  %oznaceni od Martina Tomaska
                    H.channels(ch).ass_brainAtlas = ['(' H.channels(ch).ass_brainAtlas '-' obj.H.channels(filterMatrix(:,ch)==-1).ass_brainAtlas ')']; 
                    MNI = [ H.channels(ch).MNI_x H.channels(ch).MNI_y H.channels(ch).MNI_z; ... 
                           obj.H.channels(filterMatrix(:,ch)==-1).MNI_x obj.H.channels(filterMatrix(:,ch)==-1).MNI_y obj.H.channels(filterMatrix(:,ch)==-1).MNI_z ...
                           ]; %matice MNI souradnic jednoho a druheho kanalu
                    MNI2=(MNI(1,:) + MNI(2,:))/2; % prumer MNI souradnic - nova souradnice bipolarniho kanalu
                    H.channels(ch).MNI_x =  MNI2(1); 
                    H.channels(ch).MNI_y =  MNI2(2); 
                    H.channels(ch).MNI_z =  MNI2(3); 
                    H.channels(ch).MNI_dist = sqrt( sum((MNI(1,:) - MNI(2,:)).^2));  %vzdalenost mezi puvodnimi MNI body                                                           
                    selCh_H(ch) = ch;
                end 
                H.selCh_H = selCh_H; 
%                 obj.ChangeReferenceRjEpochCh(filterMatrix); %prepocitam na bipolarni referenci i RjEpochCh                                
            end
            obj.H = H; %prepisu puvodni ulozeny header
            switch ref
                case 'h', obj.reference = 'perHeadbox'; 
                case 'e', obj.reference = 'perElectrode';  
                case 'b', obj.reference = 'Bipolar';  obj.ChannelGroups([],[],1); %forcechgroups , recompute channel groups and els
            end
            obj.SetFilterMatrix(filterMatrix); %uchovam si filterMatrix na pozdeji, kvuli prepoctu RjEpochCh
        end
        function obj = SortChannels(obj,by)
            %seradi kanaly podle vybrane MNI souradnice a ulozi do sortorder                        
            if isfield(obj.plotCh2D,'chshow') 
                chshow = obj.plotCh2D.chshow; %indexy aktualne vybranych kanalu
            else
                chshow = 1:numel(obj.H.channels); %seznam vsech kanalu
            end
            if ~exist('by','var')                
                obj.sortorder = chshow;                
                obj.sortedby = '';
            elseif strcmp(by,'x')
                [~,sortorder]=sort([obj.H.channels(chshow).MNI_x]);                
                obj.sortorder = chshow(sortorder);
                obj.sortedby = 'x';
            elseif strcmp(by,'y')
                [~,sortorder]=sort([obj.H.channels(chshow).MNI_y]);                
                obj.sortorder = chshow(sortorder);
                obj.sortedby = 'y';
            elseif strcmp(by,'z')
                [~,sortorder]=sort([obj.H.channels(chshow).MNI_z]); 
                obj.sortorder = chshow(sortorder);
                obj.sortedby = 'z';
            else
                disp(['nezname razeni podle ' by]); 
            end
        end
        function obj = NextSortChOrder(obj) 
            %zmeni na dalsi razeni kanalu v poradi, podle MNI souradnic
            switch obj.sortedby
               case ''
                   obj.SortChannels('x');
               case 'x'
                   obj.SortChannels('y');
               case 'y'
                   obj.SortChannels('z');
               case 'z'
                   obj.SortChannels();
            end 
        end
        function obj = FilterChannels(obj,chlabels,notchnlabels,selCh, chnum,label)
            %vyberu podle neurologylabel jen nektere kanaly k zobrazeni
            %selCh - jedno pismeno fghjkl podle oznaceni kanalu. 
            %chnum - primo zadam cisla kanalu k filtrovani
            %label - muzu nazvat vyber jak potrebuju
            % Pozor - funguje na zaklade obj.plotCh2D.selCh, ktere se vytvari pri volani ChannelPlot2D, takze to se musi spusti nejdriv a i po zmene vyberu kanalu
            % zatim se nedaji pouzit obe metody filtrovani dohromady
            
            filtered = false;
            if exist('chlabels','var') && ~isempty(chlabels)
                ChLabels = {obj.H.channels(:).neurologyLabel}';
                iL = contains(lower(ChLabels),lower(chlabels)); %prevedu oboji na mala pismena
                if exist('notchnlabels','var') && numel(notchnlabels) > 0
                    iLx = contains(lower(ChLabels),lower(notchnlabels));
                    iL = iL & ~iLx;
                    obj.plotCh2D.chshowstr = [ cell2str(chlabels) ' not:' cell2str(notchnlabels)];
                else
                    obj.plotCh2D.chshowstr = cell2str(chlabels);
                end
                obj.plotCh2D.chshow = find(iL)'; %vyber kanalu k zobrazeni  , chci je mit v radku     
                obj.sortorder = obj.plotCh2D.chshow; %defaultni sort order pro tento vyber - nejsou tam cisla od 1 to n, ale cisla kanalu
                disp(['zobrazeno ' num2str(numel(obj.plotCh2D.chshow)) ' kanalu']);                
                filtered = true;
            end
            if exist('selCh','var') && ~isempty(selCh)
                klavesy = 'fghjkl';
                chshow = 1:numel(obj.H.channels);
                assert (numel(selCh)<=2,'maximum of 2 letter could be in selCh ');
                if find(ismember(klavesy,selCh)) %vrati index klavesy nektereho selCh v klavesy
                    if ~isfield(obj.plotCh2D,'selCh')
                        warning('No selCh in CH object, first run the ChannelPlot2D');
                    else
                        chshow = intersect(chshow,find(obj.plotCh2D.selCh(:,ismember(klavesy,selCh)))'); %indexy kanalu se znackou f-l                    
                        filtered = true;  
                    end
                end
                if ismember('r',selCh) %rejected channels, nekde v selCh je r
                    chshow = intersect(chshow,obj.RjCh);
                    %obj.plotCh2D.chshowstr = 'rj'; 
                    filtered = true;
                elseif ismember('n',selCh) %NOT rejected channels,, nekde v selCh je r
                    chshow = intersect(chshow,setdiff(obj.H.selCh_H,obj.RjCh));
                    %obj.plotCh2D.chshowstr = 'nrj'; 
                    filtered = true;                
                end
                if filtered
                    obj.plotCh2D.chshow = chshow;
                    obj.sortorder = obj.plotCh2D.chshow;
                    if exist('label','var') && ~isempty(label)
                        obj.plotCh2D.chshowstr = label;
                    else
                        obj.plotCh2D.chshowstr = selCh; 
                    end
                    disp(['zobrazeno ' num2str(numel(obj.plotCh2D.chshow)) ' kanalu: ' obj.plotCh2D.chshowstr]); 
                end                                
            end
            if exist('chnum','var') && ~isempty(chnum)
                if size(chnum,1) > size(chnum,2), chnum = chnum'; end %chci mit cisla kanalu v radku
                obj.plotCh2D.chshow = chnum; %priradim primo cisla kanalu
                obj.sortorder = obj.plotCh2D.chshow; %defaultni sort order pro tento vyber - nejsou tam cisla od 1 to n, ale cisla kanalu
                if exist('label','var') && ~isempty(label)
                    obj.plotCh2D.chshowstr = label;
                else
                    obj.plotCh2D.chshowstr = 'chnum';
                end
                disp(['zobrazeno ' num2str(numel(obj.plotCh2D.chshow)) ' kanalu: ' obj.plotCh2D.chshowstr]);                
                filtered = true;
            end
            if ~filtered
                obj.plotCh2D.chshow = 1:numel(obj.H.channels);
                obj.plotCh2D.chshowstr = '';
                obj.sortorder = 1:numel(obj.H.channels); %defaultni sort order pro vsechny kanaly
                disp('zobrazeny vsechny kanalu');
            end
            notify(obj, 'FilterChanged');
        end
        function obj = Plot3DBoundary(obj)
            %vykresli obrys mozku ve vsech rozmerech do 3d grafu
            %pokud boundary neni vypocitana, spocita ji a ulozi do obj.plotCh3D.BrainBoundaryXYZ
            %netvori graf, kresli do existujiciho a aktivniho ChannelPlot, a predpoklada hold on;
            dimenze = [2 3; 1 3; 2 1]; %xy, xz a yz 
            load('GMSurfaceMesh.mat');            
            for d = 1:3
                stred = min(GMSurfaceMesh.node(:,d)) + range(GMSurfaceMesh.node(:,d))/2;                
                if ~isfield(obj.plotCh3D,'BrainBoundaryXYZ') || numel(obj.plotCh3D.BrainBoundaryXYZ) < 3
                    obj.plotCh3D.BrainBoundaryXYZ = cell(3,1);
                end
                if isempty(obj.plotCh3D.BrainBoundaryXYZ{d})
                    obj.plotCh3D.BrainBoundaryXYZ{d} = boundary(GMSurfaceMesh.node(:,dimenze(d,1)),GMSurfaceMesh.node(:,dimenze(d,2))); %spocitam si 3d hranici mozku                
                end
                XYZ = zeros(numel(obj.plotCh3D.BrainBoundaryXYZ{d}),3);
                XYZ(:,d) = repmat(stred,size(XYZ,1),1);
                XYZ(:,dimenze(d,1)) = GMSurfaceMesh.node(obj.plotCh3D.BrainBoundaryXYZ{d},dimenze(d,1));
                XYZ(:,dimenze(d,2)) = GMSurfaceMesh.node(obj.plotCh3D.BrainBoundaryXYZ{d},dimenze(d,2));                                    
                plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3));               
            end
            
        end
    end
    
    %  --------- privatni metody ----------------------
    methods (Access = private)
          function obj = SelChannels(obj)
            % selected channels: signal type = iEEG
            % kopie z exampleSpatialFiltering.m
            selCh_H = [];
            selSignals = {'SEEG', 'ECoG-Grid', 'ECoG-Strip'};           % select desired channel group
            for ch = 1:size(obj.H.channels,2)
                if isfield(obj.H.channels(ch), 'signalType')
                    if ismember(obj.H.channels(ch).signalType, selSignals)
                        selCh_H = [selCh_H, ch]; %#ok<AGROW> %obj.H.channels(ch).numberOnAmplifier - kamil 14.6.2016 numberOnAmplifier nefunguje v pripade bipolarni reference kde jsou jina cisla kanalu
                    end
                end
            end
            
            obj.H.selCh_H = selCh_H;
          end  
          function groups = getChannelGroups(obj,subjname,chnsel)
             %vraci cell array, kazda bunka jsou cisla kanalu ve skupine
             % subjname - jestli jsou skupiny podle subjname (napr p173) nebo elektrod (napr A)             
             
             if ~exist('subjname','var') || isempty(subjname), subjname = iff(  strcmp(obj.classname,'CHilbertMulti')'); end %defaultne podle elektod
             if ~exist('chnsel','var') || isempty(chnsel), chnsel = 1:numel(obj.H.channels); end %defaultne podle elektod
             strprev = '';
             groups = cell(1,1);
             groupN = 1;
             chgroup = [];             
             if subjname
                expr = '^\w+'; %pismena i cisla,napriklad p173, konci mezerou nebo zavorkou za tim
             elseif strcmp(obj.reference,'Bipolar') && obj.H.channels(1).name(1) == '(' % ustarych data jsou bipolarni nazvy bez zavorky
                expr = '^\(([a-zA-Z]+)';
             else
                expr = '^[a-zA-Z]+';
             end
             for ch = chnsel
                 if strcmp(obj.H.channels(ch).signalType,'SEEG')                     
                     str = regexp(obj.H.channels(ch).name,expr,'match');   %jeden nebo vice pismen na zacatku                  
                     if ~strcmp(str{1},strprev) %jiny nez predchozi pacient/elektroda
                         if ch ~= chnsel(1) %pokud to neni prvni kanal
                            groups{groupN} = chgroup; %uzavru skupinu
                            groupN = groupN + 1;                            
                         end
                         chgroup = ch;
                         strprev = str{1};    
                     else %stejny pacient/elektroda jako u minuleho kanalu
                         chgroup = [chgroup ch]; %#ok<AGROW>
                     end
                 end
             end
             groups{groupN} = chgroup;
          end
          function [bl,structure] = brainlabel(obj,label)
            %vrati jmeno struktury podle tabulky, pripadne doplni otaznik na konec             
            if isempty(label)
                bl = ''; structure={'',''}; return; 
            end
            znak = strfind(label,'?');
            if ~isempty(znak)
                label = label(1:znak-1); % label bez otazniku
            end
            iL = ''; %index=radek v BrainAtlas_zkratky
            col = 3; %slopec ve zkratkach v BrainAtlas_zkratky
            %nejdriv zkusim najit presny vyskyt
            while isempty(iL) && col < size(obj.BrainAtlas_zkratky,2)
                col = col + 1; %zkratky zacinaji v ctvrtem sloupci
                if isempty(obj.BrainAtlas_zkratky(:,col)) %pokud je tahle alternativa zkratky prazdna, nebudu hledat dal
                    break; 
                end
                iL = find(strcmpi(obj.BrainAtlas_zkratky(:,col),label)); %nezavisle na velikosti pismen
            end
            if isempty(iL) %zadne label nenalezeno
                bl = '';  structure={'',''};
            elseif(numel(iL)>1) % vic nez jedno label nalezeno
                bl = [num2str(numel(iL)) ' labels'];
                structure={'',''}; 
            else
                bl = obj.BrainAtlas_zkratky{iL,1};            
                if ~isempty(znak)
                    bl = [bl '?'];
                end
                structure={obj.BrainAtlas_zkratky{iL,2},obj.BrainAtlas_zkratky{iL,3}}; 
                structure(cellfun(@isempty,structure))={''}; %potrebuju mit prvky chararray
            end
          end         
          function obj = hybejPlot3D(obj,~,eventDat)
              switch eventDat.Key
                  case 's'    %sagital view                                       
                      obj.plotCh3D.view =  [1 0 0]; %zprava  
                      obj.plotCh3D.pohled = 's';
                      view(obj.plotCh3D.view); 
                  case {'c','f'} %coronal = predozadni, frontal                      
                      obj.plotCh3D.view =  [0 -1 0];%zepredu
                      obj.plotCh3D.pohled = 'c';
                      view(obj.plotCh3D.view); %zleva
                  case {'h','a'} %horizontal = hornodolni nebo axial                        
                      obj.plotCh3D.view =  [0 0 1]; %shora
                      obj.plotCh3D.pohled = 'h';
                      view(obj.plotCh3D.view); 
                  case 'space'
                     if isfield(obj.plotCh3D,'boundary') %prepinam v grafu cely scatter s jen hranici mozku - hlavne kvuli kopirovani do corelu
                         obj.plotCh3D.boundary  = 1 - obj.plotCh3D.boundary;
                     else
                         obj.plotCh3D.boundary  = 1;
                     end
                     obj.ChannelPlot();
                  case 'n' %names
                     if isfield(obj.plotCh3D,'labels') %prepinac neurology labels v grafu - nic, neurologyLabels, channel names
                         obj.plotCh3D.labels  = obj.plotCh3D.labels + 1;
                         if obj.plotCh3D.labels > 3, obj.plotCh3D.labels =0; end
                     else
                         obj.plotCh3D.labels  = 0;
                     end
                     obj.ChannelPlot();
                  case 'r'
                     %dialog na vlozeni souradnic roi hodnoty
                    answ = inputdlg('Enter x,y,z a edge size:','define ROIs', [2 50],{num2str(obj.plotCh3D.roi)});
                    if numel(answ)>0  %odpoved je vzdy cell 1x1 - pri cancel je to cell 0x0
                        if isempty(answ{1}) %pokud vlozim hvezdicku nebo nic, chci znovy spocitat max a min
                           obj.plotCh3D.roi = [];
                        else %jinak predpokladam 4 hodnoty
                           data = str2num(answ{:});  %#ok<ST2NM>
                           if size(data,2) == 4 && size(data,1) > 0 %pokud nejsou 4 hodnoty ve sloupcich, nedelam nic. Alespon 1 radek
                             obj.plotCh3D.roi = data;
                           end
                        end
                    end
                    obj.ChannelPlot();
                  case 'p'
                    %zobrazeni pozic vsech kanalu jako tecek
                    if isfield(obj.plotCh3D,'allpoints') 
                       obj.plotCh3D.allpoints  = 1 - obj.plotCh3D.allpoints;
                    else
                       obj.plotCh3D.allpoints  = 1;
                    end
                    obj.ChannelPlot();
                  case 'z'
                    if isfield(obj.plotCh3D,'zoom') %prepinac zoom/no zoom
                       obj.plotCh3D.zoom  = obj.plotCh3D.zoom + 1;
                       if obj.plotCh3D.zoom > 2, obj.plotCh3D.zoom = 0; end
                    else
                       obj.plotCh3D.zoom  = 1; %vykreslim zoom jen kanalu
                    end 
                    obj.ChannelPlot();
                 case 'w' %zapinani a vypinani prazdneho pozadi obrazku
                    if isfield(obj.plotCh3D,'background') 
                        obj.plotCh3D.background = 1-obj.plotCh3D.background;
                    else
                        obj.plotCh3D.background = 0;
                    end
                    obj.ChannelPlot();
                  case 'o' %reorder channels, so that highest vals will be in front
                    obj.plotCh3D.reorder = 1-obj.plotCh3D.reorder;  
                    obj.ChannelPlot();
                  case 'l' %lines - spojnice kanalu - zadne, pacienti
                    %zobrazeni pozic vsech kanalu jako tecek
                    if isfield(obj.plotCh3D,'lines') 
                       obj.plotCh3D.lines  = 1 - obj.plotCh3D.lines;
                    else
                       obj.plotCh3D.lines  = 1;
                    end
                    obj.ChannelPlot();
              end
          end
          function obj = hybejPlot2D(obj,~,eventDat) 
              iCh = find(obj.plotCh2D.ch_displayed==obj.sortorder(obj.plotCh2D.chsel)); %index v obj.plotCh2D.ch_displayed
              switch eventDat.Key
                  case {'rightarrow','c'} %dalsi kanal
                      if numel(obj.plotCh2D.ch_displayed) >= iCh + 1
                        ch = min( [obj.plotCh2D.ch_displayed(iCh + 1), obj.plotCh2D.ch_displayed(end)]);
                        obj.ChannelPlot2D( find(obj.sortorder==ch)); %#ok<FNDSB>
                      end
                  case 'pagedown' %skok o 10 kanalu dopred
                      if numel(obj.plotCh2D.ch_displayed) >= iCh + 10
                        ch = min( [obj.plotCh2D.ch_displayed(iCh + 10) , obj.plotCh2D.ch_displayed(end)]);
                        obj.ChannelPlot2D( find(obj.sortorder==ch)); %#ok<FNDSB>
                      end
                  case {'leftarrow','z'} %predchozi kanal
                      if iCh - 1 >= 1
                        ch = max( [obj.plotCh2D.ch_displayed(iCh - 1) , obj.plotCh2D.ch_displayed(1)]);
                        obj.ChannelPlot2D( find(obj.sortorder==ch)); %#ok<FNDSB>
                      end
                  case 'pageup' %skok 10 kanalu dozadu
                      if iCh - 10 >= 1
                        ch = max( [obj.plotCh2D.ch_displayed(iCh - 10), obj.plotCh2D.ch_displayed(1)]);
                        obj.ChannelPlot2D( find(obj.sortorder==ch)); %#ok<FNDSB>
                      end
                  case 'return' %zobrazi obrazek mozku s vybranych kanalem                   
                      obj.plotCh2D.plotChH(obj.plotCh2D.chsel); %vykreslim @obj.PlotResponseCh                     
                      figure(obj.plotCh2D.fh); %dam puvodni obrazek dopredu
                  case 'home' %skok na prvni kanal
                      obj.ChannelPlot2D( find(obj.sortorder==obj.plotCh2D.ch_displayed(1))); %#ok<FNDSB>
                  case 'end' %skok na posledni kanal
                      obj.ChannelPlot2D( find(obj.sortorder==obj.plotCh2D.ch_displayed(end))); %#ok<FNDSB>
                  case 'period'     % prepinani razeni kanalu
                      sortorder0 = obj.sortorder; %musi si ulozit stare razeni, abych potom nasel ten spravny kanal
                      obj.NextSortChOrder();                   
                      obj.ChannelPlot2D(find(obj.sortorder==sortorder0(obj.plotCh2D.chsel))); %#ok<FNDSB> %takhle zustanu na tom stejnem kanale 
                  case {'numpad6','d'}     % skok na dalsi oznaceny kanal   
                    if isfield(obj.plotCh2D,'selCh') 
                       selCh = find(any(obj.plotCh2D.selCh,2)); %seznam cisel vybranych kanalu
                       iselCh = find(ismember(obj.sortorder,selCh)); %indexy vybranych kanalu v sortorder
                       chn2 = iselCh(find(iselCh>obj.plotCh2D.chsel,1)); %dalsi vyznaceny kanal
                       obj.ChannelPlot2D( iff(isempty(chn2),obj.plotCh2D.chsel,chn2) ); %prekreslim grafy                        
                    end                   
                  case {'numpad4','a'}     % skok na predchozi oznaceny kanal
                    if isfield(obj.plotCh2D,'selCh')  
                       selCh = find(any(obj.plotCh2D.selCh,2)); %seznam cisel vybranych kanalu
                       iselCh = find(ismember(obj.sortorder,selCh)); %indexy vybranych kanalu v sortorder
                       chn2 =  iselCh(find(iselCh < obj.plotCh2D.chsel,1,'last')) ;
                       obj.ChannelPlot2D( iff(isempty(chn2),obj.plotCh2D.chsel,chn2) ); %prekreslim grafy
                    end
                  case 'space'
                     if isfield(obj.plotCh2D,'boundary') %prepinam v grafu cely scatter s jen hranici mozku - hlavne kvuli kopirovani do corelu
                         obj.plotCh2D.boundary  = 1 - obj.plotCh2D.boundary;
                     else
                         obj.plotCh2D.boundary  = 1;
                     end
                     obj.ChannelPlot2D();
                  case {'f','g','h'}
                      obj.plotCh2D.marks( eventDat.Key - 'f' + 1) = 1 - obj.plotCh2D.marks( eventDat.Key - 'f' + 1);
                      obj.ChannelPlot2D();
                  case {'j','k','l'}
                      obj.plotCh2D.marks( eventDat.Key - 'f' ) = 1 - obj.plotCh2D.marks( eventDat.Key - 'f' );
                      obj.ChannelPlot2D();
                  case {'tab'} %zapinani a vypinani gridu
                      if isfield(obj.plotCh2D,'grid') 
                          obj.plotCh2D.grid = 1-obj.plotCh2D.grid;
                      else
                         obj.plotCh2D.grid = 1;
                      end
                      obj.ChannelPlot2D();
                  case 'w' %zapinani a vypinani prazdneho pozadi obrazku
                      if isfield(obj.plotCh2D,'background') 
                          obj.plotCh2D.background = 1-obj.plotCh2D.background;
                      else
                         obj.plotCh2D.background = 0;
                      end
                      obj.ChannelPlot2D();
                  case 'p' %vybrany kanal je zluty na popredi /pozadi
                      obj.plotCh2D.chseltop = 1-obj.plotCh2D.chseltop;
                      obj.ChannelPlot2D();
                  case 'n' %moznost vypnout / zapnout zobrazeni jmen kanalu
                      obj.plotCh2D.names = obj.plotCh2D.names + 1; 
                      if obj.plotCh2D.names == 3, obj.plotCh2D.names =0; end % meni se postupne hodoty 0 1 2
                      obj.ChannelPlot2D();    
                  case 's'
                      obj.plotCh2D.lines = obj.plotCh2D.lines + 1;
                      if obj.plotCh2D.lines == 2, obj.plotCh2D.lines = -1; end %hodnoty -1 0 1, -1=nezobrazovat neoznacene kanaly, 0=nezobrazovat cary, 1=zobrazovat
                      obj.ChannelPlot2D();
                  case 't' %barvy oznaceni kanalu fghhjkl jsou pruhledne nebo ne
                      obj.plotCh2D.transparent = 1-obj.plotCh2D.transparent;
                      obj.ChannelPlot2D();
%                   case 'r' %zobrazi obrazek mozku s vybranych kanalem                   
%                       obj.plotCh2D.plotAUCH(obj.plotCh2D.chsel); %vykreslim @obj.PlotResponseCh    %tady to hlasi error Undefined function or variable 'obj.CS.AUCPlot'. Jak to?                  
%                       figure(obj.plotCh2D.fh); %dam puvodni obrazek dopredu     
              end              
          end
          function hybejPlot2Dclick(obj,h,~)
              mousept = get(gca,'currentPoint');
              x = mousept(1,1); y = mousept(1,2); %souradnice v grafu              
              xy = get(h, 'currentpoint'); %souradnice v pixelech
              pos = get(gcf, 'Position'); 
              width = pos(3);
              subp = xy(1) > width/2; 
              chns_mni = [];
              if isfield(obj.plotCh2D,'ch_displayed')
                  chshow = obj.plotCh2D.ch_displayed; %seznam skutecne zobrazenych kanaly
              elseif isfield(obj.plotCh2D,'chshow') && numel(obj.plotCh2D.chshow) < numel(obj.H.channels) 
                  chshow = obj.plotCh2D.chshow;
              else
                  chshow = 1:numel(obj.H.channels);
              end
              if subp == 0  %axialni graf, subplot 1
                  if x > -70 && x < 70 && y > -120 && y < 90
                    chns_mni = [[obj.H.channels(chshow).MNI_x]' , [obj.H.channels(chshow).MNI_y]'];                  
                  end
              else         %sagitalni graf, subplot 2
                  if x > -100 && x < 70 && y > -100 && y < 90
                    chns_mni = [[obj.H.channels(chshow).MNI_y]' , [obj.H.channels(chshow).MNI_z]'];   
                  end
              end
              if ~isempty(chns_mni)
                  [ch,~] = dsearchn(chns_mni,[x y]); %najde nejblizsi kanal a vzdalenost k nemu                   
                  obj.ChannelPlot2D(find(obj.sortorder==chshow(ch))); %#ok<FNDSB>   
                  if isfield(obj.plotCh2D,'plotChH')
                    obj.plotCh2D.plotChH(obj.plotCh2D.chsel); %vykreslim @obj.PlotResponseCh  
                    figure(obj.plotCh2D.fh); %dam puvodni obrazek dopredu
                  end
                  
              end
          end
    end
    methods (Access = private,Static)
         function str = joinNoDuplicates(C,delim)
              switch numel(C)
                  case 1 %nasel jsem i v bipolarni referenci neurologyLabel jen PHG, bez pomlcky a zavorek, nevim jak
                      str = C{1};
                  case 4 %jen dve struktury napriklad PHG-PHG
                      if strcmp(C{2},C{3})
                          str = strjoin(C([1 2 4]),delim([1 3]));
                      else
                          str = strjoin(C,delim); 
                      end
                  case 5 %tri struktury, treba PHG/FG-FG aj
                      if strcmp(C{2},C{3}) && strcmp(C{3},C{4}) %vsechny stejne
                          str = strjoin(C([1 2 5]),delim([1 4]));
                      elseif strcmp(C{2},C{3}) %prvni a druhy stejny
                          str = strjoin(C([1 2 4 5]),delim([1 3 4]));
                      elseif strcmp(C{3},C{4}) %druhy a treti stejny
                          str = strjoin(C([1 2 3 5]),delim([1 2 4]));
                      else
                          str = strjoin(C,delim); 
                      end
                  otherwise %ctyri struktury, napriklad PHG/FG-PHG/FG a ruzne jine divnosti
                      if strcmp(C{2},C{3}) && strcmp(C{3},C{4}) && strcmp(C{4},C{5}) 
                          str = strjoin(C([1 2 6]),delim([1 5]));
                      else      
                          if strcmp(C{4},C{5})
                            C(5)=[]; delim(4) = [];
                          end
                          if strcmp(C{2},C{3})
                            C(3)=[]; delim(2) = [];
                          end 
                          if numel(C)>=3 && isempty(C{2}), C(2) = []; delim(2) = []; end
                          if numel(C)>=4 && isempty(C{3}), C(3) = []; delim(2) = []; end
                          if numel(C)>=5 && isempty(C{4})
                              C(4) = []; delim(end-1) = []; 
                          end
                          if numel(C)>=6 && isempty(C{5}), C(5) = []; delim(end-1) = []; end                                                    
                          if numel(C)>=6 && strcmp(C{4},C{5})
                            C(5)=[]; delim(4) = [];
                          end
                          if numel(C)>=4 && strcmp(C{2},C{3})
                            C(3)=[]; delim(2) = [];
                          end
                          if numel(delim)>numel(C)+1, delim(2)=[]; end
                          delim{1} = '(';
                          delim{end}=')';
                          if(numel(delim)==3), delim{2} = '-'; end
                          str = strjoin(C,delim);                       
                      end                    
              end
          end
    end
    
end

