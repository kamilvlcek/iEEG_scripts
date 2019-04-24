classdef CHHeader < matlab.mixin.Copyable %je mozne kopirovat pomoci E.copy();
    %CHHEADER Trida na praci s headerem od Jirky Hammera
    %Kamil Vlcek, FGU AVCR, since 2016 04
    
    properties (Access = public)
        H; %data od Jirky Hammera
        %E; % unikatni jmena elektrod, napr RG, I, GC ...
        chgroups;
        els;
        RjCh;
        BrainAtlas_zkratky; %tabulka od Martina Tomaska
        filterMatrix; %kopie filterMatrix, vytvari se pri zmene reference
        sortorder; %index serazenych kanalu
        sortedby; %podle ceho jsou kanaly serazeny
        plotCh2D; %udaje o 2D grafu kanalu ChannelPlot2D, hlavne handle
        plotCh3D; %udaje o 3D grafu kanalu ChannelPlot, hlavne handle
        
    end
    %#ok<*PROPLC>
    methods (Access = public)
        function obj = CHHeader(H,filename)
            %konstruktor
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
             
        end
        
        function [obj, chgroups, els] = ChannelGroups(obj,chnsel)
            %vraci skupiny kanalu (cisla vsech channels na elekrode) + cisla nejvyssiho kanalu v na kazde elektrode v poli els
            if ~exist('chnsel','var') || isempty(chnsel)
                if isempty(obj.chgroups)
                    chgroups = getChannelGroups_kisarg(obj.H,'perElectrode');
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
            else
                if obj.plotCh3D.allpoints %pokud chci zobrazovat i ostatni kanal jako tecky
                   chgroups = {chnsel,setdiff(obj.H.selCh_H,chnsel)}; %do druhe skupiny dam ostatni kanaly
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
        function [XYZ,obj] = ChannelPlot(obj,pohled,labels,chnvals,chnsel,selch,roi)
            %zobrazi 3D obrazek elektrod v MNI prostoru. Obrazek ma rozmery podle rozmeru mozku
            %pohled muze urcti smer pohledu s-sagital,c-coronal,h-horizontal
            %chnsel jsou cisla kanalu, pokud chci jen jejich vyber
            %selch je jedno zvyraznene cislo kanalu - index v poli chnsel
            %roi je zvyraznena krychlova oblast [ x y z edge]
            if ~exist('pohled','var') || isempty(pohled), pohled = ''; end            
            
            params = {'labels','chnvals','chnsel','selch','roi'}; %zkusim hromadne zpracovani parametru touhle nedoporucovanou metodou
            for p=1:numel(params) %parametry, ktere se ukladaji do obj.plotCh3D
                if ~exist(params{p},'var') || eval(['isempty(' params{p} ')']) %pokud neni vstupni promenna nebo je prazdna
                    if isfield(obj.plotCh3D,params{p}) %pokud ale existuje ulozena hodnota
                        eval([ params{p} ' = obj.plotCh3D.' params{p} ';']); %tak ji pouziju
                    else 
                        switch params{p}
                            case 'labels'
                                labels = 0;  %defaultne se nezobrazuji neurology labels, ale jmena kanalu
                            case 'chnvals'
                                chnvals = zeros(1, numel(obj.H.channels)); %default same nuly
                            case 'chnsel'
                                chnsel = []; %default vsechny kanaly
                            case 'selch'
                                selch = []; %default vsechny kanaly
                            case 'roi'
                                roi = []; %default zadne
                        end
                        eval(['obj.plotCh3D.' params{p} ' = ' params{p} ';']); %nastavim ulozenou hodnotu na default
                    end
                else
                    eval(['obj.plotCh3D.' params{p} ' = ' params{p} ';']); %podle vstupni promenne zmeni ulozenou hodnotu
                end
            end           
            if ~isfield(obj.plotCh3D,'allpoints'), obj.plotCh3D.allpoints = 0; end
            assert(numel(chnvals) == numel(chnsel), 'unequal size of chnvals and chnsel');
            nblocks = numel(chnvals); %pocet barev bude odpovidat poctu kanalu
            cmap = parula(nblocks+1); %+1 protoze hodnoty se budou zaokrouhlovat nahoru nebo dolu
            chnvalsN = chnvals - min(chnvals);
            chnvalsN = chnvalsN/max(chnvalsN); % normalization            
            chnvals(isnan(chnvalsN)) = 0; % in case of all zeros
            clrs = cmap(round(nblocks*chnvalsN)+1, :); % color values
            sizes = 20+200*chnvalsN;
            if isfield(obj.H.channels,'MNI_x')
                if isfield(obj.plotCh3D,'fh') && ishandle(obj.plotCh3D.fh)
                    figure(obj.plotCh3D.fh); %pouziju uz vytvoreny graf
                    clf(obj.plotCh3D.fh); %graf vycistim
                else
                    obj.plotCh3D.fh = figure('Name','ChannelPlot 3D in MNI');                     
                end          
                               
                [obj,chgroups] = obj.ChannelGroups(chnsel); %rozdeli kanaly po elektrodach do skupin. 
                 %Pokud chnsel, jsou vsecny v jedne skupine. Ale pokud obj.plotCh3D.allpoints, ve druhe skupine jsou ostatni kanaly
                
                %objekt se dobre uklada i pri poradi return values XYZ,obj
                XYZ = struct('X',0,'Y',0,'Z',0);
                for chg = 1:size(chgroups,2) 
                    group = chgroups{chg};                     
                    X = [obj.H.channels(group).MNI_x];
                    Y = [obj.H.channels(group).MNI_y];
                    Z = [obj.H.channels(group).MNI_z];     
                   
                    linestyle = iff(isempty(chnsel),'-','.'); %cara bude jina pokud je pouzite chnsel
                    plot3(X,Y,Z,linestyle,'LineWidth',2);
                    if chg==1, hold on; end                         
                    if labels
                        names = {obj.H.channels(group).neurologyLabel};
                    else
                        names = {obj.H.channels(group).name};
                    end
                    iZ = mod([1:numel(Z)], 2); iZ(iZ == 0) = -1;                    
                    if chg==1 || isempty(chnsel) || ~obj.plotCh3D.allpoints %druhou skupinu chci do kulicek jen pokud zobrazuju vsechny (chnsel je prazdne) nebo pokud nejsou v druhe skupine ostatni kanaly
                        XYZ(chg) = struct('X',X,'Y',Y,'Z',Z); %export pro scatter3 nize, ktery zobrazi ruzne velke a barevne kulicky
                        text(X,Y,Z+iZ*2,names,'FontSize', 7);
                    end
                end
                % Plot with different colors and sizes based on chnvals
                if isempty(chnvals) 
                    isizes = 1:obj.H.channels;
                elseif obj.plotCh3D.allpoints 
                    isizes = find(chnsel==[chgroups{1}]); %indexy v poli chnsel pro pouziti v poli sizes, find pracuje i hromadne
                else
                    isizes = find(chnsel==[chgroups{:}]); %indexy v poli chnsel pro pouziti v poli sizes, find pracuje i hromadne
                end               
                scatter3([XYZ.X],[XYZ.Y],[XYZ.Z],sizes(isizes),clrs(isizes,:),'filled'); %ruzne velke a barevne krouzky vsech kanalu najednou
               
                if ~isempty(selch)
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
                        view([-1 0 0]); %zleva
                    case 'c' %coronal = predozadni
                        view([0 1 0]); %zepredu
                    case 'h' %horizontal = hornodolni   
                        view([0 0 1]); %shora
                end
                axis([-75 75 -120 80 -75 85]); %zhruba velikost mozku        
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
                    caxis([min(chnvals) max(chnvals)]); 
                end %barevna skala, jen pokud jsou ruzne hodnoty kanalu
                axis equal;                               
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
            size_selCh = 7;
            x_text = -100;
            if isfield(obj.plotCh2D,'fh') && ishandle(obj.plotCh2D.fh)
                figure(obj.plotCh2D.fh); %pouziju uz vytvoreny graf
                clf(obj.plotCh2D.fh); %graf vycistim
            else
                obj.plotCh2D.fh = figure('Name',['ChannelPlot2D - ' label]);                     
            end            
                   
            if isfield(obj.plotCh2D,'chshow') && numel(obj.plotCh2D.chshow) < numel(obj.H.channels) %pokud chci zobrazovat jen cast kanalu podle chshow
                els = obj.plotCh2D.chshow;
                els0 = obj.plotCh2D.chshow; %nebudu resit zactky a konec elektrod
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
                plot(x(els0(ie):els(ie)),y(els0(ie):els(ie)),'-o'); %plot kontaktu jedne elektrody
                if obj.plotCh2D.names
                for ch = els0(ie):els(ie)
                        th = text(x(ch),y(ch),num2str(ch)); %cislo kazdeho kanalu
                        th.FontSize = 8;
                end
                end                              
            end
            if ~isempty(chsel) %pokud je vybrany nejaky kanal                
                h_selection = plot(x(chselo),y(chselo),'o','MarkerSize',size_ch,'MarkerEdgeColor','y','MarkerFaceColor','y'); 
                chstr = iff(isempty(obj.sortedby),num2str(chsel), [ num2str(obj.sortorder(chsel)) '(' obj.sortedby  num2str(chsel) ')' ]);
                title( [ 'channel ' chstr ]);
                
            end
            if ~isempty(selCh) %hromadne vybrane kanaly, zobrazne cernym koleckem
                barvy = 'rbgcmk';
                for m = size(selCh,2):-1:1 %jednu znacku za druhou
                   if  obj.plotCh2D.marks(m) %pokud se ma znacka zobrazovat
                       ch = find(selCh(:,m)); %seznam cisel vybranych kanalu pro danou znacku
                       ch = intersect(chshow,ch); 
                       plot(x(ch),y(ch),'o','MarkerSize',size_selCh,'MarkerEdgeColor',barvy(m),'MarkerFaceColor',barvy(m));
                   end
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
                plot(y(els0(ie):els(ie)),z(els0(ie):els(ie)),'-o'); %plot kontaktu jedne elektrody
                if obj.plotCh2D.names
                for ch = els0(ie):els(ie)
                    th = text(y(ch),z(ch),num2str(ch));
                    th.FontSize = 8;
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
                barvy = 'rbgcmk';
                klavesy = 'fghjkl'; %abych mohl vypsat primo nazvy klaves vedle hvezdicky podle selCh
                for m = size(selCh,2):-1:1 %jednu znacku za druhou - naposled ty prvni aby byly nahore
                    if  obj.plotCh2D.marks(m) %pokud se ma znacka zobrazovat
                       ch = find(selCh(:,m)); %seznam cisel vybranych kanalu pro danou znacku
                       ch = intersect(chshow,ch); 
                       if ~isempty(ch) %pokud jsou takove nejake vybrane kanaly
                           plot(y(ch),z(ch),'o','MarkerSize',size_selCh,'MarkerEdgeColor',barvy(m),'MarkerFaceColor',barvy(m));
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
        function [mozek] = GetBrainNames(obj)
            % najde popisy mist v mozku podle tabulky od Martina
            % TODO jeste neumi dve struktury s /
            % a s otaznikem na konci
            load('BrainAtlas_zkratky.mat'); %nactu tabulku od Martina Tomaska
            obj.BrainAtlas_zkratky = BrainAtlas_zkratky; %#ok<CPROP>
            mozek = cell(numel(obj.H.channels),2);
            for ch=1:numel(obj.H.channels)
                label = obj.H.channels(ch).neurologyLabel;   % Martinovo label tohoto kanalu
                mozek{ch,1} = label;
                znak = strfind(label,'/');
                if ~isempty(znak) %pokud se jedna o dva labely oddelene lomitkem
                    label1= label(1:znak-1);
                    label2 = label(znak+1:end);
                    mozek{ch,2} = [obj.brainlabel(label1) '/' obj.brainlabel(label2)];
                else
                    mozek{ch,2} = obj.brainlabel(label) ;
                end                
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
        function obj = FilterChannels(obj,chlabels,notchnlabels)
            %vyberu podle neurologylabel jen nektere kanaly k zobrazeni
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
                obj.plotCh2D.chshow = find(iL); %vyber kanalu k zobrazeni       
                obj.sortorder = obj.plotCh2D.chshow; %defaultni sort order pro tento vyber - nejsou tam cisla od 1 to n, ale cisla kanalu
                disp(['zobrazeno ' num2str(numel(obj.plotCh2D.chshow)) ' kanalu']);
            else
                obj.plotCh2D.chshow = 1:numel(obj.H.channels);
                obj.plotCh2D.chshowstr = '';
                obj.sortorder = 1:numel(obj.H.channels); %defaultni sort order pro vsechny kanaly
                disp('zobrazeny vsechny kanalu');
            end
            
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
          function [bl] = brainlabel(obj,label)
            %vrati jmeno struktury podle tabulky, pripadne doplni otaznik na konec             
            znak = strfind(label,'?');
            if ~isempty(znak)
                label = label(1:znak-1); % label bez otazniku
            end
            iL = find(strcmp(obj.BrainAtlas_zkratky,label));   
            if isempty(iL) %zadne label nenalezeno
                bl = '';
            elseif(numel(iL)>1) % vic nez jedno label nalezeno
                bl = [num2str(numel(iL)) ' labels'];
            else
                bl = obj.BrainAtlas_zkratky{iL,2};            
                if ~isempty(znak)
                    bl = [bl '?'];
                end
            end
          end
          function obj = hybejPlot3D(obj,~,eventDat)
              switch eventDat.Key
                  case 's'                                           
                      obj.plotCh3D.view =  [1 0 0]; %zprava
                      view(obj.plotCh3D.view); 
                  case {'c','f'} %coronal = predozadni, frontal                      
                      obj.plotCh3D.view =  [0 -1 0];%zepredu
                      view(obj.plotCh3D.view); %zleva
                  case {'h','a'} %horizontal = hornodolni nebo axial                        
                      obj.plotCh3D.view =  [0 0 1]; %shora
                      view(obj.plotCh3D.view); 
                  case 'space'
                     if isfield(obj.plotCh3D,'boundary') %prepinam v grafu cely scatter s jen hranici mozku - hlavne kvuli kopirovani do corelu
                         obj.plotCh3D.boundary  = 1 - obj.plotCh3D.boundary;
                     else
                         obj.plotCh3D.boundary  = 1;
                     end
                     obj.ChannelPlot();
                  case 'l'
                     if isfield(obj.plotCh3D,'labels') %prepinac neurology labels v grafu
                         obj.plotCh3D.labels  = 1 - obj.plotCh3D.labels;
                     else
                         obj.plotCh3D.labels  = 1;
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
                    if isfield(obj.plotCh3D,'allpoints') %prepinac neurology labels v grafu
                       obj.plotCh3D.allpoints  = 1 - obj.plotCh3D.allpoints;
                    else
                       obj.plotCh3D.allpoints  = 1;
                    end
                    obj.ChannelPlot();
              end
          end
          function obj = hybejPlot2D(obj,~,eventDat) 
              switch eventDat.Key
                  case {'rightarrow','c'} %dalsi kanal
                      obj.ChannelPlot2D( min( [obj.plotCh2D.chsel + 1 , max(obj.H.selCh_H)]));
                  case 'pagedown' %skok o 10 kanalu dopred
                      obj.ChannelPlot2D( min( [obj.plotCh2D.chsel + 10 , max(obj.H.selCh_H)]));
                  case {'leftarrow','z'} %predchozi kanal
                      obj.ChannelPlot2D( max( [obj.plotCh2D.chsel - 1 , 1]));
                  case 'pageup' %skok 10 kanalu dozadu
                      obj.ChannelPlot2D( max( [obj.plotCh2D.chsel - 10, 1]));
                  case 'return' %zobrazi obrazek mozku s vybranych kanalem                   
                      obj.plotCh2D.plotChH(obj.plotCh2D.chsel); %vykreslim @obj.PlotResponseCh                     
                      figure(obj.plotCh2D.fh); %dam puvodni obrazek dopredu
                  case 'home' %skok na prvni kanal
                      obj.ChannelPlot2D(1);
                  case 'end' %skok na posledni kanal
                      obj.ChannelPlot2D( max(obj.H.selCh_H));
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
                  case 't' %vybrany kanal je zluty na popredi /pozadi
                      obj.plotCh2D.chseltop = 1-obj.plotCh2D.chseltop;
                      obj.ChannelPlot2D();
                  case 'n' %vybrany kanal je zluty na popredi /pozadi
                      obj.plotCh2D.names = 1-obj.plotCh2D.names;
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
              if isfield(obj.plotCh2D,'chshow') && numel(obj.plotCh2D.chshow) < numel(obj.H.channels) 
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
    
end

