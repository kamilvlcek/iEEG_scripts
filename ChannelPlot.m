classdef ChannelPlot < matlab.mixin.Copyable
    properties (Access = public)
        CH@CHHeader;
        plotCh3D; %udaje o 3D grafu kanalu ChannelPlot, hlavne handle
    end
    
    methods (Access = public)
        function obj = ChannelPlot(header)
            obj.CH = header;
        end
        function obj = ChannelPlotInit(obj,plotCh3D)
            %plotCh3D enables to load fields from struct and not init other fields 
            if ~exist('plotCh3D','var')
                if ~isfield(obj.plotCh3D,'names'), obj.plotCh3D.names = 0; end %by default, no labels for channels are used
                if ~isfield(obj.plotCh3D,'labels'), obj.plotCh3D.labels = 0; end   %if to show brainlabels as channel labels
                if ~isfield(obj.plotCh3D,'labesXnames'), obj.plotCh3D.labesXnames = 0; end   %if to show brainlabels (1) or channel names (0)
                if ~isfield(obj.plotCh3D,'boundary'), obj.plotCh3D.boundary = 1; end %if to plot boundary of the brain instead the 3D mesh
                if ~isfield(obj.plotCh3D,'allpoints'), obj.plotCh3D.allpoints = 0; end %if to show position of all channel, even non significant
                if ~isfield(obj.plotCh3D,'allpointnames'), obj.plotCh3D.allpointnames = 0; end %if to show labels of all the channels
                if ~isfield(obj.plotCh3D,'zoom'), obj.plotCh3D.zoom = 0; end            
                if ~isfield(obj.plotCh3D,'reorder'), obj.plotCh3D.reorder = 0; end   %defaultne se neprerazuji kanaly podle velikosti
                if ~isfield(obj.plotCh3D,'lines'), obj.plotCh3D.lines = 0; end   %defaultne se nespojuji pacienti spojnicemi            
                if ~isfield(obj.plotCh3D,'hullindex'), obj.plotCh3D.hullindex = 0; end   %index of brainlabel to plot convex hull        
                if ~isfield(obj.plotCh3D,'fontsize'), obj.plotCh3D.fontsize = 7; end   %index of brainlabel to plot convex hull 
                if ~isfield(obj.plotCh3D,'coloruse'), obj.plotCh3D.coloruse = 0; end   %which color to use 0=according to value size, 1=according to channel marks as in 2D plot, 2=according to brain labels
                if ~isfield(obj.plotCh3D,'markertype'), obj.plotCh3D.markertype = 'o'; end   %which markertype to use in scatter3 - default is ball
                if ~isfield(obj.plotCh3D,'showclusters'), obj.plotCh3D.showclusters = 1; end   %if to show clusters of channels marked by X
            end
            if exist('plotCh3D','var') && isstruct(plotCh3D)
                fields = fieldnames(plotCh3D);
                for f = 1:numel(fields)
                    obj.plotCh3D.(fields{f}) = plotCh3D.(fields{f});
                end
            end
        end
        function [XYZ,obj] = ChannelPlot3D(obj,chnvals,chnsel,selch,roi,popis,rangeZ)
            %zobrazi 3D obrazek elektrod v MNI prostoru. Obrazek ma rozmery podle rozmeru mozku
            %pohled muze urcti smer pohledu s-sagital,c-coronal,h-horizontal
            %chnsel jsou cisla kanalu, pokud chci jen jejich vyber - musi byt stejny pocet jako chvals (hodnoty k vykresleni)
            %selch je jedno zvyraznene cislo kanalu - index v poli chnsel
            %roi je zvyraznena krychlova oblast [ x y z edge]
            %popis je text k zobrazeni na obrazku        
            
            params = {'chnvals','chnsel','selch','roi','popis','rangeZ'}; %zkusim hromadne zpracovani parametru touhle nedoporucovanou metodou
            iSEEG = contains({obj.CH.H.channels.signalType},'SEEG'); %index kanalu s EEG signalem
            for p=1:numel(params) %parametry, ktere se ukladaji do obj.plotCh3D
                if ~exist(params{p},'var') || eval(['isempty(' params{p} ')']) %pokud neni vstupni promenna nebo je prazdna
                    if isfield(obj.plotCh3D,params{p}) %pokud ale existuje ulozena hodnota
                        eval([ params{p} ' = obj.plotCh3D.' params{p} ';']); %tak ji pouziju
                    else 
                        switch params{p}                      
                            case 'chnvals'
                                chnvals = zeros(1, numel(obj.CH.H.channels(iSEEG))); %default same nuly
                            case 'chnsel'
                                chnsel = 1:numel(obj.CH.H.channels(iSEEG)) ; %default vsechny kanaly
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
                   
            obj.ChannelPlotInit();
            
            assert(numel(chnvals) == numel(chnsel), 'unequal size of chnvals and chnsel');
            [clrs,sizes,rangeZ,reverse] = obj.colors4ChannelPlot(chnsel, chnvals,rangeZ);
            %if reverse, sizes = flip(sizes); end
            if isfield(obj.CH.H.channels,'MNI_x')
                if isfield(obj.plotCh3D,'fh') && ishandle(obj.plotCh3D.fh)
                    figure(obj.plotCh3D.fh); %pouziju uz vytvoreny graf
                    [caz,cel] = view(); %current view                    
                    obj.plotCh3D.view = [caz,cel];   %store it
                    clf(obj.plotCh3D.fh); %graf vycistim      
                else
                    obj.plotCh3D.fh = figure('Name','ChannelPlot 3D in MNI');                     
                    obj.plotCh3D.isColormapReversed = 0;
                    obj.plotCh3D.view = [0 90]; %default view
                end          
                               
                [obj.CH,chgroups] = obj.CH.ChannelGroups(chnsel,obj.plotCh3D.lines); %rozdeli kanaly po elektrodach do skupin. 
                 %Pokud chnsel, jsou vsecny v jedne skupine. Ale pokud obj.plotCh3D.allpoints, ve druhe skupine jsou ostatni kanaly
                
                %objekt se dobre uklada i pri poradi return values XYZ,obj
                XYZ = struct('X',0,'Y',0,'Z',0);
                for chg = 1:size(chgroups,2) 
                    chGroup = chgroups{chg};  %channel numbers in this group                   
                    X = [obj.CH.H.channels(chGroup).MNI_x];
                    Y = [obj.CH.H.channels(chGroup).MNI_y];
                    Z = [obj.CH.H.channels(chGroup).MNI_z];                                             
                    
                    linestyle = iff(numel(chgroups)>1 && ~obj.plotCh3D.allpoints,'-','.'); %cara bude jina pokud je pouzite chnsel
                    if obj.plotCh3D.allpoints, plot3(X,Y,Z,linestyle,'LineWidth',2); end
                    if chg==1, hold on; end  
                    chnnames = {};
                    if ~obj.plotCh3D.labesXnames
                        switch obj.plotCh3D.names %which channel names to show as labels for points
                            case 1 %channel numbers
                                chnnames = num2cell(chGroup);
                            case 2 %channel names
                                chnnames = {obj.CH.H.channels(chGroup).name};
                            case 3 %neurology labels
                                chnnames = {obj.CH.H.channels(chGroup).neurologyLabel};                    
                            case 4 %pacient names
                                chnnames = {obj.CH.H.channels(chGroup).name};
                                chnnames = cellstr(extractBefore(chnnames,' ')); %vsechno pred mezerou - pro CHilbertMulti
                        end
                    else
                        switch obj.plotCh3D.labels %which brainlabels to show as labels for points
                            case 1 
                                chnnames = {obj.CH.brainlabels(chGroup).lobe};
                            case 2
                                chnnames = {obj.CH.brainlabels(chGroup).label};
                            case 3
                                chnnames = {obj.CH.brainlabels(chGroup).class};
                        end
                    end
                    iZ = mod(1:numel(Z), 2); iZ(iZ == 0) = -1;                    
                    if chg==1 || ~obj.plotCh3D.allpoints %prvni skupiny do barevnych kulicek vzdy; 
                        %druhou skupinu chci jen pokud zobrazuju vsechny (chnsel je prazdne) nebo pokud nejsou v druhe skupine ostatni kanaly
                        XYZ(chg) = struct('X',X,'Y',Y,'Z',Z); %export pro scatter3 nize, ktery zobrazi ruzne velke a barevne kulicky
                    end 
                    if ~isempty(chnnames) && (chg==1 || ~obj.plotCh3D.allpoints || (obj.plotCh3D.allpoints && obj.plotCh3D.allpointnames))                         
                        text(X+abs(iZ)*0.5,Y,Z+iZ*0.5,chnnames,'FontSize', obj.plotCh3D.fontsize);     %labels=names for channels                   
                    end
                end
                % Plot with different colors and sizes based on chnvals
                if isempty(chnvals)  %indexy vsech kanalu, nemam zadne hodnoty k vykresleni
                    isizes = 1:obj.CH.H.channels;
                elseif obj.plotCh3D.allpoints  %zobrazuju pozice vsech kanalu jako tecek (dve skupiny kanalu v chgroups - barevne kulicky + tecky)
                    isizes = find(chnsel==[chgroups{1}]); %indexy v poli chnsel pro pouziti v poli sizes, find pracuje i hromadne
                else %indexy vsech kanalu ve vsech skupinach
                    isizes = find(chnsel==[chgroups{:}]); %indexy v poli chnsel pro pouziti v poli sizes, find pracuje i hromadne
                end   
                X = [XYZ.X]; Y = [XYZ.Y]; Z = [XYZ.Z]; %souradnice pres vsechny pole struct XYZ
                if obj.plotCh3D.reorder %pokud chci seradi body podle velikosti, tak aby v prislusnem pohledu byly nejvetsi v popredi
                    reordered = 1;
                    if isequal(obj.plotCh3D.view,[0 90])%'h'                        
                        Z = sortBlikeA(sizes,Z); %nejvetsi hodnoty na nejvyssich souradnicich Z
                    elseif isequal(obj.plotCh3D.view,[180 -90])%'h'                            
                        Z = sortBlikeA(-sizes,Z); %nejvetsi hodnoty na nejvyssich souradnicich Z
                    elseif isequal(obj.plotCh3D.view,[0 0]) %'c'
                        Y = sortBlikeA(-sizes,Y); %nejvetsi hodnoty na nejnizsich souradnicich y
                    elseif isequal(obj.plotCh3D.view,[180 0]) %'c'
                        Y = sortBlikeA(sizes,Y); %nejvetsi hodnoty na nejnizsich souradnicich y    
                    elseif isequal(obj.plotCh3D.view,[90 0]) % 's'
                        X = sortBlikeA(sizes,X);
                    elseif isequal(obj.plotCh3D.view,[-90 0]) % 's'
                        X = sortBlikeA(-sizes,X);    
                    else
                        reordered = 0;
                    end
                    if reordered, annotation('textbox', [.6 0.15 .2 .1], 'String', 'REORDERED', 'EdgeColor', 'none'); end
                end
%                 
                scatter3(X,Y,Z,sizes(isizes),clrs(isizes,:),'filled',obj.plotCh3D.markertype,'MarkerEdgeColor','k'); %ruzne velke a barevne krouzky vsech kanalu najednou
               
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

                view(obj.plotCh3D.view); %restore the view stored before replotting the figure                

                if obj.plotCh3D.zoom == 0
                    axis([-75 75 -120 80 -75 85]); %zhruba velikost mozku        
                else 
                    axis([ min([XYZ.X]) max([XYZ.X]) min([XYZ.Y]) max([XYZ.Y]) min([XYZ.Z]) max([XYZ.Z]) ] );
                end
                text(-70,0,0,'LEFT');        
                text(70,0,0,'RIGHT');   
                text(0,65,0,'FRONT');        
                text(0,-115,0,'BACK');                                 
       
                if ~isfield(obj.plotCh3D,'boundary') || obj.plotCh3D.boundary
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
                if obj.plotCh3D.coloruse == 2 && ~isempty(obj.CH.brainlabels) && length(obj.CH.brainlabels)==numel(obj.CH.H.channels)
                    labels = lower({obj.CH.brainlabels.label});
                    ulabels = unique(labels); %cell array of unique brainlabels
                    barvy = distinguishable_colors(numel(ulabels));
                    for ilabel = 1:numel(ulabels)
                       if isequal( [90 0],obj.plotCh3D.view) || isequal([-90 0],obj.plotCh3D.view) %sagital
                           x = 0; y = -110; z = 80-5*ilabel;
                       elseif isequal( [0 90 ],obj.plotCh3D.view) || isequal([[ 180 -90]],obj.plotCh3D.view) %axial
                           x = -70; y = 80-5*ilabel; z = 0;
                       else
                           x = -70; y = 0; z = 80-5*ilabel;
                       end
                       text(x,y,z,ulabels{ilabel},'Color',barvy(ilabel,:),'FontSize',11,'FontWeight','bold');                       
                    end
                end
                
                obj.PlotClusters(); %plot channel clusters if any exist            
                obj.plotCh3D.dispChannels = chnsel; % ulozim vyber zobrazenych kanalu (je potreba pro klikani)
                obj.highlightChannel(); %if there is any channel to be highligted, do it
                %rozhybani obrazku            
                set(obj.plotCh3D.fh,'KeyPressFcn',@obj.hybejPlot3D);
            else
                disp('No MNI data');
            end
        end
        function PlotClusters(obj)
            %PlotClusters - plots the previosly computer clusters if any exist, according to the current popis            
            iCluster = obj.CH.GetCluster(obj.plotCh3D.popis);
            if iCluster && obj.plotCh3D.showclusters
                C = obj.CH.clusters(iCluster).C;
                plot3(C(:,1),C(:,2),C(:,3),'kx','MarkerSize',20,'LineWidth',3); %clusters on right side
                text(C(:,1)+5,C(:,2)+5,C(:,3)+5, cellstr(num2str((1:size(C,1))')),'FontSize',11,'FontWeight','bold');
                hold on
                plot3(-C(:,1),C(:,2),C(:,3),'kx','MarkerSize',20,'LineWidth',3); %clusters on left side
                text(-C(:,1)+5,C(:,2)+5,C(:,3)+5, cellstr(num2str((1:size(C,1))')),'FontSize',11,'FontWeight','bold');
            end           
            
            %older plotting of complex hull
            if ~isempty(obj.CH.hull) && obj.plotCh3D.hullindex > 0
                obj.CH.HullPlot3D(obj.plotCh3D.hullindex);
            end
                
        end

        function [clrs,sizes,rangeZ,reverse] = colors4ChannelPlot(obj,chnsel,chnvals,rangeZ)
            nblocks = numel(chnvals); %the number of colors will be the same as channels
            cmap = parula(nblocks+1); %+1 as the values will be rounded up or down
            reverse = 0; %if to reverse to color map and size of the points in scatter3D 
            if isempty(rangeZ)
                rangeZ = [min(chnvals) max(chnvals)];                 
            elseif rangeZ(1) > rangeZ(2) %pokud dam minmax v obrazenem poradi, barvy i velikosti taky v obracenem poradi
                reverse = 1;
                rangeZ = flip(rangeZ);            
                cmap = flip(cmap,1);
            end
            
            chnvalsN = chnvals - rangeZ(1); %substract minimum
            chnvalsN = chnvalsN/diff(rangeZ); % normalization  - divide by maximum => values are [0;1]          
            chnvalsN(isnan(chnvalsN)) = 0; % in case of all zeros, feplace nan to 0
            chnvalsN(chnvalsN<0) = 0; chnvalsN(chnvalsN>1) = 1; %limit the range to [0;1];    
            sizes = 20+200*iff(reverse,1-chnvalsN,chnvalsN); %velikosti kulicek 
            switch obj.plotCh3D.coloruse
                case 0 %colors based on channel vals
                    clrs = cmap(round(nblocks*chnvalsN)+1, :); % rgb color values for each channel (chns x 3), prevedu na rozsah 1-nblocks a priradim barvy
                case 1 %colors based on channel markings as in ChannelPlot2D
                    clrs = repmat([.5 .5 .5],numel(chnvals),1); %default grey color for each channel
                    if ~isempty(obj.CH.plotCh2D) % Musi byt vytvoreny plotCh2D
                        barvy = [obj.CH.plotCh2D.color_def(obj.CH.plotCh2D.color_index:end,:); obj.CH.plotCh2D.color_def(1:obj.CH.plotCh2D.color_index-1,:)]; %barvy od poradi colorindexu
                        for ci = 1:numel(obj.CH.plotCh2D.color_order) %1:size(selCh,2) %jednu znacku za druhou m = size(selCh,2):-1:1 
                           m = obj.CH.plotCh2D.color_order(ci); %order of marks by color_order, similarly to 2D channel plot
                           if  obj.CH.plotCh2D.marks(m) %if to show this mark
                               ch = find(obj.CH.plotCh2D.selCh(:,m)); %number of channels for this channel mark
                               ch = intersect(chnsel,ch); %reduced for only the selected channels
                               if numel(ch) > 0
                                   ichannel = ismember(chnsel,ch)'; %index of channels in chnsel and chnvals
                                   clrs(ichannel,:)=repmat(barvy(m,:),size(ch)); %previous color is overwriten
                               end
                           end
                        end
                    end
                case 2 %colors based on brainlabels
                    clrs = repmat([.5 .5 .5],numel(chnvals),1); %default grey color for each channel                                
                    if ~isempty(obj.CH.brainlabels) && length(obj.CH.brainlabels)==numel(obj.CH.H.channels) %if brainlabels exist, otherwise channels will be all grey
                        labels = lower({obj.CH.brainlabels.label});
                        ulabels = unique(labels); %cell array of unique brainlabels
                        barvy = distinguishable_colors(numel(ulabels));
                        for j = 1:numel(ulabels) %cycle over all brainlabels
                            if j<= size(barvy,1)
                                ch = find(contains(labels,ulabels{j}));  %channels with this brainlabel
                                ch = intersect(ch,chnsel); %reduced for only the selected channels
                                if numel(ch) > 0
                                   ichannel = ismember(chnsel,ch); %index of channels in chnsel and chnvals
                                   clrs(ichannel,:)=repmat(barvy(j,:),length(ch),1); %previous color is overwriten
                                end
                            end
                        end
                    end
                case 3 %colors according to clusters                   
                    clrs = repmat([.5 .5 .5],numel(chnvals),1); %default grey color for each channel      
                    iCluster = obj.CH.GetCluster(obj.plotCh3D.popis);
                    if iCluster
                        barvy = distinguishable_colors(size(obj.CH.clusters(iCluster).C,1));  %number of clusters
                        for j = 1:size(obj.CH.clusters(iCluster).C,1)
                            ichannels = obj.CH.clusters(iCluster).idx==j; %channels plotted that are in the current cluster
                            ch = obj.CH.clusters(iCluster).channels(ichannels);
                            ch = intersect(ch,chnsel); %reduced for only the selected channels
                            if numel(ch) > 0
                               ichannel = ismember(chnsel,ch); %index of channels in chnsel and chnvals
                               clrs(ichannel,:)=repmat(barvy(j,:),length(ch),1); %previous color is overwriten
                            end
                        end
                    end
                    
            end
            obj.plotCh3D.sizes = sizes;
          end
        
        function obj = Plot3DBoundary(obj)
            %vykresli obrys mozku ve vsech rozmerech do 3d grafu
            %pokud boundary neni vypocitana, spocita ji a ulozi do obj.plotCh3D.BrainBoundaryXYZ
            %netvori graf, kresli do existujiciho a aktivniho ChannelPlot, a predpoklada hold on;
            dimenze = [2 3; 1 3; 2 1]; %xy, xz a yz 
            load('GMSurfaceMesh.mat');             %#ok<LOAD>
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
        
        function obj = hybejPlot3D(obj,~,eventDat)
              switch eventDat.Key
                  case 's'    %sagital view                                       
                      obj.plotCh3D.view =  iff( isequal( [90 0],obj.plotCh3D.view),[-90 0],[90 0]); %zprava,zleva  
                      view(obj.plotCh3D.view); 
                  case {'c','f'} %coronal = predozadni, frontal                                            
                      obj.plotCh3D.view =  iff( isequal( [0 0],obj.plotCh3D.view),[180 0],[0 0]); %zezadu, zpredu 
                      view(obj.plotCh3D.view); %zleva
                  case {'h','a'} %horizontal = hornodolni nebo axial                        
                      obj.plotCh3D.view =  iff( isequal( [0 90],obj.plotCh3D.view),[ 180 -90],[0 90]); %shora, zdola
                      view(obj.plotCh3D.view); 
                  case 'space'
                     if isfield(obj.plotCh3D,'boundary') %prepinam v grafu cely scatter s jen hranici mozku - hlavne kvuli kopirovani do corelu
                         obj.plotCh3D.boundary  = 1 - obj.plotCh3D.boundary;
                     else
                         obj.plotCh3D.boundary  = 1;
                     end
                     obj.ChannelPlot3D();
                  case 'n' %names
                     if ~isempty(eventDat.Modifier) && strcmp(eventDat.Modifier{:},'shift') 
                          obj.plotCh3D.names =0; %by the shift+n, switch off all names
                          obj.plotCh3D.labesXnames = 0; %just switch to chnames
                     else
                         if obj.plotCh3D.labesXnames %if brainlabels are plotted
                             obj.plotCh3D.labesXnames = 0; %just switch to chnames
                         else
                             obj.plotCh3D.names  = obj.plotCh3D.names + 1; %switch of channel labels - 0=nothing,1=channel no,2=channel name,3=neurologyLabels,4=pacient name
                             if obj.plotCh3D.names > 4, obj.plotCh3D.names =0; end                                             
                         end
                     end
                     obj.ChannelPlot3D();
                  case 'b' %show brain labels instead of channel names  
                     if ~obj.plotCh3D.labesXnames %if chnnames are plotted
                         obj.plotCh3D.labesXnames = 1; %just switch to brainlabels
                     else
                         obj.plotCh3D.labels = obj.plotCh3D.labels + 1; 
                         if obj.plotCh3D.labels >3, obj.plotCh3D.labels =0; end %switch of channels labels 0=nothing,1=lobe,2=brainlabel,3=class
                         obj.plotCh3D.names = 0; %switch off labeling channels by channel names
                     end
                     obj.ChannelPlot3D();
                  case 'r'
                     %dialog na vlozeni souradnic roi hodnoty
                    answ = inputdlg('Enter x,y,z a edge size:','define ROIs', [10 50],{num2str(obj.plotCh3D.roi)});
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
                    obj.ChannelPlot3D();
                  case 'p' %show all points in file, excluding the rejected, as points                                       
                    obj.plotCh3D.allpoints  = obj.plotCh3D.allpoints + 1;
                    if obj.plotCh3D.allpoints > 2 %values 0 1 and 2
                        obj.plotCh3D.allpoints = 0; %reset the value                  
                    elseif obj.plotCh3D.allpoints == 2 % this value means to show even channel labels for all points
                        obj.plotCh3D.allpointnames = 1; 
                    else
                        obj.plotCh3D.allpointnames=0;
                    end
                    obj.ChannelPlot3D();
                  case 'z' %switch of zoom 0 1(proportional) and 2 (non-proportional)                    
                    obj.plotCh3D.zoom  = obj.plotCh3D.zoom + 1;
                    if obj.plotCh3D.zoom > 2, obj.plotCh3D.zoom = 0; end                    
                    obj.ChannelPlot3D();
                 case 'w' %zapinani a vypinani prazdneho pozadi obrazku
                    if isfield(obj.plotCh3D,'background') 
                        obj.plotCh3D.background = 1-obj.plotCh3D.background;
                    else
                        obj.plotCh3D.background = 0;
                    end
                    obj.ChannelPlot3D();
                  case 'o' %reorder channels, so that highest vals will be in front
                    obj.plotCh3D.reorder = 1-obj.plotCh3D.reorder;  
                    obj.ChannelPlot3D();
                  case 'l' %lines - spojnice kanalu - zadne, pacienti
                    %zobrazeni pozic vsech kanalu jako tecek
                    if isfield(obj.plotCh3D,'lines') 
                       obj.plotCh3D.lines  = 1 - obj.plotCh3D.lines;
                    else
                       obj.plotCh3D.lines  = 1;
                    end
                    obj.ChannelPlot3D();
                  case 'd' %just reDraw the plot
                    obj.ChannelPlot3D();  
                  case {'add' ,  'equal'} %increase the font size of channel labels
                    obj.plotCh3D.fontsize = obj.plotCh3D.fontsize + 1;  
                    obj.ChannelPlot3D();  
                  case {'subtract' , 'hyphen'}    %decrease the font size of channel labels
                    obj.plotCh3D.fontsize = obj.plotCh3D.fontsize - 1; 
                    obj.ChannelPlot3D();  
                  case 'v' %toggle color code of 3D scatter                    
                    iCluster = obj.CH.GetCluster(obj.plotCh3D.popis);
                    colorusemax = iff(iCluster,4,3); %when there are the clusters computed, show them when plotCh3D.coloruse = 4
                    obj.plotCh3D.coloruse = obj.plotCh3D.coloruse + 1;
                    if obj.plotCh3D.coloruse >= colorusemax, obj.plotCh3D.coloruse=0; end %values only 0 1 or 2%                   
                    obj.ChannelPlot3D(); 
                  case 'm' %toggle markertypes to use for scatter 3D
                   markertypes = {'o','s','d','p'};
                   im = find(contains(markertypes,obj.plotCh3D.markertype)) +1;
                   if im > numel(markertypes), im = 1; end
                   obj.plotCh3D.markertype = markertypes{im}; 
                   obj.ChannelPlot3D(); 
                  case {'divide','slash'} %slash on numerical keyboard - automatic range of colorbar
                    obj.plotCh3D.rangeZ = [min(obj.plotCh3D.chnvals) max(obj.plotCh3D.chnvals)]; 
                    obj.ChannelPlot3D(); %plot the figure again
                  case 'q' %show / hides the clusters
                    obj.plotCh3D.showclusters = 1 - obj.plotCh3D.showclusters;
                    obj.ChannelPlot3D(); %plot the figure again   
              end
        end

        function highlightChannel(obj, ch)
          %channel is absolute channel number
          if ~exist('ch','var') || isempty(ch)
              if isfield(obj.plotCh3D,'ch_highlighted')
                ch = obj.plotCh3D.ch_highlighted;
              else
                ch = 0; %no highlighted channel
              end
          end
             
          if ch && isfield(obj.plotCh3D,'fh') && ishandle(obj.plotCh3D.fh) % pokud mam otevreny plot
            ax = obj.plotCh3D.fh.CurrentAxes;
            %disp(displayedChannels(closestChannel).name)
            x = obj.CH.H.channels(ch).MNI_x; y = obj.CH.H.channels(ch).MNI_y; z = obj.CH.H.channels(ch).MNI_z;
            if isfield(obj.plotCh3D, 'selHandle') % smazu predchozi oznaceni, pokud nejake bylo
                delete(obj.plotCh3D.selHandle)
            end
            if isfield(obj.plotCh3D, 'selNameHandle') % smazu predchozi oznaceni, pokud nejake bylo
                delete(obj.plotCh3D.selNameHandle)
            end
            obj.plotCh3D.selHandle = scatter3(ax, x, y, z, obj.plotCh3D.sizes(obj.plotCh3D.chnsel==ch)+60, 'ok', 'LineWidth', 2); % oznacim vybrany kanal na 3D grafu
            obj.plotCh3D.selNameHandle = annotation(obj.plotCh3D.fh, 'textbox',[0 1 0 0],'String',obj.CH.H.channels(ch).name,'FitBoxToText','on');
          else  % pokud se zadny kanal nenasel (kliknuti mimo)
             if isfield(obj.plotCh3D, 'selHandle') % smazu predchozi oznaceni, pokud nejake bylo
               delete(obj.plotCh3D.selHandle)
             end
             if isfield(obj.plotCh3D, 'selNameHandle') % smazu predchozi oznaceni, pokud nejake bylo
                delete(obj.plotCh3D.selNameHandle)
            end
          end
           obj.plotCh3D.ch_highlighted = ch;
        end
        
        function delete(obj) %destructor of a handle class
            if isfield(obj.plotCh3D,'fh') && ~isempty(obj.plotCh3D.fh) && ishandle(obj.plotCh3D.fh) 
                close(obj.plotCh3D.fh); 
            end
        end
        
    end
end