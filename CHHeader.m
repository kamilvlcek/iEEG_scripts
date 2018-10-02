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
        plotCh2D; %udaje o 2D grafu kanalu, hlavne handle
    end
    %#ok<*PROPLC>
    methods (Access = public)
        function obj = CHHeader(H)
            %konstruktor
            obj.H = H;     
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
        
        function [obj, chgroups, els] = ChannelGroups(obj)
            %vraci skupiny kanalu (cisla vsech channels na elekrode) + cisla nejvyssiho kanalu v na kazde elektrode v poli els
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
        function [XYZ,obj] = ChannelPlot(obj,pohled,labels,XYZ2)
            %zobrazi 3D obrazek elektrod v MNI prostoru. Obrazek ma rozmery podle rozmeru mozku
            %pohled muze urcti smer pohledu s-sagital,c-coronal,h-horizontal
            if ~exist('pohled','var') || isempty(pohled), pohled = ''; end
            if ~exist('labels','var') || isempty(labels), labels = 0; end
            if isfield(obj.H.channels,'MNI_x')
                figure('Name','ChannelPlot in MNI');                
                [obj,chgroups] = obj.ChannelGroups();          
                
                %objekt se dobre uklada i pri poradi return values XYZ,obj
                XYZ = struct('X',0,'Y',0,'Z',0);
                for chg = 1:size(chgroups,2) 
                    group = chgroups{chg}; 
                    X = zeros(1,numel(group)); Y = X; Z = X;
                    for ich = 1:numel(group)                        
                        X(ich) = obj.H.channels(group(ich)).MNI_x;
                        Y(ich) = obj.H.channels(group(ich)).MNI_y;
                        Z(ich) = obj.H.channels(group(ich)).MNI_z;                        
                    end
                    XYZ(chg) = struct('X',X,'Y',Y,'Z',Z);
                    plot3(X,Y,Z,'o-','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',2);
                    if chg==1, hold on; end                    
                    for ich = 1:numel(group) 
                        if labels
                            th = text(X(ich),Y(ich),Z(ich)+3,obj.H.channels(group(ich)).neurologyLabel);
                            th.FontSize = 8;
                        else
                            th = text(X(ich),Y(ich),Z(ich)+3,obj.H.channels(group(ich)).name);
                        end
                    end
                end
                if exist('XYZ2','var')
                    for chg = 1:numel(XYZ2)
                        plot3(XYZ2(chg).X,XYZ2(chg).Y,XYZ2(chg).Z,'o-','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',2);    
                    end
                end
                xlabel('MNI X'); %levoprava souradnice
                ylabel('MNI Y'); %predozadni souradnice
                zlabel('MNI Z'); %hornodolni
                
                switch pohled
                    case 's' %sagital = levoprava
                        view([-1 0 0]); %zleva
                    case 'c' %coronal = predozadni
                        view([0 1 0]); %zepredu
                    case 'h' %horizontal = hornodolni   
                        view([0 0 1]); %shora
                end
                axis([-75 75 -120 70 -50 90]); %zhruba velikost mozku        
                text(-70,0,0,'LEVA');        
                text(70,0,0,'PRAVA');   
                text(0,65,0,'VPREDU');        
                text(0,-115,0,'VZADU');                 
                load('GMSurfaceMesh.mat'); %seda hmota v MNI
                scatter3(GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3),'.','MarkerEdgeAlpha',.2);
                %load('WMSurfaceMesh.mat');
                %scatter3(WMSurfaceMesh.node(:,1),WMSurfaceMesh.node(:,2),WMSurfaceMesh.node(:,3),'.','MarkerEdgeAlpha',.1);
            else
                disp('No MNI data');
            end
        end
        function ChannelPlot2D(obj,chsel,selCh,plotChH,label)
            %vstupni promenne
            if ~exist('chsel','var')%promenna na jeden cerveny kanal
                if isfield(obj.plotCh2D,'chsel')
                    chsel = obj.plotCh2D.chsel;
                else
                    chsel = 1;  
                    obj.plotCh2D.chsel = 1;
                end
            else
                obj.plotCh2D.chsel = chsel;
            end
            chsel = obj.sortorder(chsel);
            if ~exist('selCh','var')%promenna na vic cernych kanaly, pro obj.PlotRCh.SelCh
                if isfield(obj.plotCh2D,'selCh')
                    selCh = obj.plotCh2D.selCh;
                else
                    selCh = []; 
                    obj.plotCh2D.selCh = [];
                end
            else
                obj.plotCh2D.selCh = selCh;
            end
            if exist('plotChH','var')  %handlet na funkci z CiEEGData @obj.PlotResponseCh
                obj.plotCh2D.plotChH = plotChH;
            end
            if ~isfield(obj.plotCh2D,'marks')  %handlet na funkci z CiEEGData @obj.PlotResponseCh
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
            
            %vytvoreni figure
            x = [obj.H.channels(:).MNI_x];
            y = [obj.H.channels(:).MNI_y];
            z = [obj.H.channels(:).MNI_z];            
            load('GMSurfaceMesh.mat'); %seda hmota v MNI
            if isfield(obj.plotCh2D,'boundary') && obj.plotCh2D.boundary && ~isfield(obj.plotCh2D,'BrainBoundaryXY') %trva docela dlouho nez se to spocita
                obj.plotCh2D.BrainBoundaryXY = boundary(GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2)); %vnejsi hranice mozku
                obj.plotCh2D.BrainBoundaryYZ = boundary(GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3));
            end
            
            size_ch = 10; %velikosti krouzko oznacujicich kanaly
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
            %axialni plot            
            if isfield(obj.plotCh2D,'boundary') && obj.plotCh2D.boundary
                %defaultne budu vykreslovat scatter, ale kvuli kopirovani se bude hodit i jen boundary
                plot(GMSurfaceMesh.node(obj.plotCh2D.BrainBoundaryXY,1),GMSurfaceMesh.node(obj.plotCh2D.BrainBoundaryXY,2));
            else
                scatter(GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2),'.','MarkerEdgeAlpha',.1); %seda hmota normalizovaneho mozku
            end           
            hold on;             
            
            for ie = 1:numel(els)                
                plot(x(els0(ie):els(ie)),y(els0(ie):els(ie)),'-o'); %plot kontaktu jedne elektrody
                for ch = els0(ie):els(ie)
                    th = text(x(ch),y(ch),num2str(ch)); %cislo kazdeho kanalu
                    th.FontSize = 8;
                end                              
            end
            if ~isempty(chsel) %pokud je vybrany nejaky kanal
                plot(x(chsel),y(chsel),'o','MarkerSize',size_ch,'MarkerEdgeColor','r','MarkerFaceColor','r'); 
                chstr = iff(isempty(obj.sortedby),num2str(chsel), [ num2str(obj.sortorder(chsel)) '(' obj.sortedby  num2str(chsel) ')' ]);
                title( [ 'channel ' chstr ]);
                
            end
            if ~isempty(selCh) %hromadne vybrane kanaly, zobrazne cernym koleckem
                barvy = 'bgcmky';
                for m = 1:size(selCh,2) %jednu znacku za druhou
                   if  obj.plotCh2D.marks(m) %pokud se ma znacka zobrazovat
                       ch = find(selCh(:,m)); %seznam cisel vybranych kanalu pro danou znacku
                       ch = intersect(chshow,ch); 
                       plot(x(ch),y(ch),'o','MarkerSize',size_selCh,'MarkerEdgeColor',barvy(m),'MarkerFaceColor',barvy(m));
                   end
                end
            end
            text(-70,70,'LEVA');
            text(55,70,'PRAVA');  
            axis equal;  
            if isfield(obj.plotCh2D,'grid') && obj.plotCh2D.grid==1
                grid on;
            end
                
            xticks(-70:10:70);
            yticks(-100:10:70);
            xlabel('MNI X'); %levoprava souradnice
            ylabel('MNI Y'); %predozadni souradnice
            if isfield(obj.plotCh2D,'background') && obj.plotCh2D.background==0
                set(gca,'color','none'); %zadne bile pozadi, pak ani v corelu
            end
           
            subplot(1,2,2);
            %sagitalni plot            
            if isfield(obj.plotCh2D,'boundary') && obj.plotCh2D.boundary
                plot(GMSurfaceMesh.node(obj.plotCh2D.BrainBoundaryYZ,2),GMSurfaceMesh.node(obj.plotCh2D.BrainBoundaryYZ,3));
            else
                scatter(GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3),'.','MarkerEdgeAlpha',.1);   %seda hmota normalizovaneho mozku
            end
            hold on;     
            for ie = 1:numel(els)                 
                plot(y(els0(ie):els(ie)),z(els0(ie):els(ie)),'-o'); %plot kontaktu jedne elektrody
                for ch = els0(ie):els(ie)
                    th = text(y(ch),z(ch),num2str(ch));
                    th.FontSize = 8;
                end                
            end  
            if ~isempty(chsel) %pokud je vybrany nejaky kanal
                plot(y(chsel),z(chsel),'o','MarkerSize',size_ch,'MarkerEdgeColor','r','MarkerFaceColor','r'); 
                
                text(x_text,110,[ obj.H.channels(1,chsel).name]);
                text(x_text,100,[ obj.H.channels(1,chsel).neurologyLabel ',' obj.H.channels(1,chsel).ass_brainAtlas]);
                if  isfield(obj.H.channels,'MNI_x') %vypisu MNI souradnice
                    text(x_text,90,[ 'MNI: ' num2str(round(obj.H.channels(1,chsel).MNI_x)) ', ' num2str(round(obj.H.channels(1,chsel).MNI_y )) ', ' num2str(round(obj.H.channels(1,chsel).MNI_z))]);
                else
                    text(x_text,90,'no MNI');
                end                
            end
            if ~isempty(selCh) %hromadne vybrane kanaly, zobrazne cernym koleckem                
                barvy = 'bgcmky';
                klavesy = 'fghjkl'; %abych mohl vypsat primo nazvy klaves vedle hvezdicky podle selCh
                for m = 1:size(selCh,2) %jednu znacku za druhou
                    if  obj.plotCh2D.marks(m) %pokud se ma znacka zobrazovat
                       ch = find(selCh(:,m)); %seznam cisel vybranych kanalu pro danou znacku
                       ch = intersect(chshow,ch); 
                       if ~isempty(ch) %pokud jsou takove nejake vybrane kanaly
                           plot(y(ch),z(ch),'o','MarkerSize',size_selCh,'MarkerEdgeColor',barvy(m),'MarkerFaceColor',barvy(m));
                           th = text(x_text+m*10,-90,klavesy(m), 'FontSize', 15,'Color',barvy(m)); %legenda k barvam kanalu dole pod mozkem
                           th.BackgroundColor = [.6 .6 .6];
                       end
                    end
                end
                if any(selCh(chsel,:),2)==1 %pokud je aktualni kanal jeden z vybranych                
                    klavesy = 'fghjkl'; %abych mohl vypsat primo nazvy klaves vedle hvezdicky podle selCh
                    text(x_text,80,['*' klavesy(logical(selCh(chsel,:)))], 'FontSize', 12,'Color','red');
                end
                if ~isempty(label)
                    text(x_text,-75,strrep(label,'_','\_'), 'FontSize', 10,'Color','blue' );
                end
                if isfield(obj.plotCh2D,'chshowstr') && ~isempty(obj.plotCh2D.chshowstr)
                    text(0,-90,['chshow:' obj.plotCh2D.chshowstr] ,'Color','red');
                end
            end
            axis equal;
            if isfield(obj.plotCh2D,'grid') && obj.plotCh2D.grid==1
                grid on;
            end
            yticks(-80:10:80);
            xticks(-100:10:70);
            xlabel('MNI Y'); %predozadni souradnice
            ylabel('MNI Z'); %hornodolni
            if isfield(obj.plotCh2D,'background') && obj.plotCh2D.background==0
                set(gca,'color','none'); %zadne bile pozadi, pak ani v corelu
            end
            
            %rozhybani obrazku            
            set(obj.plotCh2D.fh,'KeyPressFcn',@obj.hybejPlot2D); 
            set(obj.plotCh2D.fh, 'WindowButtonDownFcn', @obj.hybejPlot2Dclick);
        end
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
            if ~exist('by','var')
                obj.sortorder = 1:numel(obj.H.channels);
                obj.sortedby = '';
            elseif strcmp(by,'x')
                [~,obj.sortorder]=sort([obj.H.channels(:).MNI_x]);                
                obj.sortedby = 'x';
            elseif strcmp(by,'y')
                [~,obj.sortorder]=sort([obj.H.channels(:).MNI_y]);                
                obj.sortedby = 'y';
            elseif strcmp(by,'z')
                [~,obj.sortorder]=sort([obj.H.channels(:).MNI_z]); 
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
                iL = contains(ChLabels,chlabels);
                if exist('notchnlabels','var') && numel(notchnlabels) > 0
                    iLx = contains(ChLabels,notchnlabels);
                    iL = iL & ~iLx;
                    obj.plotCh2D.chshowstr = [ cell2str(chlabels) ' not:' cell2str(notchnlabels)];
                else
                    obj.plotCh2D.chshowstr = cell2str(chlabels);
                end
                obj.plotCh2D.chshow = find(iL); %vyber kanalu k zobrazeni       
            else
                obj.plotCh2D.chshow = 1:numel(obj.H.channels);
                obj.plotCh2D.chshowstr = '';
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
              if subp == 0  %axialni graf, subplot 1
                  if x > -70 && x < 70 && y > -120 && y < 90
                    chns_mni = [[obj.H.channels(:).MNI_x]' , [obj.H.channels(:).MNI_y]'];                  
                  end
              else         %sagitalni graf, subplot 2
                  if x > -100 && x < 70 && y > -100 && y < 90
                    chns_mni = [[obj.H.channels(:).MNI_y]' , [obj.H.channels(:).MNI_z]'];   
                  end
              end
              if ~isempty(chns_mni)
                  [ch,d] = dsearchn(chns_mni,[x y]); %najde nejblizsi kanal a vzdalenost k nemu                   
                  obj.ChannelPlot2D(find(obj.sortorder==ch)); %#ok<FNDSB>   
                  if isfield(obj.plotCh2D,'plotChH')
                    obj.plotCh2D.plotChH(obj.plotCh2D.chsel); %vykreslim @obj.PlotResponseCh  
                    figure(obj.plotCh2D.fh); %dam puvodni obrazek dopredu
                  end
                  
              end
          end
    end
    
end

