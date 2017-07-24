classdef CHHeader < handle
    %CHHEADER Trida na praci s headerem od Jirky Hammera
    %Kamil Vlcek, FGU AVCR, since 2016 04
    
    properties (Access = public)
        H; %data od Jirky Hammera
        %E; % unikatni jmena elektrod, napr RG, I, GC ...
        chgroups;
        els;
        RjCh; 
    end
    
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
        function obj = ChannelPlot(obj,pohled)
            %zobrazi 3D obrazek elektrod v MNI prostoru. Obrazek ma rozmery podle rozmeru mozku
            %pohled muze urcti smer pohledu s-sagital,c-coronal,h-horizontal
            if isfield(obj.H.channels,'MNI_x')
                figure('Name','ChannelPlot in MNI');                
                [obj,chgroups] = obj.ChannelGroups(); %#ok<PROP>
                for chg = 1:size(chgroups,2) %#ok<PROP>
                    group = chgroups{chg}; %#ok<PROP>
                    X = zeros(1,numel(group)); Y = X; Z = X;
                    for ich = 1:numel(group)                        
                        X(ich) = obj.H.channels(group(ich)).MNI_x;
                        Y(ich) = obj.H.channels(group(ich)).MNI_y;
                        Z(ich) = obj.H.channels(group(ich)).MNI_z;                        
                    end
                    plot3(X,Y,Z,'o-','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',2);
                    if chg==1, hold on; end                    
                    for ich = 1:numel(group) 
                        text(X(ich),Y(ich),Z(ich)+3,obj.H.channels(group(ich)).name);
                    end
                end
                xlabel('MNI X'); %levoprava souradnice
                ylabel('MNI Y'); %predozadni souradnice
                zlabel('MNI Z'); %hornodolni
                if ~exist('pohled','var'), pohled = ''; end
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
            else
                disp('No MNI data');
            end
        end
        function tag= PacientTag(obj)
            %vraci tag pacienta, napriklad p73
            if isfield(obj.H,'patientTag'), tag = obj.H.patientTag; else tag=obj.H.subjName; end
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
        function [names]=GetChNames(obj,channels)
            %vrati jmena vsech kanalu jako cell array
            if ~exist('channels','var'), channels = 1:size(obj.H.channels,2); end
            names = cell(numel(channels),1);
            for ch = 1: numel(channels)
                if strcmpi(obj.H.channels(ch).signalType,'SEEG')
                    names{ch}=obj.H.channels(ch).name;
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
            else
                obj.H.channels(tch).signalType = 'SEEG'; %nebo nastavi ze to naopak neni trigger
            end
            tch = find(strcmp({obj.H.channels.signalType}, 'triggerCh')==1);
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
    end
    
end

