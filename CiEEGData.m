classdef CiEEGData < handle
    %CEEGDATA Trida na praci s datama ve formatu ISARG od Petra Jezdika
    %   Kamil Vlcek, 2016-04
    
    properties (Access = public)
        d; %double nebo int matrix: time x channel, muze byt i time x channel x epoch
        tabs; 
        tabs_orig; %originalni tabs, ktere se zachovaji po epochaci. Downsamplovani se u nich dela
        fs; %vzorkovaci frekvence
        mults; %nepovinne
        header; %nepovinne
        
        samples; %pocet vzorku v zaznamu = rozmer v case
        channels; %pocet kanalu v zaznamu
        epochs;   %pocet epoch
        epochData; %cell array informaci o epochach; epochy v radcich, sloupce: kategorie, tab
        PsyData; %objekt ve formatu CPsyData (PPA, AEDist aj) podle prezentace KISARG
            %pole PsyData.P.data, sloupce, strings, interval, eegfile, pacientid
        epochtime; %delka eventu pre a po event v sekundach       
    end
    
    methods (Access = public)
        function obj = CiEEGData(d,tabs,fs,mults,header)
            %konstruktor         
            obj.d = d;
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
            obj.tabs = tabs;
            obj.tabs_orig = tabs;
            obj.fs = fs;
            if exist('mults','var') && ~isempty(mults)
                obj.mults = mults;
            else
                obj.mults = ones(1,size(d,2)); %defaultove jednicky pro kazdy kanal
            end
            if exist('header','var')
                obj.header = header;
            else
                obj.header = [];
            end
            
        end
        
        function [samples, channels, epochs] = DSize(obj)
            % vraci velikosti pole d
            samples = size(obj.d,1);
            channels = size(obj.d,2);
            epochs = size(obj.d,3);
        end
        
        function dd = DData(obj,ch,epoch)
            % vraci jeden kanal a jednu epochu ze zaznamu. Implementovano kvuli nasobeni mults
            if ~exist('epoch','var')
                epoch = 1;            
            end
            dd = double(obj.d(:,ch,epoch)) .* obj.mults(ch);                           
        end
        
        function ExtractEpochs(obj, PsyData,epochtime)
            % epochuje data v poli d, pridava do objektu:
            % cell array epochData, double(2) epochtime v sekundach, tridu PsyData 
            % upravuje obj.mults, samples channels epochs
            if obj.epochs > 1
                disp('already epoched data');
                return;
            end
            obj.PsyData = PsyData; %objekt CPsyData
            obj.epochtime = epochtime; %v sekundach cas pred a po udalosti  , prvni cislo je zaporne druhe kladne
            iepochtime = round(epochtime.*obj.fs); %v poctu vzorku cas pred a po udalosti
            ts_podnety = PsyData.TimeStimuli(); %timestampy vsech podnetu
            de = zeros(iepochtime(2)-iepochtime(1), size(obj.d,2), size(ts_podnety,1)); %nova epochovana data time x channel x epoch            
            tabs = zeros(iepochtime(2)-iepochtime(1),size(ts_podnety,1)); %#ok<PROP> %udelam epochovane tabs
            obj.epochData = cell(size(ts_podnety,1),3); % sloupce kategorie, cislo kategorie, timestamp
            for epoch = 1:size(ts_podnety,1) %pro vsechny eventy
                izacatek = find(obj.tabs<=ts_podnety(epoch), 1, 'last' ); %najdu index podnetu podle jeho timestampu
                    %kvuli downsamplovani Hilberta, kdy se mi muze ztratit presny cas zacatku
                    %epochy, beru posledni nizsi tabs nez je cas zacatku epochy
                [Kstring Knum] = PsyData.Category(epoch);    %jmeno a cislo kategorie
                obj.epochData(epoch,:)= {Kstring Knum obj.tabs(izacatek)}; %zacatek epochy beru z tabs aby sedel na tabs pri downsamplovani
                for ch = 1:obj.channels %pro vsechny kanaly                    
                    baseline = mean(obj.d(izacatek+iepochtime(1):izacatek-1));
                    de(:,ch,epoch) = double(obj.d( izacatek+iepochtime(1) : izacatek+iepochtime(2)-1,ch)).* obj.mults(ch) - baseline; 
                    tabs(:,epoch) = obj.tabs(izacatek+iepochtime(1) : izacatek+iepochtime(2)-1); %#ok<PROP>
                end
            end
            obj.d = de; %puvodni neepochovana budou epochovana
            obj.mults = ones(1,size(obj.d,2)); %nove pole uz je double defaultove jednicky pro kazdy kanal
            obj.tabs = tabs; %#ok<PROP>
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
        end
        function [d]= CategoryData(obj, katnum)
            %vraci epochy ve kterych podnet byl kategorie/podminky katnum
            assert(obj.epochs > 1,'data not yet epoched'); %vyhodi chybu pokud data nejsou epochovana
            iEpochy = cell2mat(obj.epochData(:,2))==katnum ; %seznam epoch v ramci kategorie ve sloupci
            d = obj.d(:,:,iEpochy);
        end
        function PlotCategory(obj,katnum,channel)
            d1=obj.CategoryData(katnum);
            d1m = mean(d1,3);
            T = (0 : 1/obj.fs : (size(obj.d,1)-1)/obj.fs) + obj.epochtime(1); %cas zacatku a konce epochy
            E = 1:obj.epochs; %vystupni parametr
            h1 = figure('Name','Mean Epoch');
            plot(T,d1m(:,channel));
            xlabel('Time [s]'); 
            title(obj.PsyData.CategoryName(katnum));
            h2 = figure('Name','All Epochs');            
            imagesc(T,E,squeeze(d1(:,channel,:))');
            colorbar;
            xlabel('Time [s]');
            ylabel('Epochs');
            title(obj.PsyData.CategoryName(katnum));
            
        end
    end
    
end

