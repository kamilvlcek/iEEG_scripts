classdef CiEEGData < handle
    %CEEGDATA Trida na praci s datama ve formatu ISARG od Petra Jezdika
    %   Kamil Vlcek, 2016-04
    
    properties (Access = public)
        d; %double nebo int matrix: time x channel, muze byt i time x channel x epoch
        tabs; 
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
            obj.fs = fs;
            if exist('mults','var') && ~isempty(mults)
                obj.mults = mults;
            else
                obj.mults = ones(size(d,2)); %defaultove jednicky pro kazdy kanal
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
            % cell array epochData, double(2) epochtime, tridu PsyData 
            % upravuje obj.mults, samples channels epochs
            obj.PsyData = PsyData; %objekt CPsyData
            obj.epochtime = epochtime; %v sekundach cas pred a po udalosti  , prvni cislo je zaporne druhe kladne
            iepochtime = round(epochtime.*obj.fs); %v poctu vzorku cas pred a po udalosti
            ts_podnety = PsyData.TimePodnety(); %timestampy vsech podnetu
            de = zeros(iepochtime(2)-iepochtime(1), size(obj.d,2), size(ts_podnety,1)); %nova epochovana data time x channel x epoch            
            obj.epochData = cell(size(ts_podnety,1),2); % sloupce kategorie, timestamp
            for epoch = 1:size(ts_podnety,1) %pro vsechny eventy
                izacatek = find(obj.tabs==ts_podnety(epoch)); %najdu index podnetu podle jeho timestampu
                obj.epochData(epoch,:)= {PsyData.Kategorie(epoch) ts_podnety(epoch)};
                for ch = 1:obj.channels %pro vsechny kanaly                    
                    baseline = mean(obj.d(izacatek+iepochtime(1):izacatek-1));
                    de(:,ch,epoch) = double(obj.d( izacatek+iepochtime(1) : izacatek+iepochtime(2)-1,ch)).* obj.mults(ch) - baseline; 
                end
            end
            obj.d = de; %puvodni neepochovana budou epochovana
            obj.mults = ones(size(obj.d,2)); %nove pole uz je double defaultove jednicky pro kazdy kanal
            [obj.samples,obj.channels, obj.epochs] = obj.DSize();
        end
    end
    
end

