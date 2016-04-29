classdef CiEEGData < handle
    %CEEGDATA Trida na praci s datama ve formatu ISARG od Petra Jezdika
    %   Kamil Vlcek, 2016-04
    
    properties (Access = public)
        d; %double nebo int matrix: time x channel, muze byt i time x channel x epoch
        tabs; 
        fs; %vzorkovaci frekvence
        mults; %nepovinne
        header; %nepovinne
        
        samples; %pocet vzorku v zaznamu
        channels; %pocet kanalu v zaznamu
        epochs;  %pocet epoch v datech
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
    end
    
end

