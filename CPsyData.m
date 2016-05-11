classdef CPsyData
    %CPSYDATA Trida na praci s behavioralnimi daty z psychopy
    %   vytvorena pomoci ppa_data.m, aedist_data.m aj
    
    properties  (Access = public)
        P; %psychopy behavioural data
    end
    
    methods (Access = public)
        function obj = CPsyData(psy)
            %konstruktor
            obj.P = psy;
        end
        function rt = ReactionTime(obj)
            %vrati matici vsech reakcnich casu roztridenych do sloupcu podle kategorii podnetu
            %kde je min hodnot v kategorii, zarovnano pomoci NaN, takze se musi pouzit nanmean
            kat = unique(obj.P.data(:,obj.P.sloupce.kategorie)); %ciselne vyjadreni kategorie podnetu 0-n
            rt = nan(size(obj.P.data,1),numel(kat)); %pole kam budu ukladat reakcni casy v sekundach
            
            for k = kat' %funguje jen pro radky, kat je sloupec
                ikat = obj.P.data(:,obj.P.sloupce.kategorie)==k;
                rt(1:sum(ikat),k+1) = 24*3600*(obj.P.data(ikat,obj.P.sloupce.ts_odpoved) - obj.P.data(ikat,obj.P.sloupce.ts_podnet));
            end
            rt = rt(any(~isnan(rt),2),:); % necha jen radky, kde je nejake ~NaN cislo           
        end
        
        function isi = InterStimulusInterval(obj)
            %vraci matici vsech interstimulus intervalu v sekundach
            isi = 24*3600*(obj.P.data(2:end,obj.P.sloupce.ts_podnet) - obj.P.data(1:end-1,obj.P.sloupce.ts_podnet));
        end
        
        function ts_podnety = TimePodnety(obj)
            %vraci matici timestampu vsech podnetu
            ts_podnety = obj.P.data(:,obj.P.sloupce.ts_podnet);
        end
        function [kat katnum] = Kategorie(obj,event)
            %vrati retezec s popisem kategorie eventu
            katnum = obj.P.data(event,obj.P.sloupce.kategorie); %cislo kategorie
            kat = obj.P.strings.podminka{katnum+1};            
        end    
            
    end
    
end

