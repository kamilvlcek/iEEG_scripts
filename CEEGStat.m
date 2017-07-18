classdef CEEGStat
    %CEEGSTAT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
       d; %matice EEG dat 
       fs;
    end
    
    methods (Access = public)
        
        function obj = CEEGStat(d,fs)
            %konstruktor
            obj.d = d;            
            obj.fs = fs;
        end
        
        function [P,ibaseline,iepochtime] = WilcoxBaseline(obj,epochtime,baseline,timewindow,iEp)
            % spocita signifikance vuci baseline
            iepochtime = round(epochtime(1:2).*obj.fs); %v poctu vzorku cas pred a po udalosti, pred je zaporne cislo           
            ibaseline = round(baseline.*obj.fs); %zaporna cisla pokud pred synchro eventem
            ibaselineA =  [ibaseline(1) ,    floor((ibaseline(2)-ibaseline(1))/2)+ibaseline(1)]; %prvni pulka baseline - 28.3.2017
            ibaselineB =  [ibaselineA(2)+1 , ibaseline(2) ]; %druha pulka baseline
            itimewindow = round(timewindow.*obj.fs); %
            
            %proc driv?  mean(obj.d( abs(iepochtime(1)-ibaselineA(1))+1 : abs(iepochtime(1)-ibaselineA(2))  - 28.3.2017            
            baselineA = mean(obj.d( ibaselineA(1)-iepochtime(1)+1 : ibaselineA(2)-iepochtime(1) , : , iEp),1); 
            baselineB = mean(obj.d( ibaselineB(1)-iepochtime(1) : ibaselineB(2)-iepochtime(1) , : , iEp),1); 
                % cas x kanaly x epochy - prumer za cas pred podnetem, pro vsechny kanaly a nevyrazene epochyt
            if numel(itimewindow) == 2  %chci prumernou hodnotu eeg odpovedi d od do
                response = mean(obj.d( (itimewindow(1) : itimewindow(2)) - iepochtime(1) , : , iEp ),1); %prumer v case                
            else %time windows je delka maximalni hodnoty p - takze tady jeste neresim
                response = obj.d(abs(iepochtime(1)-ibaseline(2))+1 : end , : , iEp ); %vsechny hodnoty po konci baseline                 
            end
            WpA = CStat.Wilcox2D(response,baselineA,1,[],'mean vs baseline A');  %1=mene striktni pdep, 2=striktnejsi dep;
            WpB = CStat.Wilcox2D(response,baselineB,1,[],'mean vs baseline B');  %1=mene striktni pdep, 2=striktnejsi dep;
            Wp = max(WpA,WpB); % %vyssi hodnota z kazde poloviny baseline
            if numel(itimewindow) == 1 %chci maximalni hodnotu p z casoveho okna
                P = CStat.Klouzaveokno(Wp,itimewindow(1),'max',1); % %pole 2D signifikanci si ulozim kvuli kresleni - cas x channels                
            else
                P = Wp; %%pole 1D signifikanci - jedna hodnota pro kazdy kanal                
            end
        end
    end
    
end

