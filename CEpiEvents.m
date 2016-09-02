classdef CEpiEvents < handle
    %CEPIEVENTS Trida na praci s epileptickymi eventy detekovanymi detektorem ISARG
    %Kamil Vlcek, FGU AVCR, since 2016 09   
    
    properties (Access = public)
        DE; %matice eventu pro snadne vyhledavani
        sloupce; %nazvy sloupcu matice pro indexovani
    end
    
    methods (Access = public)
         function obj = CEpiEvents(d,tabs,fs)
             if isstruct(d)
                DE = d; %#ok<PROP>                
             else %pokud neni zadano fs, predpoklada se, ze prvni parametr jsou uz vyhodnocene eventy v DE
                DE = spike_detector_hilbert_v16_byISARG(d,fs);%#ok<PROP>             
             end
             %chci pridat jeste tabs do objektu, to se mi pak bude hodit pri zobrazovani v epochovanych datech
             obj.DE = [DE.pos, DE.dur, DE.chan, DE.con, DE.weight, DE.pdf, tabs(round(DE.pos*fs))]; %#ok<PROP>
             obj.sloupce = {};
             obj.sloupce.pos = 1; %pozice v sekundach zaznamu
             obj.sloupce.dur = 2; %trvani episod
             obj.sloupce.chan = 3; %cislo kanalu
             obj.sloupce.con = 4; %type (1-obvious 0.5-ambiguous)
             obj.sloupce.weight = 5; %váha 0-1
             obj.sloupce.pdf = 6; %statistical significance "PDF"
             obj.sloupce.tabs = 7;
         end
         function [time,weight]= GetEvents(obj,tabs,ch,tabs0)
             % vraci casy a vahy udalosti v danem casovem intervalu (tabs(:,1) - tabs(:,2)) na danem kanalu
             % tabs muze byt i matice o vice radcich, pak to odpovida vice opocham a udalosti se
             % spojuji pod sebe
             %time = zeros(size(obj.DE,1)*size(tabs,1),1);
             %weight  = zeros(size(obj.DE,1)*size(tabs,1),1);  
             time = []; %prazdny array, abych neco vratil
             weight = []; 
             for r = 1:size(tabs,1)                 
                 iDE = obj.DE(:,obj.sloupce.tabs) >= tabs(r,1) & obj.DE(:,obj.sloupce.tabs) <= tabs(r,2) ...
                    & obj.DE(:,obj.sloupce.chan)==ch & obj.DE(:,obj.sloupce.con) == 1;
                 DEx = obj.DE(iDE,:); %nevim predem velikost, takze nemuzu pole vytvorit predem
                 if size(DEx,1)>0 
                     if size(tabs,1)>0
                        % predpokladam, ze jsou vsechny epochy stejne dlouhe 
                        % u epochovanych dat chci casy od 0 a navazujici na sebe v navazujicich epochach
                        % nejdriv odectu pocet vterin mezi zacatkem epochy a zacatkem celeho zaznamu odkud je indexovane pos
                        % a pak prictu celkovy cas predchozich epoch
                        pos = DEx(:,obj.sloupce.pos) -  (tabs(r,1)-tabs0)*24*3600 + (r-1)*(tabs(1,2)-tabs(1,1))*24*3600;                        
                     else
                        %u neepochovanych dat chci puvodni cas - sec od zacatky
                        pos =  DEx(:,obj.sloupce.pos); 
                     end
                     time = [time ; pos]; %#ok<AGROW> %cas udalosti  v sec
                     weight = [weight ; DEx(:,obj.sloupce.weight)]; %#ok<AGROW> %vaha udalosti                     
                 end
             end
         end
    end
    
end

