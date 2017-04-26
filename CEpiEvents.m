classdef CEpiEvents < handle
    %CEPIEVENTS Trida na praci s epileptickymi eventy detekovanymi detektorem ISARG
    %Kamil Vlcek, FGU AVCR, since 2016 09   
    
    properties (Access = public)
        DE; %puvodni struct vygenerovany detektorem
        sloupce; %nazvy sloupcu matice pro indexovani
        d; %matice eventu pro snadne vyhledavani
    end
    
    methods (Access = public)
         function obj = CEpiEvents(d,tabs,fs,mults)
             if exist('mults','var') && ~isempty(mults)
                    assert(size(mults,1)<= 1 || size(mults,1)==size(d,2),'d and mults have to have same number of channels');
                    d = bsxfun(@times,double(d), mults); %rovnou to roznasobim mults, nechci to resit dodatecne - 14.2.2017            
             end
             if isstruct(d)
                DE = d; %#ok<PROP>                
             else %pokud neni zadano fs, predpoklada se, ze prvni parametr jsou uz vyhodnocene eventy v DE
                DE = spike_detector_hilbert_v16_byISARG(d,fs);%#ok<PROP>             
                disp(['nalezeno ' num2str(numel(DE.pos)) ' epileptickych udalosti']); %#ok<PROP>
             end
             %chci pridat jeste tabs do objektu, to se mi pak bude hodit pri zobrazovani v epochovanych datech
             obj.DE=DE; %#ok<PROP> %ukladam jen do zalohy, abych mohl rucne ulozit na disk
             obj.d = [DE.pos, DE.dur, DE.chan, DE.con, DE.weight, DE.pdf, tabs(round(DE.pos*fs))]; %#ok<PROP> %to budu dal pouzivat
             obj.sloupce = {};
             obj.sloupce.pos = 1; %pozice v sekundach zaznamu
             obj.sloupce.dur = 2; %trvani episod
             obj.sloupce.chan = 3; %cislo kanalu
             obj.sloupce.con = 4; %type (1-obvious 0.5-ambiguous)
             obj.sloupce.weight = 5; %váha 0-1
             obj.sloupce.pdf = 6; %statistical significance "PDF"
             obj.sloupce.tabs = 7;             
         end
         function [DE] = GetDE(obj) %jen vrati originalni struct DE vypocitany detektorem; protoze ho nejak nejde ulozit
             DE = obj.DE;
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
                 iDE = obj.d(:,obj.sloupce.tabs) >= tabs(r,1) & obj.d(:,obj.sloupce.tabs) <= tabs(r,2) ...
                    & obj.d(:,obj.sloupce.chan)==ch & obj.d(:,obj.sloupce.con) == 1;
                 DEx = obj.d(iDE,:); %nevim predem velikost, takze nemuzu pole vytvorit predem
                 if size(DEx,1)>0 
                     if size(tabs,1)>1 %epochovana data
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

