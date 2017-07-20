classdef CEpiEvents < handle
    %CEPIEVENTS Trida na praci s epileptickymi eventy detekovanymi detektorem ISARG
    %Kamil Vlcek, FGU AVCR, since 2016 09   
    
    properties (Access = public)
        DE; %puvodni struct vygenerovany detektorem
        sloupce; %nazvy sloupcu matice pro indexovani
        d; %matice eventu pro snadne vyhledavani
        iDEtabs;
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
         
         function [time,weight,chans]= GetEvents(obj,tabs,ch,tabs0)
             % vraci casy a vahy udalosti v danem casovem intervalu (tabs(:,1) - tabs(:,2)) na danem kanalu
             % tabs muze byt i matice o vice radcich, pak to odpovida vice opocham a udalosti se
             % spojuji pod sebe
             % ch - cislo kanalu
             % tabs0 - timestamp prvni udalosti celeho zaznamu
             %time = zeros(size(obj.DE,1)*size(tabs,1),1);
             %weight  = zeros(size(obj.DE,1)*size(tabs,1),1);  
             time = []; %prazdny array, abych neco vratil
             weight = []; 
             chans = [];
             if numel(obj.iDEtabs)==0 %pokud zadny index pres tabs neexistuje, vytvorim ho
                 obj.iDEtabs = false([size(obj.d,1) size(tabs,1)]);
                 for r = 1:size(tabs) %pres zjistovane epochy
                     %tohle v profileru trvalo hrozne dlouho protoze se to delal pro kazdy kanala a kazdy event (kdy jsou). 
                     %ted se to dela jen pro kazdy event
                     obj.iDEtabs(:,r) = obj.d(:,obj.sloupce.tabs) >= tabs(r,1) & obj.d(:,obj.sloupce.tabs) <= tabs(r,2) & obj.d(:,obj.sloupce.con) == 1;
                     %logicke 1 pokud je udalost casove v ramci epochy a je 1=obvious
                 end
             end
             for r = 1:size(tabs,1)                 
                 iDE = obj.iDEtabs(:,r) & ismember(obj.d(:,obj.sloupce.chan),ch); %pouziju index iDEtabs vytvoreny drive - %10.7.2017 - muze byt vic kanalu najednou
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
                     time = [time ; pos repmat(r,size(pos))]; %#ok<AGROW> %cas udalosti  v sec
                     weight = [weight ; DEx(:,obj.sloupce.weight)]; %#ok<AGROW> %vaha udalosti                     
                     chans = [chans; DEx(:,obj.sloupce.chan)]; %#ok<AGROW> % cisla kanalu
                 end
             end
         end
         function obj = Clear_iDEtabs(obj)
             %vymaza index obj.iDEtabs vytvoreny ve funkci GetEvents - pouziva se pri jejichopakovanych volanich 
             obj.iDEtabs = [];
         end
         function [evts,names,epochs,evts_nonseeg]= CountEpiEvents(obj,CHH,objepochs,objtabs,objtabs_orig)
            %spocita epilepticke event pro jednotlive kanaly, vcetne poctu epoch
            evts = zeros(numel(CHH.H.selCh_H),1); %pocet epievents pro kazdy kanal
            names = cell(numel(CHH.H.selCh_H),1);  %jmena kanalu
            epochs = evts; 
            namelast = ''; %tam budu ukladat pismeno jmena elektrody
            evts_nonseeg = 0; %pocet epi udalost v nonseeg kanalech, napr. trigger
            obj.Clear_iDEtabs();
            if objepochs > 1 %epochovana data
                tabs = [ objtabs(1,:)' objtabs(end,:)' ]; % %zacatky a konce vsech epoch
            end 
            for ch = 1:CHH.ChannelsN() %obj.CH.H.selCh_H                
                if objepochs==1
                    [epitime, ~] = obj.GetEvents([objtabs(1) objtabs(end)],ch);
                else %epochovana data                   
                    [epitime, ~] = obj.GetEvents(tabs,ch,objtabs_orig(1));                                          
                end
                if ~isempty(find(CHH.H.selCh_H==ch, 1)) %pokud je kanal ve vyjmenovanych SEEG kanalech podle headeru
                    evts(ch) = size(epitime,1); 
                    if numel(epitime) > 0 && objepochs > 1
                        epochs(ch) = numel( unique(epitime(:,2)));  %v druhem sloupci jsou cisla epoch                       
                    end
                    if strcmp(namelast,CHH.H.channels(ch).name(1))
                       names{ch} = CHH.H.channels(ch).name(end); %cislo elektrody bez pismene, u druhe elektrody stejneho jmena abych usetril misto
                    else
                       names{ch} = CHH.H.channels(ch).name(1);
                       namelast = CHH.H.channels(ch).name(1);                     
                    end
                else
                   evts_nonseeg = evts_nonseeg + size(epitime,1);
                end
                
            end
         end
         
         function [RjEpoch,RjEpochCh,vyrazeno] = RejectEpochsEpi(obj,NEpi,CHH,objepochs,objtabs,objtabs_orig)
            %vyradi epochy podle poctu epileptickych udalosti - pokud >= NEpi
            %uz vyrazene epochy rucne nemeni - neoznaci jako spravne
            assert( objepochs > 1, 'nejsou epochovana data');
            vyrazeno=0;
            RjEpoch = [];
            if isempty(NEpi)
                NEpi = 0;
                RjEpochCh = false(CHH.ChannelsN(),objepochs); 
            else
                RjEpochCh = []; %musim neco vratit, i kdyz pouzivam NEpi
            end
            for e = 1:objepochs
                obj.Clear_iDEtabs(); %protoze to pole je zavisle na epochach, musim mazat pokazde
                tabs = [ objtabs(1,e)' objtabs(end,e)' ]; % %zacatky a konce zobrazenych epoch
                [epitime, ~, chans] = obj.GetEvents( tabs,CHH.GetChOK(),objtabs_orig(1));                 
                if NEpi
                    if numel(epitime) >= NEpi
                        RjEpoch = unique( [RjEpoch e]); %pridam epochu do seznamu, pokud tam jeste neni
                        vyrazeno = vyrazeno +1;
                    end
                elseif numel(epitime) > 0
                    RjEpochCh(unique(chans),e) = true;
                    vyrazeno = vyrazeno + 1;
                end
            end                            
            
            
        end
    end
    
end

