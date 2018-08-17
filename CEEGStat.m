classdef CEEGStat
    %CEEGSTAT trida na zpracovani EEG dat, oddelene od vlastni tridy CiEEGData
    %   sem budu postupne prevadet kod z CiEEGData
    % since 18.7.2017
    
    properties (Access = private)
       d; %matice EEG dat 
       fs; %vzorkovaci frekvence
    end
    
    methods (Access = public)
        
        function obj = CEEGStat(d,fs)
            %konstruktor
            obj.d = d;            
            obj.fs = fs;
        end
        
        function [P,ibaseline,iepochtime,itimewindow] = WilcoxBaseline(obj,epochtime,baseline,timewindow,iEp,RjEpCh)
            % spocita signifikance vuci baseline
            % prevedeno z CiEEGData.ResponseSearch 18.7.2017
            iepochtime = round(epochtime(1:2).*obj.fs); %v poctu vzorku cas pred a po udalosti, pred je zaporne cislo           
            ibaseline = round(baseline.*obj.fs); %zaporna cisla pokud pred synchro eventem
            ibaselineA =  [ibaseline(1) ,    floor((ibaseline(2)-ibaseline(1))/2)+ibaseline(1)]; %prvni pulka baseline - 28.3.2017
            ibaselineB =  [ibaselineA(2)+1 , ibaseline(2) ]; %druha pulka baseline
            itimewindow = round(timewindow.*obj.fs); %
            
            %proc driv?  mean(obj.d( abs(iepochtime(1)-ibaselineA(1))+1 : abs(iepochtime(1)-ibaselineA(2))  - 28.3.2017            
            baselineA = mean(obj.d( ibaselineA(1)-iepochtime(1)+1 : ibaselineA(2)-iepochtime(1) , : , iEp),1); 
            baselineB = mean(obj.d( ibaselineB(1)-iepochtime(1) : ibaselineB(2)-iepochtime(1) , : , iEp),1); 
                % cas x kanaly x epochy - prumer za cas pred podnetem, pro vsechny kanaly a vsechny nevyrazene epochy
            if numel(itimewindow) == 2  %chci prumernou hodnotu eeg odpovedi d od do
                response = mean(obj.d( (itimewindow(1) : itimewindow(2)) - iepochtime(1) , : , iEp ),1); %prumer v case                
            else %time windows je delka maximalni hodnoty p - takze tady jeste neresim
                response = obj.d(abs(iepochtime(1)-ibaseline(2))+1 : end , : , iEp ); %vsechny hodnoty po konci baseline                 
            end            
            WpA = CStat.Wilcox2D(response,baselineA,1,[],'mean vs baseline A',RjEpCh(:,iEp),RjEpCh(:,iEp));  %1=mene striktni pdep, 2=striktnejsi dep;
            WpB = CStat.Wilcox2D(response,baselineB,1,[],'mean vs baseline B',RjEpCh(:,iEp),RjEpCh(:,iEp));  %1=mene striktni pdep, 2=striktnejsi dep;
                %z vyrazenych epoch pro kazdy kanal (RjEpCh) vezmu jen ty celkove nevyrazene
            Wp = max(WpA,WpB); % %vyssi hodnota z kazde poloviny baseline
            if numel(itimewindow) == 1 %chci maximalni hodnotu p z casoveho okna
                P = CStat.Klouzaveokno(Wp,itimewindow(1),'max',1); % %pole 2D signifikanci si ulozim kvuli kresleni - cas x channels                
            else
                P = Wp; %%pole 1D signifikanci - jedna hodnota pro kazdy kanal                
            end
        end
    end
    methods (Static,Access = public)
        function [WpKat,WpKatBaseline]=WilcoxCat(kats,responsekat,baselinekat,rjepchkat,itimewindow,method,Pbaseline)
            %spocita statistiky pro ruzne kategorie - vuci sobe navzajem a vuci baseline  
            %Pbaseline urcuje jestli se ma pocitat stat pro kazdy kanal zvlast (Pbaseline=neco) nebo dohromady (Pbaseline=[])
            
            %rozdily kategorii vuci baseline - 28.3.2017
            WpKatBaseline = cell(numel(kats),1);  %rozdily kategorii vuci baseline   
            for k =  1: numel(kats)                
                 baselineall = baselinekat{k};
                 baselineA = mean(baselineall(1:floor(size(baselineall,1)/2)      ,:,:));
                 baselineB = mean(baselineall(  floor(size(baselineall,1)/2)+1:end,:,:));
                 if exist('Pbaseline','var') && ~isempty(Pbaseline) %nova verze, statistika pro kazdy kanal zvlast, tam kde je odpoved vuci baseline
                     WpKatBaseline{k,1} = ones(size(responsekat{k},1), size(responsekat{k},2)); %time x channels
                     fprintf('kat %i vs baseline Channels Wilcox2D: 1 ... ',k);
                     for ch = 1:size(responsekat{k},2) %musim baseline signif taky pocitat pro kazdy kanal zvast protoze ji pak zohlednuju mezi kategoriemi
                        WpBA = CStat.Wilcox2D(responsekat{k}(:,ch,:),baselineA(:,ch,:),0,[],['kat ' num2str(k) ' vs baseline A'],rjepchkat{k}(ch,:),rjepchkat{k}(ch,:));
                        WpBB = CStat.Wilcox2D(responsekat{k}(:,ch,:),baselineB(:,ch,:),0,[],['kat ' num2str(k) ' vs baseline B'],rjepchkat{k}(ch,:),rjepchkat{k}(ch,:));
                        WpKatBaseline{k,1}(:,ch) = max (WpBA,WpBB);
                     end
                     fprintf('... %i \n',ch);
                 else %puvodni verze, statistika pro vsechny kanaly najednou  
                     WpBA = CStat.Wilcox2D(responsekat{k},baselineA,1,[],['kat ' num2str(k) ' vs baseline A'],rjepchkat{k},rjepchkat{k});
                     WpBB = CStat.Wilcox2D(responsekat{k},baselineB,1,[],['kat ' num2str(k) ' vs baseline B'],rjepchkat{k},rjepchkat{k});
                     WpKatBaseline{k,1} = max (WpBA,WpBB);
                 end
            end
            
            %rozdily kategorii vuci sobe
            WpKat = cell(numel(kats)); %rozdily mezi kategorieme
            for k = 1:numel(kats) %budu statisticky porovnavat kazdou kat s kazdou, bez ohledu na poradi
                for j = k+1:numel(kats)
                    if strcmp(method,'wilcox')
                        if exist('Pbaseline','var') && ~isempty(Pbaseline) %nova verze, statistika pro kazdy kanal zvlast, tam kde je odpoved vuci baseline
                            Wr = ones(size(WpKatBaseline{1})); %vysoke pravdepodobnosti
                            fprintf('kat %i vs %i Channels Wilcox2D: 1 ... ',k,j);
                            for ch = 1:size(Wr,2)                                
                                if min(WpKatBaseline{k}(:,ch))<.05 || min(WpKatBaseline{j}(:,ch))<.05
                                    Wr(:,ch) = CStat.Wilcox2D(responsekat{k}(:,ch,:), responsekat{j}(:,ch,:),0,[],['kat ' num2str(k) ' vs ' num2str(j)],rjepchkat{k}(ch,:),rjepchkat{j}(ch,:)); % -------- WILCOX kazda kat s kazdou                                     
                                    %fprintf('%i,',ch);
                                end                                
                            end
                            fprintf('... %i \n',ch);
                        else  %puvodni verze, statistika pro vsechny kanaly najednou  
                            Wr = CStat.Wilcox2D(responsekat{k}, responsekat{j},1,[],['kat ' num2str(k) ' vs ' num2str(j)],rjepchkat{k},rjepchkat{j}); % -------- WILCOX kazda kat s kazdou 
                        end
                    elseif strcmp(method,'permut')
                        Wr = CStat.PermStat(responsekat{k}, responsekat{j},1,['kat ' num2str(k) ' vs ' num2str(j)],rjepchkat{k},rjepchkat{j}); % -------- Permutacni test kazda kat s kazdou 
                    else
                        disp('neznama metoda statistiky');
                        Wr = [];
                    end
                    if numel(itimewindow)==0 
                        WpKat{k,j} = Wr;  % pokud nechci pouzit klouzave okno
                    else
                        WpKat{k,j} = CStat.Klouzaveokno(Wr,itimewindow(1),'max',1); %pri nove signifikanci jen pro jeden kanal chci zase pouzivat klouzave okno
                    end
                end
            end
            
        end
        function baseline0 = Baseline(epochtime,baseline)
            %baseline0 je cast epochtime pred koncem baseline nebo pred casem 0
            baseline0 = [epochtime(1) iff(baseline(2)>epochtime(1),baseline(2),0)]; %zalezi jestli se baselina a epochtime prekryvaj
        end
    end
    
end

