classdef CPlots < handle
    %CPLOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static,Access = public)
        function PlotElectrodeEpiEvents(els,RjCh,DE,objtabs,tabs_orig,epochs,samples,epochtime,timeax,timeaxy,sec,time_n,elmaxmax,shift,iD)
            %27.4. - nahrazuju cast funkce CiEEGData.PlotElectrode toutu funkci, kvuli zkraceni, ale predavam strasne moc argumentu
            %treba bych mohl uz driv predat hodne z nich? A mit je jako properties tehle tridy?
            epieventsum = 0;  
            weights = [];
            DE.Clear_iDEtabs();
            for ch = els
                if ~ismember(ch,RjCh)
                    if epochs <= 1 %neepochovana data
                        [epitime, weight] = DE.GetEvents( [objtabs(iD(1)) objtabs(iD(2))],ch,tabs_orig(1)); 
                    else
                        epochy = sec : min(sec+ceil(time_n/samples)-1 , epochs ); %cisla zobrazenych epoch
                        tabs = [ objtabs(1,epochy)' objtabs(end,epochy)' ]; %zacatky a konce zobrazenych epoch
                        [epitime, weight] = DE.GetEvents( tabs,ch,tabs_orig(1)) ;                            
                        if numel(epitime) > 0
                            epitime(:,1) = epitime(:,1) + epochtime(1);
                        end
                    end
                    if size(epitime,1) > 0
                        plot(epitime(:,1),shift(elmaxmax-ch+1,1),'o', 'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',6);
                        epieventsum = epieventsum + size(epitime,1);
                        weights = [weights; round(weight*100)/100]; %#ok<AGROW>
                    end
                end
            end             
            text(timeax(1)+((timeax(end)-timeax(1))/8),timeaxy,['epileptic events: ' num2str(epieventsum) ' (' num2str(min(weights)) '-' num2str(max(weights)) ')']);         
        end
        function EpochsEpi(RjEpochCh,objels,CHH)
            figure('Name','Rejected epochs in individual channels');
            pocty = sum(RjEpochCh,2)/size(RjEpochCh,2); %pocty epoch v jednotlivych kanalech / celkovym poctem epoch
            plot(pocty,'.-');
            for el = 1:numel(objels)-1
                line([objels(el) objels(el)]+1,[0 1],'Color',[0.5 0.5 0.5]);
            end
            line([1 CHH.H.selCh_H(end)],[0.30 0.30],'Color','red');
            spatne = find(pocty >= 0.30); 
            disp(['kanaly s pocet vyrazenych epoch >= 0.3:' mat2str(spatne)]);                
        end
    end
    
end

