classdef CPsyDataMulti < CPsyData
    %CPSYDATAMULTI shromazduje data z CPsyData z vice subjektu dohromady
    %   kvuli CHilbertMulti
    
    properties (Access = public)
        nS; %pocet nactenych subjektu
        iS; %aktualni index subjektu
        Pmulti; %shromazduje P data od subjektu 
        iStat; %index statistiky
        nStat; %pocet nactenych statistiky
    end
    
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS
        function obj = CPsyDataMulti(psy)
            %konstruktor            
            obj@CPsyData(psy);            
            obj.iS = 1;
            obj.Pmulti = psy; %nyni dvourozmerne Subject x Statistika
            obj.nS = 1;
            obj.iStat = 1;
            obj.nStat = 1;
        end
        function [obj] = SubjectChange(obj,iS)
            %aktivuje data z udaneho subjektu, pro stejnou statistiku
            if iS ~= obj.iS && iS <= obj.nS
                obj.iS = iS;
                obj.P = obj.Pmulti(iS, obj.iStat);                
                %disp(['subject changed to ' num2str(iS)]);
                obj.blocks = []; %to plati jen pro aktualni subjekt
            end
        end
        function [obj] = StatChange(obj,iStat)
            %aktivuje data z nove statistiky, pro stejny subjekt
            if iStat <= obj.nStat
                obj.iStat = iStat;
                obj.P = obj.Pmulti(obj.iS, iStat);
                %disp(['stat changed to ' num2str(iStat)]);
            end
        end
        function [obj]= GetPsyData(obj,psy)
            %nacte data z dalsiho subjektu a rovnou ho aktivuje; ulozi do aktualni statistiky (druhy rozmer Pmulti)
            testname = obj.GetTestName(inputname(2));
            assert(strcmp(testname,obj.testname),['data ze dvou ruznych testu: ' testname ' x ' obj.testname]);
            obj.P = psy;
            obj.DoplnZpetnavazba(); %pokud neni v puvodnich datech, doplnim sloupec zpetnavazba
            obj.Pmulti(obj.iS+1,obj.iStat) = obj.P;
            obj.nS = obj.nS + 1;
            obj.SubjectChange(obj.iS + 1, obj.iStat); %inkrementuju subjekt a aktivuju jeho data
        end   
        
        function [obj]= GetStatData(obj,psy)
            %nacte data z dalsi statistiky, pro aktualni subjekt, a rovnou ji aktivuje
            obj.Pmulti(obj.iS, obj.iStat + 1) = psy;
            obj.nStat = obj.nStat + 1;
            obj.StatChange(obj.iStat + 1);
        end
        
        function [obj]= UpdateStatData(obj,psy)
            %prepise psy data pro aktualni statistiku a aktualni subjekt
            obj.Pmulti(obj.iS, obj.iStat) = psy;
            % mozno este sem pridat statchange?
        end
    end
    
end

