classdef CPsyDataMulti < CPsyData
    %CPSYDATAMULTI shromazduje data z CPsyData z vice subjektu dohromady
    %   kvuli CHilbertMulti
    
    properties (Access = public)
        nS; %pocet nactenych subjektu
        iS; %aktualni index subjektu
        Pmulti; %shromazduje P data od subjektu 
        blocksMulti; %zalohuju spocitane bloky, protoze trva hodne dlouho je spocitat
    end
    
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS
        function obj = CPsyDataMulti(psy)
            %konstruktor            
            obj@CPsyData(psy);            
            obj.iS = 1;
            obj.Pmulti = psy;
            obj.nS = 1;
            obj.blocksMulti = {[]}; %jedna prazdna bunka
        end
        function [obj] = SubjectChange(obj,iS)
            %aktivuje data z udaneho subjektu
            if iS ~= obj.iS && iS <= obj.nS
                if ~isempty(obj.blocks) 
                    obj.blocksMulti{obj.iS} = obj.blocks; %ulozim bloky,ktere byly zatim spocitany
                end
                obj.iS = iS;
                obj.P = obj.Pmulti(iS);                
                %disp(['subject changed to ' num2str(iS)]);
                if obj.iS <= numel(obj.blocksMulti)
                    obj.blocks = obj.blocksMulti{obj.iS}; %nactu spocitane bloky noveho subjektu
                else
                    obj.blocks = []; %nebo necham prazdne
                end
            end
        end
        function [obj]= GetPsyData(obj,psy)
            %nacte data z dalsiho subjektu a rovnou ho aktivuje
            testname = obj.GetTestName(inputname(2));
            assert(strcmp(testname,obj.testname),['data ze dvou ruznych testu: ' testname ' x ' obj.testname]);
            obj.P = psy;
            obj.DoplnZpetnavazba(); %pokud neni v puvodnich datech, doplnim sloupec zpetnavazba
            obj.Pmulti(obj.iS+1) = obj.P;
            obj.nS = obj.nS + 1;
            obj.SubjectChange(obj.iS + 1);
            obj.blocksMulti = [obj.blocksMulti, {[]}]; %pridam dalsi prazdny cell na konec
        end        
    end
    
end

