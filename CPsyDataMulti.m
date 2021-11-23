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
        function Responses2XLS(obj,xlslabel, Wp)
            %export xls table for responses of all patients in CHilbertMulti file
            %similar to psydataavg(), but this function works over all patients without CHilbertMulti file
            %Wp - use for trialtypes
            if ~exist('xlslabel','var') || isempty(xlslabel) , xlslabel = ''; end
            if ~exist('Wp','var'), Wp = []; end
            iS_backup =obj.iS; %backup the current active subject
            [katnum, katstr] = obj.Categories(0,Wp); %assume same categories in all subjects
            varnames = {'n','subjname'}; %columnames
            varn0 = numel(varnames); %number of beginning columns before the actual measures rt and resp        
            for ikat = 1:numel(katnum) %measures - rt + resp, mean + stderr, for all stimulus categories
                varnames = horzcat(varnames,{ ['rt_' katstr{ikat} '_mean'],['rt_' katstr{ikat} '_stderr'],['rt_' katstr{ikat} '_num']...
                        ['resp_' katstr{ikat} '_mean'],['resp_' katstr{ikat}  '_stderr'],['resp_' katstr{ikat}  '_num']}); %#ok<AGROW>
            end
            katn0 = (numel(varnames) - varn0) / numel(katnum); %number of columns per stimulus category - always 6 ?
            output = cell(obj.nS,numel(varnames)); %xls table data
            for iS=1:obj.nS %#ok<*PROPLC,*PROP> %over all subjects
                obj.SubjectChange(iS);                
                [resp,rt,kat,test] = obj.GetResponses(Wp);       
                %TODO obj.GetTrialType
                output(iS,1:varn0) = {num2str(iS),obj.P.pacientid};
                for ikat=1:numel(katnum) %over all categories
                    iresp = kat==katnum(ikat) & test==1; %index for this category for all test responses
                    irt = iresp & resp==1; %rt compute only from correct responses                    
                    means = [mean(rt(irt)) mean(resp(iresp))]; %mean rt and resp for this subject and this category
                    stderr = [std(rt(irt))/sqrt(length(rt(irt)))  std(resp(iresp))/sqrt(length(resp(iresp)))]; %stderr of rt and resp
                    num = [sum(irt) sum(iresp) ]; %number of values for rt / resp
                    output(iS,varn0+(ikat-1)*katn0+1 : varn0+(ikat-1)*katn0+katn0) = {double2str(means(1),3) , double2str(stderr(1),3),num2str(num(1)) , double2str(means(2),3) , double2str(stderr(1),3),num2str(num(2)) };                    
                end                
            end
            if isfield(Wp,'trialtypes') && ~isempty(Wp.trialtypes) && iscell(Wp.trialtypes)
                ttname = strrep(cell2str(Wp(1).trialtypes,1),' ',''); %short strig of the trialtypes               
            else
                ttname = '';
            end
            katname = ['kats' strrep(cell2str(Wp(1).kats,1),'  ','-')]; %short string of the stimulus kategories
            xlsfilename = ['./logs/Responses2XLS_' xlslabel '_' katname '_' ttname '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')];
            xlswrite(xlsfilename ,vertcat(varnames,output)); %write to xls file
            disp([xlsfilename '.xls, with ' num2str(size(output,1)) ' lines saved']);
            obj.SubjectChange(iS_backup); %activate the orignal subject
        end
        function epochs = GetEpochsMax(obj)
            %returns the maximum number of epochs over all subjects
            epochs = 0;
            for s=1:obj.nS
                epochs=max(epochs,size(obj.Pmulti(s).data,1));
            end
        end
    end
    
end

