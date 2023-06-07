classdef CPsyDataMulti < CPsyData
    %CPSYDATAMULTI shromazduje data z CPsyData z vice subjektu dohromady
    %   kvuli CHilbertMulti
    
    properties (Access = public)
        nS; %pocet nactenych subjektu
        iS; %aktualni index subjektu
        Pmulti; %shromazduje P data od subjektu 
        blocksMulti; %zalohuju spocitane bloky, protoze trva hodne dlouho je spocitat
        fhRxls; %figure handle from Responses2XLS
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
        function [obj]= FilterEpochs(obj,iepochs)
            assert(isempty(obj.iepochs),'epochs are already filtered');             
            for iS=1:obj.nS %#ok<*PROPLC,*PROP> %over all subjects
                %obj.SubjectChange does not write to obj.Pmulti, so it cannot be used for this purpose
                obj.Pmulti(iS).data = obj.Pmulti(iS).data(iepochs,:); %filter directly
                %leave only the epochs matching the filter 
            end
            if ~isempty(obj.trialtypes)
                obj.trialtypes = obj.trialtypes(iepochs,:);%leave only the trialtypes for epochs matching the filter 
            end
            obj.iepochs = iepochs;
        end
        function Responses2XLS(obj,Wp,makefile,xlslabel)
            %export xls table for responses of all patients in CHilbertMulti file
            %similar to psydataavg(), but this function works over all patients without CHilbertMulti file
            %Wp - used to select specific stats, as Wp from the main object is not accessible
            if ~exist('xlslabel','var') || isempty(xlslabel) , xlslabel = ''; end
            if ~exist('Wp','var'), Wp = []; end
            assert(isempty(Wp) || isstruct(Wp) && numel(Wp) == 1, 'Wp needs to be struct 1x1');
            if ~exist('makefile','var'), makefile = 1; end
            iS_backup =obj.iS; %backup the current active subject
            [katnum, katstr] = obj.Categories(0,Wp); %assume same categories in all subjects
            varnames = {'n','subjname'}; %columnames
            titleline = {'',''};
            varn0 = numel(varnames); %number of beginning columns before the actual measures rt and resp        
            for ikat = 1:numel(katnum) %measures - rt + resp, mean + stderr, for all stimulus categories
                varnames = horzcat(varnames,{ ['rt_' katstr{ikat} '_mean'],['rt_' katstr{ikat} '_stderr'],['rt_' katstr{ikat} '_num']...
                        ['resp_' katstr{ikat} '_mean'],['resp_' katstr{ikat}  '_stderr'],['resp_' katstr{ikat}  '_num']}); %#ok<AGROW>
                titleline = horzcat(titleline,{katstr{ikat}},repmat({''},1,5)); %#ok<AGROW>
            end
            katn0 = (numel(varnames) - varn0) / numel(katnum); %number of columns per stimulus category - always 6 ?
            varnames = horzcat(varnames,{'rt_p','resp_p'}); %ttest resuts as last columns in the xls table
            titleline =  horzcat(titleline,{'STAT',''});
            output = cell(obj.nS,numel(varnames)); %xls table data
            outS = struct(); %all values from all subjects for the first-level statistics
            prumery = zeros(obj.nS,numel(katnum),3,2); %subjects x kats  x [mean,strerr,p] x [rt,resp]- averages and std errs a ttest pvalues for all subjects, to be plotted  
            
            %cycle over subject and stimulus categories / trialtypes
            for iS=1:obj.nS %#ok<*PROPLC,*PROP> %over all subjects
                obj.SubjectChange(iS);                
                [resp,rt,kat,test] = obj.GetResponses(Wp);                       
                %TODO obj.GetTrialType
                output(iS,1:varn0) = {num2str(iS),obj.P.pacientid};
                for ikat=1:numel(katnum) %over all categories
                    iresp = any(kat==cellval(katnum,ikat),2) & test==1; %index for this category for all test responses
                    irt = iresp & resp==1; %rt compute only from correct responses                    
                    means = [mean(rt(irt)) mean(resp(iresp))]; %mean rt and resp for this subject and this category                    
                    stderr = [std(rt(irt))/sqrt(length(rt(irt)))  std(resp(iresp))/sqrt(length(resp(iresp)))]; %stderr of rt and resp
                    prumery(iS,ikat,:,:) = [means; stderr; 1 1];
                    num = [sum(irt) sum(iresp) ]; %number of values for rt / resp
                    output(iS,varn0+(ikat-1)*katn0+1 : varn0+(ikat-1)*katn0+katn0) = {double2str(means(1),3) , double2str(stderr(1),3),num2str(num(1)) , double2str(means(2),3) , double2str(stderr(1),3),num2str(num(2)) };                    
                    %save values for statistics
                    outS(iS,ikat).pacientid = obj.P.pacientid;
                    outS(iS,ikat).resp = resp(iresp); 
                    outS(iS,ikat).rt = rt(irt);
                    outS(iS,ikat).kat = katnum(ikat);                    
                end 
                if numel(katnum) == 2 %ttest if there are two stimulus categories
                    [~,rtP] = ttest2(outS(iS,1).rt,outS(iS,2).rt);                  
                    [~,respP] = ttest2(outS(iS,1).resp,outS(iS,2).resp);                   
                    if isnan(respP) %when both arrays contain only 1, the ttest2 returns pvalue = nan
                        respP = 1; %in this case we consider the arrays not to be different
                    end
                else
                    respP = 1; rtP = 1; %for one or more than two values, no stat is computed
                end
                output(iS,end-1:end) = {double2str(rtP,3),double2str(respP,3)};
                prumery(iS,ikat,3,:) = [rtP,respP];
            end  
            
            % export XLS file
            if isfield(Wp,'trialtypes') && ~isempty(Wp.trialtypes) && iscell(Wp.trialtypes)
                ttname = strrep(cell2str(Wp(1).trialtypes,1),' ',''); %short string of the trialtypes               
            else
                ttname = '';
            end
            katname = ['kats' strrep(cell2str(Wp(1).kats,1),'  ','-')]; %short string of the stimulus kategories
            titleline(1:2) = {katname,ttname};
            if makefile
                xlsfilename = ['./logs/Responses2XLS_' xlslabel '_' katname '_' ttname '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')];
                xlswrite(xlsfilename ,vertcat(titleline,varnames,output)); %write to xls file
                disp([xlsfilename '.xls, with ' num2str(size(output,1)) ' lines saved']);
            end
            
            %plot figure of all subjects - private function
            obj.Responses2Plot(prumery,Wp,katname,ttname,katnum,katstr);
            
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
    methods (Access = private)
        function Responses2Plot(obj,prumery,Wp,katname,ttname,katnum,katstr)
            %plot figure of behavioral data of all subjects
            % called from Responses2XLS, which collects the data to plot
            % prumery:subjects x kats  x [mean,strerr,p] x [rt,resp]- averages and std errs a ttest pvalues for all subjects, to be plotted 
            % ttname:short string of the trialtypes
            if isprop(obj,'fhRxls') && ~isempty(obj.fhRxls) && ishandle(obj.fhRxls) 
                figure(obj.fhRxls); %pokud uz graf existuje, nebudu tvorit znova
                clf; %smazu aktualni figure
            else                
                obj.fhRxls = figure('Name','PsyData.Responses2XLS');
            end            
            
            subplot(2,1,1); %reaction times
            hold all;
            for ikat=1:numel(katnum)                
                errorbar(1:obj.nS,prumery(:,ikat,1,1),prumery(:,ikat,2,1),'.'); %,'color',colorskat{2,k});                 
            end            
            signif = find(prumery(:,ikat,3,1)<=0.05); % subjects with signif difference between two categories
            if ~isempty(signif)
                plot(signif,0.5,'*','color',[1 0 0]); %mark subjects with significant differences
            end
            xticks(1:obj.nS);
            xticklabels({obj.Pmulti.pacientid});
            set(gca,'TickLabelInterpreter','none'); % plot _ as _ in the subject names 
            xtickangle(90);
            legend(katstr);
            if numel(Wp.kats)==1                
                titlename = [katname ' ' obj.CategoryName(Wp.kats) ', ' ttname  ];
            elseif ~isempty(Wp.trialtypes) && numel(Wp.trialtypes)<=2                
                titlename = [katname  ', ' ttname ' ' obj.TrialTypeName(Wp.trialtypes) ];
            else
                titlename = [katname ', ' ttname  ];
            end
            title(titlename, 'Interpreter', 'none');
            ylabel('reaction time');
            
            subplot(2,1,2); %response accuracy
            hold all;
            for ikat=1:numel(katnum)                
                errorbar(1:obj.nS,prumery(:,ikat,1,2),prumery(:,ikat,2,2),'.'); %,'color',colorskat{2,k}); 
            end
            signif = find(prumery(:,ikat,3,2)<=0.05);
            if ~isempty(signif)
                 plot(signif,0.5,'*','color',[1 0 0]); %mark subjects with significant differences
            end
            xticks(1:obj.nS);
            xticklabels({obj.Pmulti.pacientid});
            set(gca,'TickLabelInterpreter','none'); % plot _ as _ in the subject names
            xtickangle(90);
            ylabel('response accuracy');
             
        end
    end
    
end

