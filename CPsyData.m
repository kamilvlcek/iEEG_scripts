classdef CPsyData < matlab.mixin.Copyable %je mozne kopirovat pomoci E.copy();
    %CPSYDATA Trida na praci s behavioralnimi daty z psychopy
    %   vytvorena pomoci ppa_data.m, aedist_data.m aj
    % Kamil Vlcek, FGU AVCR, since 2016 04
    
    properties  (Access = public)
        P; %psychopy behavioural data
        fhR; %figure handle from PlotResponses
        warning_rt=false; %jestli uz byl warning o reakcnich casech
        testname; %jmeno testu, ze ktereho jsou data
        blocks; %struktura, ktera uklada vysledky GetBlocks, kvuli uspore casu
        trialtypes; %table with trialtypes specific for a tset
    end
    
    methods (Access = public)
        %% ELEMENTAL FUNCTIONS
        function obj = CPsyData(psy)
            %konstruktor
            obj.P = psy;
            if ~isfield(psy,'pacientid'), obj.P.pacientid = ''; end
            if ~isfield(psy,'eegfile'), obj.P.eegfile = ''; end
            obj.DoplnZpetnavazba();
            obj.testname = obj.GetTestName(inputname(1));
        end
        function [obj] = SubjectChange(obj,~)
            %jen prazdna funkce, ale pretezuju ji v CPsyDataMulti, a naplnuju obsahem
            % v tehle trida mam jen jeden subjekt
        end
        function [testname2] = GetTestName(obj,testname)  
            %zjisti jmeno testu, ze ktereho jsou data            
            if strcmp(testname,'aedist') || strcmp(testname,'menrot') || strcmp(testname,'ppa')
                testname2 = testname;
            elseif strcmp(obj.P.strings.podminka{1,1},'cervena')
                testname2 = 'aedist';
            elseif strcmp(obj.P.strings.podminka{1,1},'vy-2D')
                testname2 = 'menrot';
            elseif strcmp(obj.P.strings.podminka{1,1},'Ovoce')
                testname2 = 'ppa';
            else
                testname2 = '';
            end
        end
        function rt = ReactionTime(obj,kategories,trialtypes)
            %return the matrix of all reaction times split to columnts by stimulus category (if kategories is empty, ie. default)
            %using NaN where there are less values in a columns, it is therefore necessary to use nanmean or similar functions
            %the RT are from the time of sunchronization pulse (so not from the PsychoPy). 
            %if kategories = -1, only one column is returned without considering the category
            %if kategories is array with numbers of categories, return only these
            %trialtypes: cellarray e.g. {'tt' [2 0] [2 1]} or {'rep' 1 2}
            if ~exist('kategories','var') || isempty(kategories)
                kat = unique(obj.P.data(:,obj.P.sloupce.kategorie))'; % all category numbers - ciselne vyjadreni kategorie podnetu 0-n + transpose to make it row
            elseif cell2double(kategories(1)) >= 0
                kat = cell2double(kategories); %selected categories                
            else
                kat = -1;
            end   
            if kat(1)>=0
                
                if exist('trialtypes','var') && ~isempty(trialtypes) && size(obj.P.data,1) < size(obj.trialtypes,1) %epochs are missing for this subject due to some error
                    nepochs =     size(obj.trialtypes,1) -     size(obj.P.data,1);
                    obj.P.data = cat(1,obj.P.data,repmat(-1, nepochs, size(obj.P.data,2))); %add empty epochs to behavioral data, filled with -1
                end
                if ~exist('trialtypes','var') || isempty(trialtypes) || numel(trialtypes) <=2 %no or just one trialtype                   
                    rt = nan(size(obj.P.data,1),numel(kat)); %pole kam budu ukladat reakcni casy v sekundach
                    for k = 1:numel(kat) %funguje jen pro radky, kat je sloupec
                        ikat = obj.P.data(:,obj.P.sloupce.kategorie)==kat(k);
                        if exist('trialtypes','var') && ~isempty(trialtypes) %filter epochs for just this trialtype
                            itt = obj.GetTrialType(1,trialtypes)'; %1 x epochs, index of epochs with this trialtype/repetition                         
                        else
                            itt = true(size(obj.P.data,1),1);  %all epochs , in one column             
                        end
                        rt(1:sum(itt & ikat),k) = 24*3600*(obj.P.data(itt & ikat,obj.P.sloupce.ts_odpoved) - obj.P.data(itt & ikat,obj.P.sloupce.ts_podnet));
                        %rt: trials x kats - when kategories have different no of trials, the rest for the columnt is NaN
                    end                               
                else % get reaction times for different trialtypes
                    rt = nan(size(obj.P.data,1),numel(trialtypes)-1); %pole kam budu ukladat reakcni casy v sekundach
                    ikat = ismember(obj.P.data(:,obj.P.sloupce.kategorie),kat);
                    for t = 1:numel(trialtypes)-1
                        if strcmp(trialtypes{1},'tt')
                            itt = obj.trialtypes{:,trialtypes{t+1}(1)} == trialtypes{t+1}(2); %eg. {'tt' [1 0] [1 1]}
                        elseif strcmp(trialtypes{1},'rep')
                            itt = obj.P.data(:,obj.P.sloupce.opakovani)==trialtypes{t+1}; %eg. {'rep' 1 2} 
                        else
                            itt = true(size(ikat));
                        end
                        rt(1:sum(itt & ikat),t) = 24*3600*(obj.P.data(itt & ikat,obj.P.sloupce.ts_odpoved) - obj.P.data(itt & ikat,obj.P.sloupce.ts_podnet));                                                    
                    end
                end
                rt = rt(any(~isnan(rt),2),:); % necha jen radky, kde je nejake ~NaN cislo
            else
                rt = 24*3600*(obj.P.data(:,obj.P.sloupce.ts_odpoved) - obj.P.data(:,obj.P.sloupce.ts_podnet));
            end
            if any(sum(rt<0)>=1) && isempty(obj.warning_rt)
                warning('Pozor: %d zaporne reakcni casy',sum(rt<0));
                obj.warning_rt = true;
            end
        end
        
        function isi = InterStimulusInterval(obj)
            %vraci matici vsech interstimulus intervalu v sekundach
            isi = 24*3600*(obj.P.data(2:end,obj.P.sloupce.ts_podnet) - obj.P.data(1:end-1,obj.P.sloupce.ts_podnet));
        end
        
        function ts_podnety = TimeStimuli(obj,response)
            %returns the array of timestamps of all stimuli, if response=1 returns array of responses
            %ts_podnety, size n_epochs x 1
            if exist('response','var') && response==1
                ts_podnety = obj.P.data(:,obj.P.sloupce.ts_odpoved); %responses
            else
                ts_podnety = obj.P.data(:,obj.P.sloupce.ts_podnet); %stimuli
            end
        end
                       
        function [kat, katnum] = Category(obj,event)
            %vrati 1 retezec s popisem kategorie eventu a 2. cislo kategorie
            katnum = obj.P.data(event,obj.P.sloupce.kategorie); %cislo kategorie
            kat = obj.P.strings.podminka{katnum+1};            
        end 
        
        function [katnum, katstr] = Categories(obj,tisk,Wp)
            %return numbers and names of all categories
            %katnum are category numbers (0-n), equivalent to PsyData.P.strings.podminka
            %katstr are category names, also from PsyData.P.strings.podminka, in the same order
            %if Wp.trialtypes is not empty, returns trialtypes instead of categories
            
            if ~exist('tisk','var'), tisk = 0; end %dfaultne netisknu kategorie
            if ~exist('Wp','var') 
                katnum = cell2mat(obj.P.strings.podminka(:,2))'; %kategorie v radku
                katstr = cell(size(katnum));
                for k = 1:numel(katnum)
                    katstr{k} = obj.P.strings.podminka{k,1};
                    if tisk
                        disp([ num2str(katnum(k)) ': ' katstr{k}]);
                    end
                end     
            else %categories accordig to Wp 
                if isfield(Wp,'trialtypes') && numel(Wp.trialtypes) >= 3 %for more than two trialtypes, make contrasts between them
                    katnum = 0:numel(Wp.trialtypes)-2; %numbers starting from 0
                    katstr = cell(size(katnum));                    
                    for k = 1:numel(katnum)
                        katstr{k} = [obj.trialtypes.Properties.VariableNames{Wp.trialtypes{k+1}(1)} '=' num2str(Wp.trialtypes{k+1}(2))];
                    end                   
                else %for only one or no trialtype, make contrasts according to wp.kats
                    katnum = Wp.kats; %katnum are numbers in obj.P.strings.podminka{:,2}, not indexes of obj.P.strings.podminka{:,1} !
                    katstr = cell(size(katnum));
                    for k = 1:numel(katnum)                                                
                        katstr{k} = obj.CategoryName(cellval(katnum,k));
                    end
                end
            end
        end
        function [kat] = CategoryName(obj,katnum,concat)
            %vraci jmeno kategorie z cisla, jmeno kategorie se pocita od 0
            %muze byt i vic cisel kategorii
            % pak se jmena spoji do jednoho retezce pomoci parametru concat (default +)
            % pokud je concat prazdne, vrati se cell array
            
            if ~exist('concat','var'), concat = '+'; end %defaulte se vice jmen kategorii spojuje pomoci +            
            if numel(katnum) == 1               
               kat = obj.P.strings.podminka{katnum+1};               
            elseif ~isempty(concat)
                kat = '';
                for k = katnum
                    if numel(kat) == 0
                        kat = obj.P.strings.podminka{k+1};
                    else
                        kat = [ kat '+' obj.P.strings.podminka{k+1}]; %#ok<AGROW>
                    end
                end
            else %chci jmena kategorii vratit jako cell array
                kat = cell(1,numel(katnum));
                for k  = 1:numel(katnum)
                    if numel(cellval(katnum,k))>1
                        katkat = cell(numel(katnum{k}),1);
                        for kk = 1:numel(katnum{k})
                            katkat{kk} = obj.P.strings.podminka{katnum{k}(kk)+1};
                        end
                        kat{k} = cell2str(katkat);
                    else
                        kat{k} = obj.P.strings.podminka{cellval(katnum,k)+1};
                    end
                end
            end
        end
        function [katnum] = CategoryNum(obj,katname)
            %vraci cislo kategorie podle jmena
            %muze byt i nekolik jmen kategorii v cell array
            katnum = nan(size(katname));
            for k = 1:numel(katname)
                ikat = not(cellfun('isempty',strfind(obj.P.strings.podminka(:,1),katname{k}))); %logical index, ale s jen jednim true
                katnum(k) = obj.P.strings.podminka{ikat,2};
            end
        end
        function [blocks, srate, test, kategorie]= GetBlocks(obj)
            %vraci promenne pro bloky o stejne kategorii
            b1 = [find(obj.P.data(2:end,obj.P.sloupce.kategorie) ~= obj.P.data(1:end-1,obj.P.sloupce.kategorie)); size(obj.P.data,1)]; %konce bloku
            if min(diff(b1)) == 1
                b1 = (1:size(obj.P.data,1))'; %pro PPA test - pokud jsou nektere bloky velikosti 1, tak je udelam vsechy tak mal
            end %jinak to dela problemy pri prohazeni a pocitam srate
            b0 = [ 1; (b1(1:end-1)+1)]; %zacatky bloku
            kategorie = obj.P.data(b0,obj.P.sloupce.kategorie);            
            blocks = [b0 b1];
            if isempty(obj.blocks) || size(blocks,1) ~= numel(obj.blocks.srate)
                obj.blocks.srate = zeros(numel(b0),1); % prumerna uspesnost za blok
                obj.blocks.test = ones(numel(b0),1); %jestli by blok testovy, tady nastavim ze vsechny
                for block = 1:numel(b0)
                    obj.blocks.srate(block,1)=mean(obj.P.data(b0(block) : b1(block),obj.P.sloupce.spravne));                       
                    obj.blocks.test(block) = obj.blocks.test(block) - max(obj.P.data(b0(block) : b1(block),obj.P.sloupce.zpetnavazba));                        
                        %test=1 vsechny epochy v bloku testove, test=0 alespon jedna epocha treningova
                end                
            end
            srate = obj.blocks.srate;
            test = obj.blocks.test;
            
        end
        
        function [resp,rt,category,test] = GetResponses(obj,Wp)
            %returns the subjects responses - correct=1/incorrect=0, reaction time (from PsychoPy), stimulus category, test=1/training=0
            if ~exist('Wp','var'), Wp = struct(); end
            S = obj.P.sloupce;
            if ~isfield(Wp,'trialtypes') || isempty(Wp.trialtypes) || ~iscell(Wp.trialtypes)
                resp = obj.P.data(:,S.spravne); % spravnost odpovedi
                if strcmp(obj.testname,'ppa') %correction of correct responses for PPA test
                    klavesa = obj.P.data(:,S.klavesa); % spravnost odpovedi
                    iovoce = obj.P.data(:,S.kategorie) == 0; %ovoce
                    iklavesa = klavesa == 1; %pressed space bar
                    if sum(iovoce & iklavesa)>0 %if he responded at least ones for ovoce
                        if mean(resp(iovoce & iklavesa) < 0.5) %if correct response is marked by 0
                           resp(iovoce) = 1- resp(iovoce); %reverse all ovoce reponses: 1=correct
                        end
                    else
                        iklavesa = klavesa == -1; %index of no responses
                        if mean(resp(iovoce & iklavesa) > 0.5) %if incorrect reps is marked by 1
                            resp(iovoce) = 1- resp(iovoce); %reverse all ovoce reponses: 1=correct
                        end
                    end
                    %tab = [resp(iovoce & iklavesa), klavesa(iovoce & iklavesa)];
                end
                rt = obj.P.data(:,S.rt); % reakcni casy
                category = obj.P.data(:,S.kategorie); %kategorie
                test = obj.P.data(:,S.zpetnavazba)==0; %vratim index testovych trialu
            elseif strcmp(Wp.trialtypes{1}, 'tt') 
                assert(~strcmp(obj.testname,'ppa'),'CPsyData.GetResponses cannot now return PPA test trialtypes'); %todo - similar correction as above
                lines = 0;
                if numel(Wp.trialtypes) <= 2 %only one trialtype => make contrast between stimulus categoires
                        kats = Wp.kats; %kategories in rows, for each one cycle below
                    else  %at least two trials types + 'tt'. It means we are contrasting these trialtypes
                        kats = Wp.kats'; %kategorie in column, so treat them as one, only one cycle below 
                end   
                if size(obj.P.data,1) < size(obj.trialtypes,1) %epochs are missing for this subject due to some error
                    nepochs =     size(obj.trialtypes,1) -     size(obj.P.data,1);
                    obj.P.data = cat(1,obj.P.data,repmat(-1, nepochs, size(obj.P.data,2))); %add empty epochs to behavioral data, filled with -1. They are not used later, as their condition is -1
                end
                categnum = 1; %maximal number of categories in each (combined) category
                for t = 2:numel(Wp.trialtypes) 
                    for kat = 1:size(kats,2) %cycle over columns of kats
                        itt = obj.trialtypes{:,Wp.trialtypes{t}(1)} == Wp.trialtypes{t}(2); %index of epochs with this trialtype ((column) == val)                    
                        ikat = ismember(obj.P.data(:,S.kategorie),cellval(kats(:,kat)));  %index of epochs with any of these categories in this column
                            %is member seems to work with the whole matrix, so independed on columns/rows for second argument
                        lines =  lines + sum(itt & ikat ); %number of trials/epochs with this trialtype
                        categnum = max(categnum,numel(cellval(kats(:,kat))));
                    end                    
                end
                resp = zeros(lines,1);
                rt = zeros(lines,1);
                category = zeros(lines,categnum); %more than one column for combined categories
                test = zeros(lines,1);
                lines = 0;
                for t = 2:numel(Wp.trialtypes) 
                    for kat = 1:size(kats,2) %cycle over columns of kats
                        %there are only test trials in CM.PsyData.Pmulti
                        itt = obj.trialtypes{:,Wp.trialtypes{t}(1)} == Wp.trialtypes{t}(2);
                        ikat = ismember(obj.P.data(:,S.kategorie),cellval(kats(:,kat)));  %index of epochs with any of these categories
                        resp(lines+1:lines+sum(itt & ikat)) = obj.P.data(itt & ikat,S.spravne); %correcness
                        rt(lines+1:lines+sum(itt & ikat)) = obj.P.data(itt & ikat,S.rt); %reaction time
                        if numel(Wp.trialtypes) == 2 %if there is only one trialtype, so contrast categories
                            categ = cellval(kats(kat));
                        else 
                            categ = t-2; %contrast trialtypes
                        end
                        category(lines+1:lines+sum(itt & ikat),1:numel(categ)) = repmat(categ,sum(itt & ikat),1); %not original stimulus category, but number of trialtype category, should start from 0
                        test(lines+1:lines+sum(itt & ikat)) = obj.P.data(itt & ikat,S.zpetnavazba)==0; %feedback 1=training,0=test
                        lines = lines + sum(itt & ikat);
                    end
                end
            else
                error(['CPsyData.GetResponses: not defined trialtype ' Wp.trialtypes{1}]);
            end
        end
        
        function chyby = GetErrorTrials(obj)
            %vrati pole indikujici chybu/vyrazeni pro kazdy trial=radky
            %sloupce: jednotlive chyby, chybne bloky, treningovy trial, prilis kratky reakcni cas
            %chybny blok ma < 75% uspesnost
            S = obj.P.sloupce;          
            
            chyby = zeros(size(obj.P.data,1),4); %ctyri sloupce - chybne trials a chybne bloky, trening, prilis rychle reakcni casy
            chyby(:,1) = obj.P.data(:,S.spravne)==0; %pro PPA jsou vsechny ovoce spatne. Sloupec spravne je u ovoce vzdy 0, chyba v PHP asi
            rt = obj.ReactionTime(-1); %reakcni casy podle Sychropulsu - do not distinguish stimulus categories
            rtPsy = obj.P.data(:,S.rt); %reakcni cas podle psychopy
            chyby(:,4) = rt(:,1) < 0.1 | (rtPsy(:,1) < 0.1 & rtPsy(:,1) > 0);  %v PPA clovek nereaguje spravne, takze 0 jako cas odpovedi me nezajima 
                    %chyba, pokud je reakcni cas prilis kratky (0 v PsychoPy znamena, ze nereagoval, to je taky chyba)
            [blocks,srate,blocktest]=obj.GetBlocks();            %#ok<PROP>
            for b = 1:size(blocks,1) %#ok<PROP>
                if srate(b) < 0.75  %chybny blok
                    chyby(blocks(b,1) : blocks(b,2) , 2) = ones( blocks(b,2) - blocks(b,1) +1,1);%#ok<PROP> %vyplnim jednickami
                end 
                if  blocktest(b)==0 %treningovy blok
                    chyby(blocks(b,1) : blocks(b,2) , 3) = ones( blocks(b,2) - blocks(b,1) +1,1); %#ok<PROP> %vyplnim jednickami
                end
            end
        end
        function [iTrialTypeCh, TrialTypeCh] = GetTrialType(obj,els,trialtype)
            %returns index of epochs for each channel with this repetition/trialtype
            %els is copy of E.els from CiEEGData    
            %trialtype: cellarray e.g. {'tt' [2 0]} or {'rep' 1}, ignores trialtype{3} 
            %TrialTypeCh - trialtype or stimulus repetition number for each channels x epoch 
            %iTrialTypeCh - bool values for trialtype{2} channels x epoch 
            if ~exist('els','var') || isempty(els), els = 1; end %get trialtypes for only one subject
            assert(numel(trialtype)>=2,'GetTrialType: trialtype needs to be cell array with 2 values');            
            if strcmp(trialtype{1}, 'tt')
                assert(~isempty(obj.trialtypes),'no trialtypes loaded');                    
                assert(size(obj.trialtypes,2)>= trialtype{2}(1) & isa(obj.trialtypes{1,trialtype{2}(1)},'double'), 'wrong trialtype column number or column type');
            end
            S = obj.P.sloupce;
            TrialTypeCh = zeros(els(end),obj.GetEpochsMax()); %channels x epochs
            if ~isa(obj,'CPsyDataMulti') 
                els = els(end); % only one pacient, the repeat is same for all channels
            end
            els0 = [1 , els(1:end-1)+1]; %start channel of each subject, or only [1] when there is only one subject
            iSbak = obj.iS; %backup the selected pacient
            for s = 1:numel(els) %over all subjects
                if numel(els) > 1, obj.SubjectChange(s); end %for one subject use the current one = do not change it
                if strcmp(trialtype{1},'rep') %repetitions can be different for different subjects
                    if isfield(S,'opakovani_obrazku')
                        opakovani = obj.P.data(:,S.opakovani_obrazku); % ppa test - array of repetitions number for each epoch
                    elseif isfield(S,'opakovani')
                        opakovani = obj.P.data(:,S.opakovani);  %aedist test
                    else
                        opakovani = zeros(size(obj.P.data,1),1); %opakovani 0
                    end
                    TrialTypeCh(els0(s):els(s),1:numel(opakovani)) = repmat(opakovani',els(s)-els0(s)+1,1);
                elseif strcmp(trialtype{1}, 'tt')                    
                    tt = obj.trialtypes{:,trialtype{2}(1)};                  
                    TrialTypeCh(els0(s):els(s),1:numel(tt)) = repmat(tt',els(s)-els0(s)+1,1); %repetition of tt for all channels of this pacient
                else 
                    error('GetTrialType: unknown trialtype type');
                end                
            end
            obj.SubjectChange(iSbak);
            if strcmp(trialtype{1},'rep') && numel(trialtype)>=2
               iTrialTypeCh = ismember(TrialTypeCh , trialtype{2}); 
            elseif strcmp(trialtype{1},'tt') && numel(trialtype{2})>=2   
               iTrialTypeCh = ismember(TrialTypeCh, trialtype{2}(2));  %channels x epochs
            else
               iTrialTypeCh = []; 
            end            
        end
        function filtered = FilteredIn(obj,epochs, filter)      
            %returns array of 0/1 for each epoch: 1 means the epochs meets the filter condition, ie there is the looked for value in a given column
            %epochs - numbers of epochs to be checked
            %filter - cellarray filter{1} - column number in obj.P.data, filter{2} - looked for values in this column
            %filtered - array numel(epochs) x 1;
            filtered = 0;
            if isempty(filter)
                filtered = ones(size(epochs)); %when no filter given, all epochs are filtered in 
            elseif max(epochs) <= size(obj.P.data,1) && filter{1} > 0 && filter{1}<=size(obj.P.data,2)
                filtered = ismember(obj.P.data(epochs,filter{1}),filter{2});                    
            end
            
        end
        function agree = CategoriesAgreeWithStat(obj,kategories,Wp)
            %return true if the behavioral kategories agree with kategories in stat
            kategoriesVector = cell2double(kategories); %kategories in one vector, even for combined categories
            katnum = obj.Categories([],Wp);            
            katnumVector =  cell2double(katnum);
            agree = sum(ismember(kategoriesVector,katnumVector))==numel(kategoriesVector);    
        end
        function obj = Cond2Epochs(obj)
            %predela conditions na epochy, kvuli ITPC
            katnum = obj.Categories();
            data2 = zeros(numel(katnum),size(obj.P.data,2));
            for k = 1:numel(katnum)
                idata = obj.P.data(:,obj.P.sloupce.kategorie)==katnum(k) & obj.P.data(:,obj.P.sloupce.zpetnavazba)==0 ;
                data2(k,:) = [0 0 1 mean(obj.P.data(idata,4)) 0 0 katnum(k) min(obj.P.data(idata,8)) min(obj.P.data(idata,9))];
            end
            obj.P.data = data2;
        end
        function obj = LoadTrialTypes(obj)
            if strcmp(obj.testname,'aedist')
                Tname = 'AedistTrialTypes'; %d:\eeg\motol\scripts\AedistTrialTypes.mat                 
            elseif strcmp(obj.testname,'menrot')
                Tname = 'MenrotTrialTypes';                               
            else
                disp(['trial types for ' obj.testname 'not defined']);                
            end    
            if exist('Tname','var')
                T = load(Tname);
                assert(size(T.(Tname),1) >= size(obj.P.data,1),['wrong number of epochs in ' Tname ', should be ' num2str(size(obj.P.data,1))]);                                   
                trialsdiff = size(T.(Tname),1) - size(obj.P.data,1); %how larger is trialtype table than trial data table
                obj.trialtypes = T.(Tname)(trialsdiff+1:end,:);
                disp(['last ' num2str(size(obj.P.data,1)) ' trialtypes loaded from ' Tname]);
            end
            
        end
        %% PLOT FUNCTIONS
        function [obj, chyby] = PlotResponses(obj)
            %plots reponses times of all trials and blocks,including errors and success rate
            %blue circles - response time,green line=conditions, 
            %nakresli graf rychlosti vsech odpovedi a bloku, vcetne chyb a uspesnosti za blok
            
            S = obj.P.sloupce;
            test = obj.P.data(:,:); %vyberu vsechny trialy
            if isprop(obj,'fhR') && ~isempty(obj.fhR) && ishandle(obj.fhR) 
                figure(obj.fhR); %pokud uz graf existuje, nebudu tvorit znova
                clf; %smazu aktualni figure
            else                
                obj.fhR = figure('Name','PsyData.PlotResponses');
            end
            plot(test(:,S.rt),'-o'); %plot responses time of all trials
            hold on;
            treningtrials = find(obj.P.data(:,S.zpetnavazba)==1);            
            bar(treningtrials,test(treningtrials,S.rt),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]); %oznacim treningove pokusy
            
            chyby = find(test(:,S.spravne)==0);
            plot(chyby,test(chyby,S.kategorie),'sr','MarkerSize',10,'MarkerFaceColor','y'); %vykreslim chyby na vysce kategorie
            plot(chyby,test(chyby,S.rt),'or','MarkerFaceColor','y'); %vykreslim reakcni cas cervene
            
            plot(test(:,S.kategorie),'g','LineWidth',2); %vykreslim kategorie podnetu, jako zelenou caru
            for j = 1:size(obj.P.strings.podminka,1)
                text (10,obj.P.strings.podminka{j,2}+0.05,obj.P.strings.podminka{j,1},'Color','r','FontSize',15);
            end
                        
            [blocks,srate,blocktest]=obj.GetBlocks();
            for b = 1:size(blocks,1)
                if blocktest(b) == 1 %pokud se jedna o testovy blok
                    procent = round(srate(b)*100);
                    if procent < 75, color = 'r'; else color = 'b'; end %cervenou barvou, pokud je uspesnost pod 75%
                    if procent < 75 || size(blocks,1)<100 %vypisuju uspesnost u vsech bloku jen pokud jich je min nez 100; jinak jen chybne
                        text( blocks(b,1) , -0.1, [ num2str(procent) '%'],'FontSize',8,'Color',color);
                    end
                elseif size(blocks,1) < 100
                    text( blocks(b,1) , -0.1, '*','FontSize',8,'Color','k'); %treningove bloky jako hvezdicky
                end
            end
            resp = obj.P.data(:,S.spravne); %vratim i reakcni casy            
            text(0,-0.2, ['pocet chyb: ' num2str(numel(resp)-sum(resp)) ]);
        end
        function [obj] = PlotITI(obj)
            %graf intervalu mezi obrazky
            %TODO korelace mezi RT podle Synchropulsu a PsychoPy
            figure('Name','ITI');
            subplot(1,2,1);
            isi = obj.InterStimulusInterval();
            plot(isi,'.');
            ylim([1 prctile(isi,99)]); %rozsah do 99  percentilu - odstranim odlehle hodnoty
            title('ITI');
            ylabel('sec');
            
            subplot(1,2,2);
            rt = obj.ReactionTime(); %ctyri sloupce reakcnich casu podle synchropulsu - all categories
            numkatvals = size(rt,1);
            numkats = size(rt,2);
            rt = rt(:);      % jeden sloupec, v poradi ovoce, 
            plot(rt,'o');
            for j = 1:numkats
                line([numkatvals*j numkatvals*j],[0.4 1.4]);
                text(numkatvals*(j-0.5),0.1,obj.P.strings.podminka{j},'Color','blue');
            end
            title('RT');
            ylabel('sec');
            hold on;
            rtPsy = obj.P.data(:,obj.P.sloupce.rt);
            plot(rtPsy,'xm'); %reakcni casy podle PsychoPy - tohle neni serazene podle kategorii je tam 0 pokud zadna reakce
            
            %druhy graf - 21.6.2017
            figure('Name','ITI - Synchropuls vs PsychoPy');
            rt = obj.ReactionTime(-1); %no categories distinguished
            plot(rt,rtPsy,'o');
            xlabel('RT Synchropulse');
            ylabel('RT PsychoPy');
       
        end
        function [id] = PacientID(obj,full)
            %funkce ktera bude vracet id pacienta
            %jestli chci plne id nebo jen obj.P.pacientid
            id = obj.P.pacientid;
            if ~exist('full','var'), full = true; end
            if full && length(id) <= 5 %pokud id obsahuje jen napr p132
                podtrzitko = strfind(obj.P.eegfile,'_');
                id = [id '-' obj.P.eegfile(1 : podtrzitko-1)]; %pritam jeste zacatke jmena eeg souboru, naprilad VT18
             end
        end
        function [obj]= RemoveTraining(obj)
            %odstrani treningove trialy - ty kde byla zpetna vazba
            %since 15.12.2017
            test = obj.P.data(:,obj.P.sloupce.zpetnavazba)==0;
            obj.P.data = obj.P.data(test,:);
        end
        function epochs = GetEpochsMax(obj)
            %returns the number of epochs for current subject. Created only for overloading by CPsyDataMulti            
            epochs=size(obj.P.data,1);
        end
        function [name] = TrialTypeName(obj,trialtype)
            %returns name of trialtype/repetition
            %trialtype is cellarray, e.g. {'tt' [2 0]} or {'rep' 1}
            if strcmp(trialtype{1},'rep')
                name = ['Repeat_' num2str(trialtype{2})];
            elseif strcmp(trialtype{1},'tt')
                name = ['TT_' strrep(obj.trialtypes.Properties.VariableNames{trialtype{2}(1)},'_','_') '=' num2str(trialtype{2}(2))];
            else
                name = cell2str(trialtype);
            end
        end
        
    end
    methods (Access = protected)
        function [obj] = DoplnZpetnavazba(obj)
            %doplni sloupec zpetnavazba, pokud neexistuje a naplni ho nulama
            if isstruct(obj.P) && ~isfield(obj.P.sloupce,'zpetnavazba')                
                obj.P.data = [ obj.P.data zeros(size(obj.P.data,1),1)]; %doplnim dalsi sloupec do dat
                obj.P.sloupce.zpetnavazba = size(obj.P.data,2); %pojmenuju ho zpetnavazba
            end
        end
    end 
end

