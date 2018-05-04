classdef CFiles < handle
    %CFILES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Access = public)  
        function filenames = FindFiles(obj,testname, mask)
            %najde vsechny CHilb soubory pro konkretni test, a udela xls tabulku s parametry ve sloupcich
            if ~exist('mask','var') || isempty(mask) , mask = '*'; end 
            [pacienti,setup] = pacienti_setup_load( testname );           
            filenames = cell(0,14); %budu vrace vsechny nalezene udaje, nejen filename, kvuli prehledu
            filenames(1,:) = {'filename','folder','datetime','filesize','pacientno','pacient','datatype','frekvence',... 
                              'epochtime','reference','name','align','brainlabel','suffix'};
            for p = 1:numel(pacienti)
               path = [setup.basedir pacienti(p).folder filesep setup.subfolder filesep];
               files = dir([path mask '.mat']);
               files_cell = struct2cell(files)';
               files_cell = cat(2,files_cell(:,1:4),cell(size(files_cell,1),10));                
               for f = 1:size(files_cell,1)
                   files_cell{f,5} = p;
                   files_cell{f,6} = pacienti(p).folder;
                   parts = obj.FileParts(files_cell{f,1});                   
                   files_cell(f,7:14) = parts;
               end
               filenames = cat(1,filenames,files_cell);
            end
            disp(['nalezeno ' num2str(size(filenames,1)-1) ' souboru pro ' num2str(numel(pacienti)) ' pacientu']);
            outputfilename = [setup.basedir testname '_FindFilesOutput.xlsx'];
            xlswrite(outputfilename,filenames);
            disp(['zapsano do ' outputfilename]);            
        end    
        
    end
    methods (Static,Access = public) 
        function parts = FileParts(filename) 
            %vraci cell array s udaji o souboru podle casti filename, celkem 8 prvku
            parts = cell(1,8); %DataType freq time ref jmeno zarovnani label zbytek            
            ispace = find(filename==' ');            
            idot = find(filename == '.',1,'last');
            isub = find(filename == '_');
            if numel(isub)==0, isub = idot(1); end 
            if numel(ispace)<4
                for p = 1:numel(ispace)
                    if p < numel(ispace)
                        parts{p} = filename(ispace(p)+1: ispace(p+1)-1); %CHilbert, CMorlet nebo CiEEG
                    else
                        parts{p} = filename(ispace(p)+1: idot(1)-1); %ERP nebo rozsah frekvenci
                    end
                end                
            else
                timelimits =  filename(ispace(3)+1: ispace(4)-1); %casove rozmezi
                if numel(find(timelimits=='-'))==0 %pokud se nejedna o casove rozmezi
                    ispace (end+1)= ispace(end); %prodlouzim o jedno
                    ispace(4:end) = ispace(3:end-1); % rozsah casoveho rozmezi udelam 0 = bude prazdne
                end                            
                for s = 1:numel(ispace)
                    parts{1} = filename(ispace(1)+1: ispace(2)-1); %CHilbert, CMorlet nebo CiEEG
                    parts{2} = ['''' filename(ispace(2)+1: ispace(3)-1)]; %ERP nebo rozsah frekvenci, na zacatku apostrof aby to excel cetl dobre
                    parts{3} = ['''' filename(ispace(3)+1: ispace(4)-1)]; % castove rozmezi -1.4-0.3, na zacatku apostrof aby to excel cetl dobre
                    if numel(ispace)>4
                        parts{4} = filename(ispace(4)+1: ispace(5)-1); % refEle aj
                    end
                    if numel(ispace)>5
                        name = filename(ispace(5)+1: ispace(6)-1); %treba Ep2018-04RES
                        label = filename(ispace(6)+1: isub(1)-1); %treba RSCpariet
                    else
                        name = filename(ispace(5)+1: isub(1)-1); %treba Ep2018-04RES
                        label = '';
                    end
                    if numel(name)>=3 && (strcmp(name(end-2:end),'RES') || strcmp(name(end-2:end),'PCO'))
                        parts{5} = name(1:end-3);
                        parts{6} = name(end-2:end);                    
                    else
                        parts{5} = name(1:end);
                        parts{6} = 'stimul'; %zarovnani podle podnetu
                    end
                    parts{7} = label;
                    parts{8} = filename(isub(1)+1:idot(1)-1);    
                end
            end
        end
    end
    
end

