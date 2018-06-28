files = {   'AEdist CHilbert 50-150 -0.5-1.2 refBipo Ep2018-04_CHilb.mat',...
            'AEdist CMorlet 1-10M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat',...
            'AEdist CMorlet 4-8M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat',...
            'AEdist CMorlet 1-4M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat'};       
 
kategorie = {4,  'EgoXControl';5,  'AlloXControl';6,  'AlloXEgo'};  %jmena kategorii (nazvy extraktu) a jejich cisla v CB.PAC
testname = 'aedist';
[ pacienti, setup  ] = pacienti_setup_load( testname );
pocetcyklu = numel(files) * size(kategorie,1) ;
kontrast = 1; %statisticky kontrast
stat = setup.stat_kats{kontrast}; %resp setup.stat_kats{1} {2} nebo {3} pro menrot 
overwrite_extracts = 1; %jestli se maji prepisovat extrakty pro kazdeho pacienta
overwriteCM = 0; %jestli se maji prepisovat soubory CHilbertMulti
FFfilenames = cell(numel(files),size(kategorie,1));
FFFilenames_XLS = cell(1+pocetcyklu*numel(pacienti),5);
FFFilenames_XLS(1,:)={'file','fileno','kat','katno','extract'};
pocetextracts = 1;
cyklus = 1;

logfilename = ['logs\BatchExtract_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.log'];
FFFilenames_logname = ['logs\BatchExtractFilenames_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.xls'];
[fileID,~] = fopen(logfilename,'wt'); %soubor na logovani prubehu

tablelog = cell(pocetcyklu+1,6); %frekvence, soubor, reference, status, chyba -  z toho bude vystupni xls tabulka s prehledem vysledku
tablelog(1,:) = {'file','fileno','kategorie','result','file','datetime'}; %hlavicky xls tabulky

for f = 1:numel(files)
    try
        msg = [' --- ' files{f} ': IntervalyResp *********** '];
        disp(msg); fprintf(fileID,[ msg '\n']);
        
        CB = CBrainPlot;    
        CB.IntervalyResp(testname,[0.1 1],files{f},kontrast); %ziskam signif rozdily pro kategorie a mezi kategoriemi pro vsechny pacienty       
        for kat = 1:size(kategorie,1)
            try
                outfilename = ['d:\eeg\motol\pacienti\0sumarne\CM ' kategorie{kat,2} ' ' files{f}]; %jmeno souboru CHilbertMulti
                if exist(outfilename,'file')==2 && overwriteCM == 0
                    msg = [ ' --- ' strrep(outfilename,'\','\\') ' NEULOZENO, preskoceno '  datestr(now)];
                    disp(msg); fprintf(fileID,[ msg '\n']);                     
                    tablelog(cyklus+1,:) = { files{f}, num2str(f), kategorie{kat,2}, 'preskoceno', outfilename,datestr(now) };                     
                else
                    msg = [ ' --- ' strrep(outfilename,'\','\\') ' zpracovavam '  datestr(now)];
                    disp(msg); fprintf(fileID,[ msg '\n']);  
                    
                    CM = CHilbertMulti;
                    %vytvorim extrakty podle tabulky PAC, pro vsechny pacienty a pro tuto kategorii
                    filenames_extract = CM.ExtractData(CB.PAC{1,kategorie{kat,1}},'aedist',files{f},kategorie{kat,2},overwrite_extracts);

                    FFfilenames{f,kat} = filenames_extract; 
                    FFFilenames_XLS(pocetextracts:pocetextracts+numel(filenames_extract)-1,:) = ...
                        cat(2,repmat({files{f},f,kategorie{kat,2},kat},numel(filenames_extract),1),filenames_extract);
                    pocetextracts = pocetextracts + numel(filenames_extract);                

                    CM.ImportExtract(filenames_extract);
                    CM.ResponseSearch(0.1,stat);  

                    CM.Save(outfilename);

                    msg = [ ' --- ' files{f} ': ' kategorie{kat,2} ' OK '  datestr(now)];
                    disp(msg); fprintf(fileID,[ msg '\n']);            
                    tablelog(cyklus+1,:) = { files{f}, num2str(f), kategorie{kat,2}, 'saved', outfilename,datestr(now) }; 
                end
            catch exception 
                errorMessage = sprintf('** Error in function %s() at line %d.\nError Message:\n%s', ...
                    exception.stack(1).name, exception.stack(1).line, exception.message);                            
                disp(errorMessage);  fprintf(fileID,[errorMessage '\n']);  %#ok<DSPS> %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal                                        
                tablelog(cyklus+1,:) = { files{f}, num2str(f), kategorie{kat,2}, 'error', exception.message , datestr(now)}; 
                clear CM; 
            end 
            cyklus = cyklus + 1;
            xlswrite([logfilename '.xls'],tablelog); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
            xlswrite(FFFilenames_logname,FFFilenames_XLS); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
        end
    catch exception
        errorMessage = sprintf('** Error in function %s() at line %d.\nError Message:\n%s', ...
                exception.stack(1).name, exception.stack(1).line, exception.message);                            
        disp(errorMessage);  fprintf(fileID,[errorMessage '\n']);  %#ok<DSPS> %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal                                        
        tablelog(cyklus+1,:) = { files{f}, num2str(f), 'no kat', 'error', exception.message , datestr(now)}; 
        clear CB;        
        cyklus = cyklus + 1;
        xlswrite([logfilename '.xls'],tablelog); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
        xlswrite(FFFilenames_logname,FFFilenames_XLS); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
    end
end

%system('shutdown -h') 