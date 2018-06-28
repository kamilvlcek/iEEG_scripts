files = {   'AEdist CHilbert 50-150 -0.5-1.2 refBipo Ep2018-04_CHilb.mat',...
            'AEdist CMorlet 1-10M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat',...
            'AEdist CMorlet 4-8M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat',...
            'AEdist CMorlet 1-4M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat'};       
 
%kategorie = {4,  'EgoXControl';5,  'AlloXControl';6,  'AlloXEgo'};  %jmena kategorii (nazvy extraktu) a jejich cisla v CB.PAC
testname = 'aedist';
[ pacienti, setup  ] = pacienti_setup_load( testname );
kontrasts = 1:numel(setup.stat_kats); %statisticky kontrast

pocetcyklu = 0;
for kontrast = 1:numel(kontrasts) %cyklus jen na vypocet celkoveho poctu cyklu pres vsechny kontrasty ve statistice
    stat = setup.stat_kats{kontrast}; %resp setup.stat_kats{1} {2} nebo {3} pro menrot 
    kombinace_kat = combinator(length(stat),2,'p'); %permutace bez opakovani z poctu kategorii
    kombinace_kat = kombinace_kat(kombinace_kat(:,1)>kombinace_kat(:,2),:); %vyberu jen permutace, kde prvni cislo je vetsi nez druhe     
    pocetcyklu = pocetcyklu + numel(files) * size(kombinace_kat,1) ;
end

overwrite_extracts = 1; %jestli se maji prepisovat extrakty pro kazdeho pacienta
overwrite_brainplots = 0;
overwriteCM = 0; %jestli se maji prepisovat soubory CHilbertMulti
cyklus = 1;
pocetextracts = 1;

%log soubory

%1. seznam vsech extraktu
logfilename = ['logs\BatchExtract_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.log'];
FFFilenames_logname = ['logs\BatchExtractFilenames_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.xls'];
FFFilenames_XLS = cell(1+pocetcyklu*numel(pacienti),5); 
FFFilenames_XLS(1,:)={'file','fileno','kat','katno','extract'};
%2. soubor na logovani prubehu
[fileID,~] = fopen(logfilename,'wt'); 
%3. tabulka vyslednych CHilbertMulti souboru
tablelog = cell(pocetcyklu+1,6); %frekvence, soubor, reference, status, chyba -  z toho bude vystupni xls tabulka s prehledem vysledku
tablelog(1,:) = {'file','fileno','kategorie','result','file','datetime'}; %hlavicky xls tabulky

for f = 1:numel(files) %cyklus pres vsechny soubory
    for kontrast = 1:numel(kontrasts) %cyklus pres vsechny kontrasty
        stat = setup.stat_kats{kontrast};    
        try
            msg = [' --- ' files{f} ': IntervalyResp *********** '];
            disp(msg); fprintf(fileID,[ msg '\n']);

            CB = CBrainPlot;     %brainplot na ziskani signif odpovedi
            CB.IntervalyResp(testname,[0.1 1],files{f},kontrast); %ziskam signif rozdily pro kategorie a mezi kategoriemi pro vsechny pacienty       
            kategorie = find(~cellfun('isempty',strfind(CB.katstr,'X'))); %strfind je jenom case sensitivni
            for kat = 1:numel(kategorie)
                katstr = CB.katstr{kategorie(kat)}; %jmeno kombinace kategorii z CB, naprikad znackaXvy
                try
                    outfilename = ['d:\eeg\motol\pacienti\0sumarne\CM ' katstr ' ' files{f}]; %jmeno souboru CHilbertMulti
                    CM = CHilbertMulti;
                    if exist(outfilename,'file')==2 && overwriteCM == 0
                        msg = [ ' --- ' strrep(outfilename,'\','\\') ' NEULOZENO, preskoceno '  datestr(now)];
                        disp(msg); fprintf(fileID,[ msg '\n']);                     
                        tablelog(cyklus+1,:) = { files{f}, num2str(f), katstr, 'preskoceno', outfilename,datestr(now) };                     
                        CM.Load(outfilename);
                    else
                        msg = [ ' --- ' strrep(outfilename,'\','\\') ' zpracovavam '  datestr(now)];
                        disp(msg); fprintf(fileID,[ msg '\n']);  

                        %vytvorim extrakty podle tabulky PAC, pro vsechny pacienty a pro tuto kategorii
                        filenames_extract = CM.ExtractData(CB.PAC{1,kategorie(kat)},'aedist',files{f},katstr,overwrite_extracts);

                        FFFilenames_XLS(pocetextracts:pocetextracts+numel(filenames_extract)-1,:) = ...
                            cat(2,repmat({files{f},f,katstr,kat},numel(filenames_extract),1),filenames_extract);
                        pocetextracts = pocetextracts + numel(filenames_extract);                

                        CM.ImportExtract(filenames_extract,katstr);
                        CM.ResponseSearch(0.1,stat);  

                        CM.Save(outfilename);

                        msg = [ ' --- ' files{f} ': ' katstr ' OK '  datestr(now)];
                        disp(msg); fprintf(fileID,[ msg '\n']);            
                        tablelog(cyklus+1,:) = { files{f}, num2str(f), katstr, 'saved', outfilename,datestr(now) }; 
                    end
                    CBo = CBrainPlot; %brainplot na generovani obrazku mozku
                    BPD = CM.ExtractBrainPlotData(); %vytvori data pro import do CBrainPlot
                    CBo.ImportData(BPD); %naimportuje data z CHilbertMulti
                    CBo.PlotBrain3D([],[],[],overwrite_brainplots); %vykresli obrazek mozku
                catch exception 
                    errorMessage = sprintf('** Error in function %s() at line %d.\nError Message:\n%s', ...
                        exception.stack(1).name, exception.stack(1).line, exception.message);                            
                    disp(errorMessage);  fprintf(fileID,[errorMessage '\n']);  %#ok<DSPS> %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal                                        
                    tablelog(cyklus+1,:) = { files{f}, num2str(f), katstr, 'error', exception.message , datestr(now)}; 
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
end

%system('shutdown -h') 