function [ ] = BatchExtracts( testname,files,kontrasts )

% files = {   'AEdist CHilbert 50-150 -0.5-1.2 refBipo Ep2018-04_CHilb.mat',...
%             'AEdist CMorlet 1-10M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat',...
%             'AEdist CMorlet 4-8M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat',...
%             'AEdist CMorlet 1-4M -0.5-1.2 refBipo Ep2018-06 FE_CHilb.mat'};       
% 
% testname = 'aedist';

[ pacienti, setup  ] = pacienti_setup_load( testname );
if ~exist('kontrasts','var'), kontrasts = 1:numel(setup.stat_kats); end %statisticky kontrast, pokud nezadam zvnejsku, udelam vsechny

pocetcyklu = 0;
for kontrast = 1:numel(kontrasts) %cyklus jen na vypocet celkoveho poctu cyklu pres vsechny kontrasty ve statistice
    stat = setup.stat_kats{kontrast}; %resp setup.stat_kats{1} {2} nebo {3} pro menrot 
    kombinace_kat = combinator(length(stat),2,'p'); %permutace bez opakovani z poctu kategorii
    kombinace_kat = kombinace_kat(kombinace_kat(:,1)>kombinace_kat(:,2),:); %vyberu jen permutace, kde prvni cislo je vetsi nez druhe     
    pocetcyklu = pocetcyklu + numel(files) * size(kombinace_kat,1) ;
end

overwrite_extracts = 0; %jestli se maji prepisovat extrakty pro kazdeho pacienta
overwrite_brainplots = 1;
overwriteCM = 0; %jestli se maji prepisovat soubory CHilbertMulti
doIntervalyResp = 0; %jestli se maji hledaty signif soubory pres vsechny pacienty pomoci CN.IntervalyResp, pokud ne, potrebuju uz mit hotove CHilbertMulti soubory
cyklus = 1;
pocetextracts = 1;
%dirCM = 'd:\eeg\motol\CHilbertMulti\Menrot\'; %musi koncit \
%fileCS = 'd:\eeg\motol\CHilbertMulti\Menrot\CSelCh_Menrot.mat';
dirCM = 'd:\eeg\motol\CHilbertMulti\Aedist\'; %musi koncit \
fileCS = 'd:\eeg\motol\CHilbertMulti\Aedist\CSelCh_AEdist.mat';
brainplots_onlyselch = 1; %generovat CBrainPlot3D jedine ze souboru, kde jsou selected channels
plotallchns = 0; %jestli generovat obrazky mozku i se vsema kanalama (bez ohledu na signifikanci)

%LOG SOUBORY
%1. seznam vsech extraktu
logfilename = ['logs\BatchExtract_' testname '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.log'];
FFFilenames_logname = ['logs\BatchExtractFilenames_' testname '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.xls'];
FFFilenames_XLS = cell(1+pocetcyklu*numel(pacienti),5); 
FFFilenames_XLS(1,:)={'file','fileno','kat','katno','extract'};
%2. soubor na logovani prubehu
[fileID,~] = fopen(logfilename,'wt'); 
%3. tabulka vyslednych CHilbertMulti souboru
tablelog = cell(pocetcyklu+1,7); 
tablelog(1,:) = {'file','fileno','kategorie','stat','result','file','datetime'}; %hlavicky xls tabulky
if exist('fileCS','var') && exist(fileCS,'file')==2
    CS = CSelCh(fileCS);
end

for f = 1:numel(files) %cyklus pres vsechny soubory
    for kontrast = 1:numel(kontrasts) %cyklus pres vsechny kontrasty
        stat = setup.stat_kats{kontrasts(kontrast)};      %#ok<NASGU>
        try
            if doIntervalyResp
                msg = [' --- ' files{f} ': IntervalyResp *********** ']; %#ok<UNRCH>
                disp(msg); fprintf(fileID,[ msg '\n']);

                CB = CBrainPlot;     %#ok<USENS> %brainplot na ziskani signif odpovedi
                CB.IntervalyResp(testname,[0.1 1],files{f},kontrasts(kontrast)); %ziskam signif rozdily pro kategorie a mezi kategoriemi pro vsechny pacienty       
                kategorie = find(~cellfun('isempty',strfind(CB.katstr,'X'))); %strfind je jenom case sensitivni
                katsnames = CB.katstr;
            else
                msg = [' --- ' files{f} ': Load *********** '];
                disp(msg); fprintf(fileID,[ msg '\n']);
                E = pacient_load(pacienti(2).folder,testname,files{f},[],[],[],0); %nejspis objekt CHilbert, pripadne i jiny; loadall = 0
                E.SetStatActive(kontrasts(kontrast));
                katsnames = E.GetKatsNames();
                kategorie = find(~cellfun('isempty',strfind(katsnames,'X'))); %strfind je jenom case sensitivni
            end
            for kat = 1:numel(kategorie)
                katstr = katsnames{kategorie(kat)}; %jmeno kombinace kategorii z CB, naprikad znackaXvy
                try
                    outfilename = [dirCM 'CM ' katstr ' ' files{f}]; %jmeno souboru CHilbertMulti
                    CM = CHilbertMulti;
                    if exist(outfilename,'file')==2 && overwriteCM == 0
                        msg = [ ' --- ' strrep(outfilename,'\','\\') ' nacteno '  datestr(now)];
                        disp(msg); fprintf(fileID,[ msg '\n']);                     
                        tablelog(cyklus+1,:) = { files{f}, num2str(f), katstr,cell2str(stat), 'nacteno', outfilename,datestr(now) };                     
                        CM.Load(outfilename);                        
                    elseif doIntervalyResp
                        msg = [ ' --- ' strrep(outfilename,'\','\\') ' zpracovavam '  datestr(now)]; %#ok<UNRCH>
                        disp(msg); fprintf(fileID,[ msg '\n']);  

                        %vytvorim extrakty podle tabulky PAC, pro vsechny pacienty a pro tuto kategorii
                        filenames_extract = CM.ExtractData(CB.PAC{1,kategorie(kat)},testname,files{f},katstr,overwrite_extracts);

                        FFFilenames_XLS(pocetextracts:pocetextracts+numel(filenames_extract)-1,:) = ...
                            cat(2,repmat({files{f},f,katstr,kat},numel(filenames_extract),1),filenames_extract);
                        pocetextracts = pocetextracts + numel(filenames_extract);                
                        
                        %FILES = CM.TestExtract(filenames_extract);
                        CM.ImportExtract(filenames_extract,katstr);
                        CM.ResponseSearch(0.1,stat);  

                        CM.Save(outfilename);

                        msg = [ ' --- ' files{f} ': ' katstr ' OK '  datestr(now)];
                        disp(msg); fprintf(fileID,[ msg '\n']);            
                        tablelog(cyklus+1,:) = { files{f}, num2str(f), katstr,cell2str(stat), 'saved', outfilename,datestr(now) }; 
                    else
                        msg = [ ' --- ' strrep(outfilename,'\','\\') ' nevytvoreno, nenacteno '  datestr(now)];
                        disp(msg); fprintf(fileID,[ msg '\n']);   
                        tablelog(cyklus+1,:) = { files{f}, num2str(f), katstr,cell2str(stat), 'nothing to do', outfilename,datestr(now) }; 
                    end
                    if exist('CS','var')                        
                        selCh = CS.GetSelCh(CM); %pokud se tam nazev souboru najde, vlozi se selected channels, jinak ne
                    else
                        selCh = [];
                    end
                    if ~brainplots_onlyselch || ~isempty(selCh) %pokud negenerovat jen pro selch, nebo pokud nejsou prazne selch
                        CBo = CBrainPlot; %brainplot na generovani obrazku mozku
                        BPD = CM.ExtractBrainPlotData([],kategorie(kat)); %vytvori data pro import do CBrainPlot
                        CBo.ImportData(BPD); %naimportuje data z CHilbertMulti
                        CBo.PlotBrain3D(iff(plotallchns,[1 2],2),[],[],overwrite_brainplots); %vykresli obrazek mozku
                    end
                catch exception 
                    errorMessage = exceptionLog(exception);                            
                    disp(errorMessage);  fprintf(fileID,[errorMessage '\n']);   %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal                                        
                    tablelog(cyklus+1,:) = { files{f}, num2str(f), katstr, 'error', exception.message , datestr(now)}; 
                    clear CM; 
                end 
                cyklus = cyklus + 1;
                xlswrite([logfilename '.xls'],tablelog); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
                xlswrite(FFFilenames_logname,FFFilenames_XLS); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
            end
        catch exception
            errorMessage = exceptionLog(exception);                          
            disp(errorMessage);  fprintf(fileID,[errorMessage '\n']);   %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal                                        
            tablelog(cyklus+1,:) = { files{f}, num2str(f), 'no kat', 'error', exception.message , datestr(now)}; 
            clear CB;        
            cyklus = cyklus + 1;
            xlswrite([logfilename '.xls'],tablelog); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
            xlswrite(FFFilenames_logname,FFFilenames_XLS); %budu to psat znova po kazdem souboru, abych o log neprisel, pokud se program zhrouti
        end
    end
end

%system('shutdown -h') 