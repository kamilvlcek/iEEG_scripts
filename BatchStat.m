function [ ] = BatchStat( testname,filenames )
%BATCHSTAT vypocita statistiku pro vsechny soubory pro jeden test
%   puvodne hlavne kvuli permutacni statistice a potom kvuli doplneni statistiky o dalsi 
if ~iscell(filenames), filenames = {filenames}; end %kdyz zadam jen jedno jmeno souboru
cfg=struct;
cfg.srovnejresp = 0;
cfg.hybernovat = 0;
%cfg.statmethod = 'permut'; 
cfg.statmethod = 'wilcox';

if ~isfield(cfg,'hybernovat'), cfg.hybernovat = 0; end %jestli chci po konci skriptu pocitac uspat - ma prednost
if ~isfield(cfg,'vypnout'), cfg.vypnout = 0; end %jestli chci po konci skriptu pocitac vypnout (a nechci ho hybernovat) 
if ~isfield(cfg,'srovnejresp'), cfg.srovnejresp = 0; end %jestli se maji epochy zarovnava podle odpovedi

if strcmp(testname,'aedist')
    pacienti = pacienti_aedist(); %nactu celou strukturu pacientu    
    setup = setup_aedist( cfg.srovnejresp );
elseif strcmp(testname,'ppa')
    pacienti = pacienti_ppa(); %nactu celou strukturu pacientu    
    setup = setup_ppa( cfg.srovnejresp );
elseif strcmp(testname,'menrot')
    pacienti = pacienti_menrot(); %nactu celou strukturu pacientu    
    setup = setup_menrot( cfg.srovnejresp ); 
else
    error('neznamy typ testu');
end

logfilename = ['logs\BatchStat_' setup.prefix '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.log'];
[fileID,~] = fopen(logfilename,'wt'); %soubor na logovani prubehu
assert(fileID>=0,['nemohu otevrit soubor pro zapis: ' logfilename ]);
setuptext = setup2text(setup,cfg);
fprintf(fileID,setuptext); %ulozi setup do log souboru

stat_kats = setup.stat_kats;
stat_opak = setup.stat_opak;
batchtimer = tic; %zacnu merit celkovy cas
for f = 1:numel(filenames) %muzu zpracovavat vic souboru za sebou
    filename = filenames{f};
    msg = [' +++++ ' filename ' +++++ '];
    disp(msg); fprintf(fileID,[msg '\n']);
    for p = 1:numel(pacienti) % cyklus pacienti
        if pacienti(p).todo
            msg = ['***   ' pacienti(p).folder '   ***'];
            disp(msg); fprintf(fileID,[msg '\n']);
            try
                E = pacient_load(pacienti(p).folder,testname,filename); %nejspis objekt CHilbert, pripadne i jiny
                if isempty(E)
                    msg = 'no data';
                    disp(msg); fprintf(fileID,[msg '\n']);
                    pacienti(p).todo = 0; %nechci ho dal zpracovavat
                    continue;
                end
                if iscelldeep(stat_kats) %pokud mam nekolik ruznych statistik na spocitani
                    for WpA = 1:numel(stat_kats)
                        E.SetStatActive(WpA);
                        disp(['pocitam kontrast' cell2str(stat_kats{WpA}) ]);
                        E.ResponseSearch(0.1,stat_kats{WpA},stat_opak,cfg.statmethod);
                    end
                else
                    E.ResponseSearch(0.1,stat_kats, stat_opak,cfg.statmethod);
                end
                E.Save();
                msg = ' file saved OK';
                disp(msg); fprintf(fileID,[msg '\n']);
                cas = toc(batchtimer); %zjistim celkovy cas v sec
                msg = sprintf(' cas zatim: %.1f min',cas/60); %celkovy cas v minutach
                disp(msg); fprintf(fileID,[msg '\n']); %#ok<DSPS>
            catch exception 
                errorMessage = sprintf('** Error in function %s() at line %d.\nError Message:\n%s', ...
                    exception.stack(1).name, exception.stack(1).line, exception.message);                            
                disp(errorMessage);  fprintf(fileID,[errorMessage '\n']);  %#ok<DSPS> %zobrazim hlasku, zaloguju, ale snad to bude pokracovat dal                            
                clear E ans; 
            end  
        end
    end
end

fclose(fileID);

if cfg.hybernovat
    system('shutdown -h') 
elseif cfg.vypnout            
    system('shutdown -s')
end

end

