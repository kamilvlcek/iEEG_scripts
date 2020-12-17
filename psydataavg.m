function [ RP,RT,RPe,RTe,names,katnames ] = psydataavg( testname,filename )
%PSYDATAAVG vrati prumery a stderr odpovedi a reakcnich casu pro vsechny pacienty
%   napr [RP,RT,RPe,RTe] = psydataavg({'p073','p079'},'aedist','AEdist CHilbert 50-120 refBipo Ep2017-11_CHilb.mat');
    RP = []; RPe = [];
    RT = []; RTe = [];
    names = {}; %jmena pacientu
    if strcmp(testname,'aedist')
        pacienti = pacienti_aedist(); %nactu celou strukturu pacientu    
    elseif strcmp(testname,'ppa')
        pacienti = pacienti_ppa(); %nactu celou strukturu pacientu    
    elseif strcmp(testname,'menrot')
        pacienti = pacienti_menrot(); %nactu celou strukturu pacientu    
    else
        error('neznamy typ testu');
    end
    firstpacient = true;
    for p = 1:numel(pacienti) % cyklus pacienti
        if pacienti(p).todo 
            E = pacient_load(pacienti(p).folder,testname,filename,[],[],[],0); %nejspis objekt CHilbert, pripadne i jiny; loadall = 0
            if isa(E,'CiEEGData')
                [resp,rt,kat,test] = E.PsyData.GetResponses();
                kats = unique(kat);
                if firstpacient
                    RP = nan(numel(pacienti),numel(kats));
                    RPe = nan(numel(pacienti),numel(kats));
                    RT = nan(numel(pacienti),numel(kats));
                    RTe = nan(numel(pacienti),numel(kats));
                    firstpacient = false;
                    katnames = E.PsyData.CategoryName(kats,[]);
                end
                iTest = test==1;          
                for k = 1:numel(kats)
                    iKat = kat==kats(k);
                    RT(p,k) = mean(rt(iTest & iKat)); %#ok<AGROW>
                    RTe(p,k) = std(rt(iTest & iKat))/sqrt(length(rt(iTest & iKat))); %#ok<AGROW>
                    RP(p,k) = mean(resp(iTest & iKat)); %#ok<AGROW>
                    RPe(p,k) = std(resp(iTest & iKat))/sqrt(length(resp(iTest & iKat))); %#ok<AGROW>
                end
            end
            names{p,1} = pacienti(p).folder; %#ok<AGROW>
        end
    end


end


