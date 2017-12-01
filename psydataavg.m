function [ RP,RT,RPe,RTe ] = psydataavg( pacients, testname,filename )
%PSYDATAAVG vrati prumery a stderr odpovedi a reakcnich casu pro vsechny pacienty
%   napr [RP,RT,RPe,RTe] = psydataavg({'p073','p079'},'aedist','AEdist CHilbert 50-120 refBipo Ep2017-11_CHilb.mat');
    RP = []; RPe = [];
    RT = []; RTe = [];
    for p = 1:numel(pacients)
        E = pacient_load(pacients(p),testname,filename);
        [resp,rt,kat,test] = E.PsyData.GetResponses;
        kats = unique(kat);
        if p==1
            RP = zeros(numel(pacients),numel(kats));
            RT = zeros(numel(pacients),numel(kats));
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


end

