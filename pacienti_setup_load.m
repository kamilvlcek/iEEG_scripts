function [ pacienti, setup,frekvence,reference ] = pacienti_setup_load( testname,alignresponse )
%PACIENTI_SETUP_LOAD nacte seznam pacientu a nastaveni pro konkretni test
%   [ pacienti, setup,frekvence,reference  ] = pacienti_setup_load( testname )
%   3.5.2018 kvuli castemu opakovani
    if(~exist('alignresponse','var')) || isempty(alignresponse) , alignresponse = 0; end
    if strcmp(testname,'aedist')
        pacienti = pacienti_aedist(); %nactu celou strukturu pacientu   
        setup = setup_aedist( alignresponse );
        [ frekvence,reference ] = freqref_aedist(); %nactu frekvence a reference k vyhodnoceni         
    elseif strcmp(testname,'ppa')
        pacienti = pacienti_ppa(); %nactu celou strukturu pacientu    
        setup = setup_ppa(alignresponse);
        [ frekvence,reference ] = freqref_ppa(); %nactu frekvence a reference k vyhodnoceni 
    elseif strcmp(testname,'menrot')
        pacienti = pacienti_menrot(); %nactu celou strukturu pacientu    
        setup = setup_menrot(alignresponse);
        [ frekvence,reference ] = freqref_menrot(); %nactu frekvence a reference k vyhodnoceni 
    else
        error(['neznamy typ testu: ' testname]);
    end

end

