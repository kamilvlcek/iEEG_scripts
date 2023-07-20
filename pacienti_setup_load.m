function [ pacienti, setup,frekvence,reference ] = pacienti_setup_load( testname,typeEpochs )
%PACIENTI_SETUP_LOAD loads a list of patients and settings for a specific test
%   [ pacienti, setup,frekvence,reference  ] = pacienti_setup_load( testname )
%   3.5.2018 
    if(~exist('typeEpochs','var')) || isempty(typeEpochs) , typeEpochs = 0; end
    if strcmp(testname,'aedist')
        pacienti = pacienti_aedist(); %nactu celou strukturu pacientu   
        setup = setup_aedist( typeEpochs );
        [ frekvence,reference ] = freqref_aedist(); %nactu frekvence a reference k vyhodnoceni         
    elseif strcmp(testname,'ppa')
        pacienti = pacienti_ppa(); %nactu celou strukturu pacientu    
        setup = setup_ppa(typeEpochs);
        [ frekvence,reference ] = freqref_ppa(); %nactu frekvence a reference k vyhodnoceni 
    elseif strcmp(testname,'menrot')
        pacienti = pacienti_menrot(); %nactu celou strukturu pacientu    
        setup = setup_menrot(typeEpochs);
        [ frekvence,reference ] = freqref_menrot(); %nactu frekvence a reference k vyhodnoceni 
    elseif strcmp(testname,'memact')
        pacienti = pacienti_memact(); %load the entire structure of the patients   
        setup = setup_memact(typeEpochs); %load the setup for this test
        [ frekvence,reference ] = freqref_memact(); %load frequencies and references
    else
        error(['unknown test type: ' testname]);
    end

end

