function [ setup ] = setup_ppa( alignresponse )
%SETUP_AEDIST funkce vraci nastaveni vyhodnoceni testu Aedist
%   Detailed explanation goes here
if(~exist('alignresponse','var')) || isempty(alignresponse) , alignresponse = 0; end
setup = {};
setup.basedir = 'd:\eeg\motol\pacienti\';
if alignresponse %zarovnani epoch podle odpovedi %to u PPA nema smysl
    setup.epochtime =  [-0.8 0.3 1];  % hranice epochy: [-0.3 0.8] PPA, zarovnani podle odpovedi/podnetu [-1 1]; [-0.2 1.2] AEdist [-0.2 2.0] pro ResampleEpochs
    setup.baseline = [-0.8 -0.5]; %baseline [-1 0.8]; [-0.5 -0.2] Aedist 2017. 2017/11 - zase [-.2 0]
else %zarovnani epoch podle podnetu (tj. normalne)
    setup.epochtime =  [-0.2 0.8];  % hranice epochy: [-0.2 0.8] PPA (epochy nekdy cele 1s)
    setup.baseline = [-.05 0]  ; %baseline [-1 0.8]; 
    
end
setup.suffix = 'Ep2017-11'; %Ep
if(alignresponse)
   setup.suffix = [setup.suffix 'Resp']; %pokud zarovnavam podle odpovedi, pridavam priponu
end
setup.prefix = 'PPA'; %musi byt bud AlloEgo, PPA, AEdist
setup.stat_kats = {[2 3 1],[1 3 2],[1 2 3]};  % PPA [2 3 1] 2=Face x 3=Object x 1=Scene ; SceneX, FaceX, ObjectX
setup.stat_opak = {}; %{[1 2],[4 5]}; %PPA opakovani 12 vs 45
setup.subfolder = 'PPA'; %podadresar, specificky pro test, muze byt prazdne pokud se nepouzivaji podadresare
setup.alignresponse = alignresponse;
end

