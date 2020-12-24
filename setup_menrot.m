function [ setup ] = setup_menrot( alignresponse )
%SETUP_MENROT returns scripts setup for the test Menrot
%   alignresponse = zarovnani epoch podle odpovedi

if(~exist('alignresponse','var')) || isempty(alignresponse) , alignresponse = 0; end
setup = {};
setup.basedir = 'd:\eeg\motol\pacienti\';
if alignresponse %zarovnani epoch podle odpovedi
    setup.epochtime =  [-1.2 0.3 1];  % hranice epochy: [-0.3 0.8] PPA, zarovnani podle odpovedi/podnetu [-1 1]; [-0.2 1.2] AEdist [-0.2 2.0] pro ResampleEpochs
    setup.baseline = [-1.2 -1.0]; %baseline [-1 0.8]; [-0.5 -0.2] Aedist 2017. 2017/11 - zase [-.2 0]
else %zarovnani epoch podle podnetu (tj. normalne)
    setup.epochtime =  [-1 1.5 0];  % hranice epochy: [-0.3 0.8] PPA, zarovnani podle odpovedi/podnetu [-1 1]; [-0.2 1.2] AEdist [-0.2 2.0] pro ResampleEpochs
    setup.baseline = [-.2 0]  ; %baseline [-1 0.8]; [-0.5 -0.2] Aedist 2017. 2017/11 - zase [-.2 0]    
end
setup.prefix = 'Menrot'; %musi byt bud AlloEgo, PPA, AEdist
setup.stat_kats = {[0 1 2 3], ...
        {[0 1],[2 3]}, ... %vy x znacka
        {[0 2],[1 3]}      %3D vs 2D
        };  %  Menrot 'vy-2D' 'vy-3D' 'znacka-2D' 'znacka-3D'; 
setup.stat_opak = {}; %{[1 2],[4 5]}; %PPA opakovani 12 vs 45
setup.subfolder = 'menrot'; %podadresar, specificky pro test, muze byt prazdne pokud se nepouzivaji podadresare
setup.alignresponse = alignresponse;
end

