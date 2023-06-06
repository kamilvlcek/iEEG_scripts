function [ setup ] = setup_memact( alignresponse )
%SETUP_MEMACT returns scripts setup for the test MemoryActions
%   alignresponse = align epochs relative to the behavioral response time

if(~exist('alignresponse','var')) || isempty(alignresponse) , alignresponse = 0; end
setup = {};
setup.basedir = 'd:\eeg\motol\pacienti\';
if alignresponse % align epochs to the response 
    setup.epochtime =  [-2 0.3 1];  % hranice epochy: [-0.3 0.8] PPA, zarovnani podle odpovedi/podnetu [-1 1]; [-0.2 1.2] AEdist [-0.2 2.0] pro ResampleEpochs
    setup.baseline = [-1.9 -1.6]; %baseline [-1 0.8]; [-0.5 -0.2] Aedist 2017. 2017/11 - zase [-.2 0]
else % align epochs to the stimulus 
    setup.epochtime =  [-0.5 2.0 0];  % only for immediate trials - 2s response time
    setup.baseline = [-.2 0]  ; %baseline [-.2 0], similarly to menrot or aedist   
end
setup.prefix = 'MemAct'; %has to be AlloEgo, PPA, AEdist or MemAct
setup.stat_kats = {[0 1 2 3], ...  % 'immed_same';'immed_diff';'del_same';'del_diff'
        {[0 1],[2 3]}, ... % immed vs del
        {[0 2],[1 3]}      % same vs diff
        };   
setup.stat_opak = {}; % contrasts for repetitions. {[1 2],[4 5]}; %PPA opakovani 12 vs 45
setup.subfolder = 'memact'; %subdirectory, specific to the test, can be empty if no subdirectories are used
setup.alignresponse = alignresponse;
end

