function [ setup ] = setup_memact(typeEpochs)
%SETUP_MEMACT returns scripts setup for the test MemoryActions
%   typeEpochs = to use one of four different epochs with different epoch time and baseline and alignment to stimulus:
% 0 - immediate epochs
% 1 - epochs before delay
% 2 - epochs after delay
% 3 - epochs within delay
% 4 - baseline activity before encoding phase, requires for normalization and computing stats for epochs within delay

if(~exist('typeEpochs','var')) || isempty(typeEpochs), typeEpochs = 0; end % by default analyses immediate epochs
setup = {};
setup.basedir = 'd:\eeg\motol\pacienti\';

if typeEpochs == 0  % immediate epochs 

    % The setup.epochtime(3) has 3 possible values
    % 0 - align to stimuli 
    % 1 - align to responses 
    % 2 - align to stimuli - delay 
    setup.epochtime =  [-0.5 2.0 0];  % immediate trials - 2s response time
    setup.baseline = [-.2 0]; % baseline [-.2 0], similarly to menrot or aedist
    setup.index = 1; % serves as index to rjepoch struct 
    setup.suffix = 'imm'; % short name of epoch type to distinguish the CHilbert file created; used in BatchHilbert
    setup.filter = {7,[0 1]}; % used in CiEEGData.ExtractEpochs: filter{1} - column number in obj.P.data, filter{2} - searches for kats values in this column
    % 0 - 'immed_same'; 1 - 'immed_diff'
    setup.stat_kats = {[0 1],[1 0]}; % immed_same x immed_diff
    
elseif typeEpochs == 1 % epochs before delay 
    setup.epochtime =  [-0.5 3.95 2]; % encoding phase (2 sec) + first part of delay (1.95 sec); % 2 - align to stimuli - delay (new)- stimulus which is presented during the encoding phase
    setup.baseline = [-0.5 -0.1]; % isn't used in BatchHilbert; use it only after appending 2 objects in function CM.NormalizeEpochs([-0.5 -0.1]);
    setup.index = 2;
    setup.suffix = 'bdel';
    setup.filter = {7,[2 3]}; % 2 - 'del_same'; 3 - 'del_diff'
    setup.stat_kats = {[2 3],[3 2]}; % del_same x del_diff
    
elseif typeEpochs == 2 % epochs after delay
    setup.epochtime =  [-1.95 2 0];  % second part of delay (1.95 sec) + action phase (2 sec); stimulus - start of action phase
    setup.baseline = [-0.5 -0.1]; % baseline activity before encoding phase (ITI); 
        % isn't used in BatchHilbert; use it only after appending 2 objects in function CM.NormalizeEpochs([-0.5 -0.1]);
    setup.index = 3;
    setup.suffix = 'adel';
    setup.filter = {7,[2 3]}; % 2 - 'del_same'; 3 - 'del_diff'
    setup.stat_kats = {[2 3],[3 2]}; % del_same x del_diff
    
elseif typeEpochs == 3 % epochs within delay   
    % setup.epochtime =  [-3.9 0.3 0]; % stimulus - start of action phase; how to normalize such epochs?
    setup.epochtime =  [2.8 5.3 2]; % stimulus - start of encoding phase, in this way we can exctract the whole delay and normalize by baseline acvitity before encod phase
%     setup.epochtime =  [2.0 5.9 2];
    setup.baseline = [-0.5 -0.1]; % baseline activity before encoding phase 
        % isn't used in BatchHilbert; use it only after appending 2 objects in function CM.NormalizeEpochs([-0.5 -0.1]);
    setup.index = 4;
    setup.suffix = 'del';
    setup.filter = {7,[2 3]}; % 2 - 'del_same'; 3 - 'del_diff'
    setup.stat_kats = {[2 3],[3 2]}; % del_same x del_diff
    
elseif typeEpochs == 4 % baseline activity before encoding phase, used for normalization and computing stats for epochs within delay
    setup.epochtime =  [-0.5 0 2];
    setup.baseline = [-0.5 -0.1]; % isn't used in BatchHilbert; use it only after appending with epochs within delay in function CM.NormalizeEpochs([-0.5 -0.1]);
    setup.index = 4; % the same as for within delay epochs
    setup.suffix = 'bs';
    setup.filter = {7,[2 3]}; % 2 - 'del_same'; 3 - 'del_diff'
    setup.stat_kats = {[2 3],[3 2]}; % del_same x del_diff
else
    disp('wrong value of typeEpochs, use only values: 0, 1, 2, 3 or 4')
end
setup.prefix = 'MemAct'; %has to be AlloEgo, PPA, AEdist or MemAct
% setup.stat_kats = {[0 1 2 3], ...  % 'immed_same';'immed_diff';'del_same';'del_diff'
%         {[0 1],[2 3]}, ... % immed vs del
%         {[0 2],[1 3]}      % same vs diff
%         };   
setup.stat_opak = {}; % contrasts for repetitions. {[1 2],[4 5]}; %PPA opakovani 12 vs 45
setup.subfolder = 'memact'; %subdirectory, specific to the test, can be empty if no subdirectories are used
setup.typeEpochs = typeEpochs;
end

