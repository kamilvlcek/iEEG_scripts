function [files,zmenseno] = cvutBatch(spojit,ratio,yratio,testrozdil)
%spojit = {}; %then insert list of files to join with the full path
if ~exist('ratio','var'), ratio = 4; end %z 2048 na 512
if ~exist('yratio','var'), yratio = 0.01; end  % pro novy system Quantum
if ~exist('testrozdil','var'), testrozdil = 1; end %defaultne chci testovat casovy rozdil mezi daty pri spojovani.

zmenseno = zmensidata(spojit,ratio,yratio); %vraci se seznam souboru ke spojeni vcetne plne cesty

[files] = concatCVUT(zmenseno,[],testrozdil); %vraci files: {filename, spojeno, delka, rozdil_sec}

k = strfind(spojit{1},'\');
folder = spojit{1}(1:k(end));
logfilename = [folder 'cvutBatch' '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')];
warning('off','MATLAB:xlswrite:AddSheet'); %[msg,msgID] = lastwarn;
xlswrite([logfilename '.xls'],[{'filename','spojeno','delka','rozdil_sec'}; files],'files'); %files v xls tabulce, sheed 'files'
xlswrite([logfilename '.xls'],zmenseno,'zmenseno'); %files v xls tabulce, sheed 'files'


