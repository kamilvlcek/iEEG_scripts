function [files,zmenseno] = cvutBatch(spojit,ratio,yratio,testrozdil)

if ~exist('testrozdil','var'), testrozdil = 1; end %defaultne chci testovat casovy rozdil mezi daty pri spojovani.
%spojit = {};
%ratio = 4; %z 2048 na 512
%yratio = 0.01; % pro novy system Quantum
 

zmenseno = zmensidata(spojit,ratio,yratio); %vraci se seznam souboru ke spojeni vcetne plne cesty

[files] = concatCVUT(zmenseno,[],testrozdil);
