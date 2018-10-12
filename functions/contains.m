function tf = contains(s, pattern)
%CONTAINS True if pattern is found in text.
%   TF = CONTAINS(S,PATTERN) returns true if CONTAINS finds PATTERN in any
%   element of string array S. TF is the same size as S.
%
%   S can be a string array, a character vector, or a cell array of
%   character vectors. So can PATTERN. PATTERN and S need not be the same
%   size. If PATTERN is nonscalar, CONTAINS returns true if it finds any
%   element of PATTERN in S.
% 
%   TF = CONTAINS(S,PATTERN,'IgnoreCase',IGNORE) ignores case when searching 
%   for PATTERN in S if IGNORE is true. The default value of IGNORE is false.
% 
%   Examples
%       S = string('data.tar.gz');
%       P = string('tar');
%       contains(S,P)                   returns  1
%
%       S = string({'abstracts.docx','data.tar.gz'});
%       P = 'tar';         
%       contains(S,P)                   returns  [0 1]
%
%       S = string('data.tar.gz');
%       P = {'docx','tar'};
%       contains(S,P)                   returns  1
%
%       S = string({'DATA.TAR.GZ','SUMMARY.PPT'});
%       P = string('tar');
%       contains(S,P,'IgnoreCase',true) returns  [1 0]
%
%   See also endsWith, startsWith.

%   Copyright 2015-2016 The MathWorks, Inc.
%   Kamil 2018-10-05 kvuli matlab 2016a

    
    assert(ischar(s) || iscellstr(s),'prvni argument musi byt char nebo cellstr');
    index = false(size(s));
    if ~iscell(pattern), pattern = {pattern}; end
    for jj = 1:numel(pattern) %protoze contains umi hledat vic patternu najednou a ja to pouzivam
        indexjj =  find(~cellfun('isempty',strfind(lower(s),lower(pattern{jj})))); %nahrazuju contains pomoci strfind
        if numel(indexjj)>0
            index(indexjj) = true;
        end
    end
    tf = index; %vraci logical index - rozmer jako s
end