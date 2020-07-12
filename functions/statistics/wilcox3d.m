% TODO - this behaves weridly for the unequal sizes in the first axis
% Say if the matrices are time x channel x frequency and you are testing
% agains the baseline. Therefore matrix A is 50 x 10 x 10 and B is say 
% 10 x 10 x 10, as the baseline is shorter. A(1:10, :, :) will be compared
% to B(1:10,:,:) in B, but A(11:50, :, :) will ALSO be comapred to B(10, :, :)
% Basically all the above B size are compared to the last B

function wp = wilcox3d(A, B, varargin)
% WILCOX3D compares two 3d matrices in their 3rd dimension.
% Original author: Lukáš Hejtmánek
%
% A: 3d matrix of interest. The repeated measures need to be in the
%   3rd axix
% B: 3d matrix. Needs to either have the same size in  Can have size either of the first two dimensions == 1.
%   Then it is compared to all values in A. 
% fdr: logical if fdr should be calculated
% fdrMethod: string defining 'pdep' (default) or 'dep'
% paired: logical, if paired test to be used. Default 0
is3dmatrix = @(x)length(size(x)) == 3 && isnumeric(x);
p = inputParser;
addRequired(p, 'A', is3dmatrix);
addRequired(p, 'B', is3dmatrix);
addParameter(p, 'paired', false, @islogical);
addParameter(p, 'fdr', false, @islogical);
addParameter(p, 'fdrMethod', 'pdep',...
    @(x)isstring(x) && any(strcmp({'pdep' 'dep'},x)));
parse(p, A, B, varargin{:});
A = p.Results.A; B = p.Results.B;

%% validations
% B needs to either have size 1 or it needs to be the same as A in the first
% two dimensions
sizeValid = (size(B,1) == size(A,1)) || (size(B,1) == 1);
sizeValid = sizeValid && ((size(B,2) == size(A,2)) || (size(B,2) == 1));
if ~sizeValid, warning('Passed matrices do not have acceptable sizes'); end
% need to have at least two values
observationsValid = (size(A,3) > 1 || size(B,3) > 1);
if ~observationsValid, warning('Passed matrices do not have acceptable sizes'); end
if(any(~[sizeValid observationsValid]))
    wp = [];
    return
end

%% calculation
wp = NaN(size(A,1), size(A,2));

for jA = 1:size(A, 1)
    jB = min(size(B,1), jA); %either jA, or 1
    for kA = 1:size(A, 2)
        kB = min(size(B,2), kA);%either kA, or 1
        aa = squeeze(A(jA, kB,:));
        bb = squeeze(B(jB, kB, :));
        if p.Results.paired
            wp(jA,kA) = signrank(aa,bb);
        else
            wp(jA,kA) = ranksum(aa,bb);
        end
    end
end
if p.Results.fdr
    [~, ~, adj_p] = fdr_bh(wp, 0.05, p.Results.fdrMethod, 'no');
    wp = adj_p;
end
end