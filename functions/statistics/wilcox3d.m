function wp = wilcox3d(A, B, varargin)
    % A: 3d matrix of interest. The repeated measures need to be in the
    %   3rd axix
    % B: 3d matrix. Can have either of the first two dimensions = 1.
    %   Then it is compared to all values in A
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
    wp = NaN(size(A,1), size(A,2));
    for j = 1:size(A, 1)
        for k = 1:size(A, 2)
            aa = squeeze(A(j,k,:));
            bb = squeeze(B(min(j,size(B,1)), min(k,size(B,2)),:));
            % need to have at least two values
            if any([numel(aa) < 2, numel(bb) < 2])
                wp(j,k) = NaN;
            else
                if p.Results.paired
                    wp(j,k) = signrank(aa,bb);
                else
                    wp(j,k) = ranksum(aa,bb);
                end
            end
        end
    end
    if p.Results.fdr
        [~, ~, adj_p] = fdr_bh(wp, 0.05, p.Results.fdrMethod, 'no');
        wp = adj_p;
    end
end