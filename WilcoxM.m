function [ W ] = WilcoxM( A,B,fdr )
%WILCOXM vrati matici p z Wilcoxova testu ranksum
%   A a B jsou 3d matice, posledni dimenze jsou epochs neboli opakovani.
%   fdr urcuje, jestli se ma pocitat korekce signifikanci podle fdr
W = zeros(size(A,1),size(A,2));
fprintf('Wilcox Test M - proti matrix: ');
for el = 1:size(A,1)
    fprintf('%d ', el);
    %disp([num2str(el) ' ']);
    for s = 1:size(A,2) %cas, test se dela pro kazdy casovy okamzik zvlast
        aa = squeeze (A(el,s,:));
        bb = squeeze (B(el,s,:));
        W(el,s)= ranksum(aa,bb); %Wilcoxon rank sum test
        %[~,p] = ttest2(aa,bb);        W(el,s)= p;
    end
end
if ~exist('fdr','var') , fdr = 1; end
if fdr
    [~, ~, adj_p]=fdr_bh(W,0.05,'pdep','no'); %dep je striktnejsi nez pdep
    W = adj_p; %prepisu puvodni hodnoty korigovanymi podle FDR
end

fprintf(' .. done\n');
end
 

