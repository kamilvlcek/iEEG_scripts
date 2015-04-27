function [ W ] = WilcoxA( A,C )
%WILCOXA vrati matici p z Wilcoxova testu ranksum
%   A a C jsou 3d matice, C ma ale jen jednu hodnotu v druhe dimenzi
%   rozmery A: elektrody, cas, epochy=opakovani

%W = rand(size(A,1),size(A,2));
W = zeros(size(A,1),size(A,2));
fprintf('Wilcox Test A - proti Array:');
for el = 1:size(A,1)
    fprintf('%d ', el);
    %disp([num2str(el) ' ']);
    for s = 1:size(A,2)
        aa = squeeze (A(el,s,:));
        bb = squeeze (C(el,1,:));
        W(el,s)= ranksum(aa,bb);
    end
end
[~, ~, adj_p]=fdr_bh(W,0.05,'dep','no'); %dep je striktnejsi nez pdep
W = adj_p; %prepisu puvodni hodnoty korigovanymi podle FDR

fprintf(' done\n ');
end
