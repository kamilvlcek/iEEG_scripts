function p_corrected=CmpPerm(X,Y,Nperm1,Nperm2)

% CmpPerm compares distributions of the arithmetic means using permutation tests.
%   p_corrected=CmpPerm(X,Y,Nperm1,Nperm2) performs comparison of the
%   distribution of the arithmetic means of multiple signals in the arrays 
%   X and Y, and computes corrected p values controlling the family-wise error.
%
%   The null hypothesis H0 is that the distributions of the arithmetic means
%   of corresponding signals in X and Y are identical. First, permutations 
%   are used to create empirical distributions of the difference of the 
%   arithmetic means under H0. Then, these empirical distributions are used to 
%   compute the uncorrected p values. The corrected p values are then computed 
%   using a permutation algorithm with a step down correction.
%
%   Arrays X and Y can be N dimensional, where the realizations are always 
%   in the last dimension.
%
%   Nperm1 is the number of permutations used for the computation of the
%   empirical distribution under H0.
%
%   Nperm2 is the number of permutations used to compute the corrected p
%   values.
%
%
%   Example 1:
%
%   % generate signals under H0
%   X=randn(10,12,40); 
%   Y=randn(10,12,50);
%   % compute corrected p values
%   p_corrected=CmpPerm(X,Y,30000,2000);
%   
%   Example 2:
%
%   % generate signals where one signal violates H0
%   X=randn(10,12,40); 
%   Y=randn(10,12,50);
%   X(2,3,:)=X(2,3,:)+1; % force different mean
%   % compute corrected p values
%   p_corrected=CmpPerm(X,Y,30000,2000);
%
if isempty(X) || isempty(Y)
  p_corrected = []; %kamil - jinak nastane chyba  
  return
end
Xsize=size(X);
Ysize=size(Y);
X=reshape(X,[],Xsize(end));
Y=reshape(Y,[],Ysize(end));

Nsig=size(X,1);
Lsig1=size(X,2);
Lsig2=size(Y,2);

X=X'; Y=Y';

mx=mean(X);
my=mean(Y);
t0=abs(mx-my); %./(std(X).*std(Y));

Z=[X; Y];

Fempirical=empirical_distribution(Z,Lsig1,Lsig2,Nperm1);
p_uncorr=(sum(Fempirical>repmat(t0,[Nperm1 1]))+1)/(Nperm1+1);
p_corr=permcorrect(p_uncorr,Z,Lsig1,Lsig2,Fempirical,Nperm2);
p_corrected=reshape(p_corr,[Xsize(1:end-1) 1]);
