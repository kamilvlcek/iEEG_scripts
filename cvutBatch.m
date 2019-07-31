%% inicializace 
%spojit = {};
ratio = 4; %z 2048 na 512
yratio = 0.01; % pro novy system Quantum
 
%%
%ted si  naplnim spojit adresarema

%% zmensim data
zmenseno = zmensidata(spojit,ratio,yratio);

[files] = concatCVUT(zmenseno);
