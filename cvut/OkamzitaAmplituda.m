function A=OkamzitaAmplituda(x)

%OkamzitaAmplituda vypocita okamzitu amplitudu signalu
%   A = OkamzitaAmplituda(X), kde X je jednorozmerny vektor, vypocita
%   okamzitu amplitudu signalu X. Okamzita amplituda je vypocitana ako
%   absolutna hodnota analitickyeho signalu Xa
%
%   Xa=Xr+i*Xi,
%
%   kde Xi je Hilbertova transformacia signalu X.
%
%   A = OkamzitaAmplituda(X), kde X je viacrozmerne pole, prevedie vypocet 
%   okamzitej amplitudy pre kazdy signal prvej dimenzie pola X.
%   
%
%   Priklad pouzitia c. 1:
%   
%     % vygenerovanie signalov
%     fs=100;              % vzorkovacia frekvencia
%     n=1:2000;
%     x=zeros(length(n),3);
%     % amplitudovo modulovany signal
%     x(:,1)=sin(2*pi*2*n/fs).*blackman(length(n))';  
%     % uzkopasmovy sum
%     x(:,2)=filter(1,poly([0.999*exp(j*0.2),0.999*exp(-j*0.2)]),randn(size(n)));
%     % impulzova odozva Butterworthovho filtru
%     [b,a]=butter(5,[0.05 0.055]);
%     x(:,3)=impz(b,a,length(n));
% 
%     % vypocet okamzitej amplitudy
%     A=OkamzitaAmplituda(x);
% 
%     % vykreslenie vysledkov
%     figure
%     for n=1:size(x,2)
%        subplot(size(x,2),1,n)
%        plot([x(:,n) A(:,n)]);
%     end;
% 
%   Priklad pouzitia c. 2: 
% 
%     % vygenerovanie signalov
%     x=zeros(2000,3,3);
%     for k1=1:size(x,2)
%         for k2=1:size(x,3)
%             [b,a]=butter(5,[0.05 0.05+0.002*k2]*k1);
%             x(:,k1,k2)=impz(b,a,2000);
%         end
%     end;
% 
%     % vypocet okamzitej amplitudy
%     A=OkamzitaAmplituda(x);
% 
%     % vykreslenie vysledkov
%     figure
%     for k1=1:size(x,2)
%         for k2=1:size(x,3)
%             subplot(size(x,2),size(x,3),(k1-1)*size(x,2)+(k2-1)+1);
%             plot([x(:,k1,k2),A(:,k1,k2)]);
%         end
%     end;



% Rev160629

x_size=size(x);
if min(x_size)==1 && length(x_size)==2
    x=x(:); 
else 
    x=reshape(x,x_size(1),[]);
end;

A=abs(hilbert(x,size(x,1)*2));
A=A(1:size(x,1),:);
% A=abs(hilbert(x));

A=reshape(A,x_size);
