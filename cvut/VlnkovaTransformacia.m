function [S,f,t]=VlnkovaTransformacia(x,scales,fs)

%VlnkovaTransformacia vypocita aproximaciu spojitej vlnkovej transformacie
%
%   S = VlnkovaTransformacia(X,SKALY,Fs), kde X je jednorozmerny vektor,
%   vrati aproximaciu spojitej vlnkovej transformacie signalu X.
%   Vlnkova transformacia je vypocitana pomocou komplexneho Morletovho
%   okna 
%     
%   exp(-t^2/2)*(exp(j*5*t)-exp(-(5^2)/2))/(pi^0.25)
%
%   skalovaneho pomocou skal vo vektore SKALY. Fs urcuje vzorkovaciu
%   frekvenciu signalu. Vysledkom je matica S, ktorej riadky, resp. stlpce
%   obsahuju hodnoty spojitej vlnkovej transformacie pre rozne skaly,
%   resp. casy. 
%
%   S = VlnkovaTransformacia(X,SKALY,Fs), kde X je viacrozmerne pole,
%   prevedie vypocet vlnkovej transformacie pre kazdy signal prvej
%   dimenzie X.
%
%   [S,f,t] = VlnkovaTransformacia(...) vrati vektor pseudofrekvencii f
%   a vektor casov t, pre ktore bola Vlnkova transformacia napocitana.
%
%   Priklad pouzitia c. 1:
%   
%     % vygenerovanie signalov
%     fs=100;              % vzorkovacia frekvencia
%     n=1:2000;
%     x=zeros(length(n),3);
%     % signal obsahujuci harmonicku zlozku a dirakov impulz
%     x(:,1)=sin(2*pi*20*n/fs);  x(500,1)=10;  
%     % chirp signal
%     x(:,2)=sin(2*pi*n.^2/(fs.^2)); % chirp signal
%     % frekvencne modulovany signal
%     xm=10*sin(2*pi*0.2*n/fs); 
%     ph=filter(1,[1 -1],2*pi*(xm+25)/fs);
%     x(:,3)=sin(ph);
% 
%     % vypocet vlnkovej transformacie
%     scales=logspace(-0.1,-1.8,100);
%     [S,f,t]=VlnkovaTransformacia(x,scales,fs);
% 
%     % vykreslenie vysledkov
%     figure
%     for k=1:size(x,2)
%        subplot(size(x,2),1,k)
%        pcolor(t,f,abs(S(:,:,k)))
%        set(gca,'yscale','log')
%        shading flat
%     end;
% 
%     Priklad pouzitia c.2: 
% 
%     % vygenerovanie signalov
%     fs=100;              % vzorkovacia frekvencia
%     n=1:2000;
%     x=zeros(length(n),3,3);
%     for k1=1:size(x,2)
%         for k2=1:size(x,3)
%             x(:,k1,k2)=sin(2*pi*k1*k2*2*n/fs);
%             x(k1*k2*50,k1,k2)=10;
%         end
%     end;
% 
%     % vypocet vlnkovej transformacie
%     scales=logspace(-0.1,-1.8,100);
%     [S,f,t]=VlnkovaTransformacia(x,scales,fs);
% 
%     % vykreslenie vysledkov
%     figure
%     for k1=1:size(x,2)
%         for k2=1:size(x,3)
%             subplot(size(x,2),size(x,3),(k1-1)*size(x,2)+(k2-1)+1);
%             pcolor(t,f,abs(S(:,:,k1,k2)));
%             set(gca,'yscale','log')
%             shading flat
%         end
%     end;


% Rev160629

x_size=size(x);
if min(x_size)==1 && length(x_size)==2
    x=x(:); 
else 
    x=reshape(x,x_size(1),[]);
end;


N=size(x,1); %delka vektoru = cas
Nsig=size(x,2); %pocet kanalu

S=zeros(length(scales),N,Nsig); %rozmery frekvence x delka x kanaly

fprintf('kanal ze %i: ', Nsig );
for n=1:Nsig    
    for k=1:length(scales)
        a=scales(k); %sirka waveletu 
        t=(-5*a:1/fs:a*5); %cas - delka waveletu v sec
        w=conj(  1/sqrt(a*sqrt(pi)) *   exp(-(t/a).^2/2) .* (exp(1i*5*(t/a))-exp(-(5^2)/2)));
        % complex conjugate (  (A=frequency band-specific scaling factor) * gaussian * complex sin  )
        y=conv(x(:,n),w,'same')/fs; %convolution - kvuli tomu conj predtim je vysledek rovnou power?
        S(k,:,n)=y;
    end
    fprintf('%i,',n);
end

if min(x_size)~=1 || length(x_size)~=2
    S=reshape(S,[size(S,1), size(S,2), x_size(2:end)]);
end;

f=5./(scales*2*pi);
t=(0:N-1)/fs;
