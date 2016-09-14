function [S,f,t]=STFT(x,w,Noverlap,fs)

%STFT Kratkodoba Fourierova transformacia (Short Time Fourier Transform - STFT)
%
%   S = STFT(X,OKNO,PREKRYV), kde X je jednorozmerny vektor, vrati
%   kratkodobu Fourierovu transformaciu vypocitanu zo signalu X. Pre
%   vypocet STFT je signal nasegmentovany do segmentov s dlzkou vektoru
%   OKNO s prekryvom PREKRYV. Kazdy segment je navahovany hodnotami vo
%   vektore OKNO. Vysledkom je matica S, ktorej riadky, resp. stlpce
%   odpovedaju hodnotam STFT pre rozne frekvencie, resp. casy.
%
%   S = STFT(X,OKNO,PREKRYV), kde X je viacrozmerne pole, prevedie vypocet 
%   STFT pre kazdy signal prvej dimenzie pola X.
%   
%   [S,f,t] = STFT(...) vrati vektor frekvencii f a vektor casov t,
%   pre ktore bola STFT napocitana.
%
%   Pokial signal nemoze byt rozdeleny na cely pocet segmentov, signal je
%   skrateny.
%
%   Priklad pouzitia c. 1
%   
%     % vygenerovanie signalov
%     fs=100;              % vzorkovacia frekvencia
%     n=1:2000;
%     x=zeros(length(n),3);
%     % signal obsahujuci harmonicku zlozku a dirakov impulz
%     x(:,1)=sin(2*pi*20*n/fs);  x(500,1)=20;  
%     % chirp signal
%     x(:,2)=sin(2*pi*n.^2/(fs.^2)); % chirp signal
%     % frekvencne modulovany signal
%     xm=10*sin(2*pi*0.2*n/fs); 
%     ph=filter(1,[1 -1],2*pi*(xm+25)/fs);
%     x(:,3)=sin(ph);
% 
%     % vypocet STFT
%     [S,f,t]=STFT(x,hamming(100),95,fs);
% 
%     % vykreslenie vysledkov
%     figure
%     for n=1:size(x,2)
%        subplot(size(x,2),1,n)
%        pcolor(t,f,abs(S(:,:,n)))
%        shading flat
%     end;
% 
%   Priklad pouzitia c. 2
% 
%     % vygenerovanie signalov
%     fs=100;              % vzorkovacia frekvencia
%     n=1:2000;
%     x=zeros(length(n),3,3);
%     for k1=1:size(x,2)
%         for k2=1:size(x,3)
%             x(:,k1,k2)=sin(2*pi*k1*k2*2*n/fs);
%             x(k1*k2*50,k1,k2)=20;
%         end
%     end;
% 
%     % vypocet STFT
%     [S,f,t]=STFT(x,hamming(100),75,fs);
% 
%     % vykreslenie vysledkov
%     figure
%     for k1=1:size(x,2)
%         for k2=1:size(x,3)
%             subplot(size(x,2),size(x,3),(k1-1)*size(x,2)+(k2-1)+1);
%             pcolor(t,f,abs(S(:,:,k1,k2)));
%             shading flat
%         end
%     end;


% Rev160629

w=w(:);
x_size=size(x);
if min(x_size)==1 && length(x_size)==2
    x=x(:); 
else 
    x=reshape(x,x_size(1),[]);
end;




N=size(x,1);
Nsig=size(x,2);
Lseg=length(w);

Nseg=floor((N-Lseg)/(Lseg-Noverlap))+1;

S=zeros(ceil(Lseg/2),Nseg,Nsig);

for n=1:Nsig

    [xb,~]=buffer(x(:,n),Lseg,Noverlap,'nodelay');

    xb=xb.*repmat(w,[1 size(xb,2)]);
    
    X=fft(xb);

    S(:,:,n)=X(1:ceil(Lseg/2),:);
    
end;

if min(x_size)~=1 || length(x_size)~=2
    S=reshape(S,[size(S,1), size(S,2), x_size(2:end)]);
end;

t=((0:Nseg-1)*(Lseg-Noverlap)+Lseg/2)/fs;
f=(0:ceil(Lseg/2)-1)/Lseg*fs;