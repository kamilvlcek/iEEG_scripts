function [S,f,t]=KamilovaTransformacia(x,freq,fs)

%KamilovaTransformacia rozdeli signal na vybrane frekvencne pasma 
%
%   S = KamilovaTransformacia(X,FREKV,Fs), kde X je jednorozmerny vektor,
%   rozdeli signal na N-1 frekvencnych pasiem, s hranicami
%   <FREKV(1),FREKV(2)>, <FREKV(2),FREKV(3)>, ... ,<FREKV(N-1),FREKV(N)>,
%   kde N je dlzka vektoru FREKV. Vektor FREKV musi byt rastuca postupnost
%   realnych cisel v intervale <0,Fs/2>, kde Fs je vzorkovacia frekvencia
%   signalu. Filtracia je prevedena pomocou FIR filtrov s nulovou fazou
%   navrhnutych pomocou metody okien s pouzitim Hammingovho okna. Vysledkom
%   je matica S, ktorej riadky obsahuju vystupy jednotlivych filtrov.
%     
%   S = KamilovaTransformacia(X,FREKV,Fs), kde X je viacrozmerne pole,
%   prevedie vypocet filtracie pre pre kazdy signal prvej dimenzie X.
%
%   [S,f,t] = KamilovaTransformacia(...) vrati vektor pseudofrekvencii f
%   a vektor casov t, pre ktore bola waveletova transformacia napocitana.
%
%   Priklad pouzitia c. 1
%
%     % vygenerovanie signalov
%     fs=100;              % vzorkovacia frekvencia
%     n=1:2000;
%     x=zeros(length(n),3);
%     % signal obsahujuci harmonicku zlozku a dirakov impulz
%     x(:,1)=sin(2*pi*20*n/fs);  x(500,1)=10;  
%     % chirp signal
%     x(:,2)=sin(2*pi*n.^2/(fs.^2)); 
%     % frekvencne modulovany signal
%     xm=10*sin(2*pi*0.2*n/fs); 
%     ph=filter(1,[1 -1],2*pi*(xm+25)/fs);
%     x(:,3)=sin(ph);
% 
%     % vypocet "transformacie"
%     freq=0:1:fs/2;
%     [S,f,t]=KamilovaTransformacia(x,freq,fs);
%     S=OkamzitaAmplituda(permute(S,[2,1,3]));
% 
%     % vykreslenie vysledkov
%     figure
%     for k=1:size(x,2)
%        subplot(size(x,2),1,k)
%        pcolor(t,f,abs(S(:,:,k))')
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
%             x(k1*k2*50,k1,k2)=10;
%         end
%     end;
% 
%     % vypocet "transformacie"
%     freq=0:1:fs/2;
%     [S,f,t]=KamilovaTransformacia(x,freq,fs);
%     S=OkamzitaAmplituda(permute(S,[2 1 3 4]));
% 
%     % vykreslenie vysledkov
%     figure
%     for k1=1:size(x,2)
%         for k2=1:size(x,3)
%             subplot(size(x,2),size(x,3),(k1-1)*size(x,2)+(k2-1)+1);
%             pcolor(t,f,abs(S(:,:,k1,k2))');
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


N=size(x,1);
Nsig=size(x,2);

S=zeros(length(freq)-1,N,Nsig);

for n=1:Nsig    
    for k=1:length(freq)-1
        Nord=round(fs/(freq(k+1)-freq(k))/2);
        if (freq(k)==0)
           w=fir1(Nord,freq(k+1)/fs*2);
        elseif (freq(k+1)==fs/2)
           w=fir1(Nord,freq(k)/fs*2,'high');
        else 
           w=fir1(2*Nord,[freq(k)/fs*2,freq(k+1)/fs*2]);
        end;
        S(k,:,n)=conv(x(:,n),w,'same');
    end
end

if min(x_size)~=1 || length(x_size)~=2
    S=reshape(S,[size(S,1), size(S,2), x_size(2:end)]);
end;

f=freq(1:end-1)+diff(freq)/2;
t=(0:N-1)/fs;