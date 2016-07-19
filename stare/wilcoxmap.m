function [ W ] = wilcoxmap( ALLEEG, pCfg,testname,blok,W,range)
% blok: cisla datasetu EEGlabu k porovnani. Wilcox test.
% Pokud jeden blok, dela se stat vuci prumeru pred podnete
% Pokud dva bloky, dela se stat dvou matic vuci sobe
% pCfg: konfigurace pacienta
% testname, jmeno testu v konfiguraci pacienta
% W: pokud se neuda zadny blok a muze se pouzit a pouze vykreslit drive vypocitane W
% range: cisla elektrod k vypoctu

els = pCfg.els;
cas = pCfg.(testname).cas;

if numel(blok)==1
    CondA = ALLEEG(blok).data; %dimensions: channel x time x epochs 
else
    CondA =  ALLEEG(blok(1)).data;   
    CondB =  ALLEEG(blok(2)).data;
end

stimul = ceil(cas(1)/1000*-512); %vzorek frekvence 512 Hz

if ~exist('range','var') 
    range = 1:ALLEEG(blok(1)).nbchan; %rozsah elektrod na vykresleni
end
   
if numel(blok)==1
    %obrazek hrubych dat - 3.11.2014
    figure('Name',['prumerne odpovedi ' pCfg.name ', ' testname]);
    M = mean(CondA,3); %mean over all epochs
    imagesc(cas,range,M(range,:));
    colorbar;
    caxis([-30,30]);
    carydografu(range,cas,els);
end

%bipolarni reference
reference = 0; %0=bipolarni, -1 zadna
disp(['reference ' iff(reference==1,'prumerna El',iff(reference==2,'prumerna celkove',iff(reference==0,'bipolarni','zadna')))]);
if reference >= 0
    CondA = bipolarRef(CondA,els,reference);  
end
if exist('CondB','var') , CondB = bipolarRef(CondB,els,reference);   end %12.11.2014


% if(sum(blok==AlloI)>0),     Allo = bipolarRef(Allo,els,reference);          end
% if(sum(blok==EgoI)>0),      Ego = bipolarRef(Ego,els,reference);            end
% if(sum(blok==ControlI)>0),  Control = bipolarRef(Control,els,reference);    end

if numel(blok)==1
    %obrazek hrubych dat - 3.11.2014
    figure('Name',['prumerne odpovedi bipolarni ' pCfg.name ', ' testname])
    M = mean(CondA,3);
    imagesc(cas,range,M(range,:));
    colorbar;
    caxis([-30,30]);
    carydografu(range,cas,els);
end

%pridam funkci na orezani ALLO podle casu odpovedi -600 a +900ms kolem ni
%Allo = vyrovnejresponse(Allo,responsetime(ALLEEG,AlloI),512,-cas(1)/1000);
%Ego = vyrovnejresponse(Ego,responsetime(ALLEEG,EgoI),512,-cas(1)/1000);
%Control = vyrovnejresponse(Control,responsetime(ALLEEG,ControlI),512,-cas(1)/1000);
%cas = [-200 1500]; %600+900

%mapu signifikanci vuci casu pred stimulem
if numel(blok)==1
    fdr = 1; % pouzivam fdr korekci
    C = mean(CondA(:,1:stimul-1,:),2);
    W = WilcoxA(CondA,C); 
    
%     if blok == 0 && nargin >=4
%         %nic, pouziju W z argumentu
%         disp('pozivam zadane W');
%     elseif blok == 1
%         C = mean(Allo(:,1:stimul-1,:),2);
%         W = WilcoxA(Allo,C); %Wa = W; 
%     elseif blok == 2
%         C = mean(Ego(:,1:stimul-1,:),2);
%         W = WilcoxA(Ego,C); %We = W; 
%     elseif blok == 3
%         C = mean(Control(:,1:stimul-1,:),2);
%         W = WilcoxA(Control,C); %Wc = W;  
%     else
%         error('nezname cislo bloku');
%         %exit;
%     end

elseif numel(blok)==2
    fdr = 1; %pro tento kontrast budu/nebudu pouzivat fdr korekci
    W = WilcoxM(CondA,CondB,fdr);
%     if isequal(blok,[1 2])
%         W = WilcoxM(Allo,Ego,fdr);
%     elseif isequal(blok,[1 3])
%         W = WilcoxM(Allo,Control,fdr);
%     elseif isequal(blok,[2 3])
%         W = WilcoxM(Ego,Control,fdr);
%     else
%         error('neznama dvojice bloku');
%         %exit;
%     end
   
elseif ~exist('W','var') %pokud existuje W, pouziju tu
    error('trojice bloku nejdou');
    %exit;
end

WK = klouzaveokno(W,25,'max'); %25 nebo 50
%WK = W;

%nakreslim 2D obrazek signifikanci
figure('Name',['W map ' pCfg.name ', ' testname]);
imagesc(cas,range, 1-WK(range,:),[ iff(fdr,0.95,0.99) 1]);%mapa, od p>0.05 bude modra barva 
colorbar;
WKmin = find(min(WK,[],2)< iff(fdr,0.05,0.01) );
disp (['elektrody se signifikancemi: ' num2str(WKmin')]) 
disp (['p>' num2str(min(min(WK)))]);

for m = WKmin'
    text(800,m,num2str(m),'Color',[0.5 0.5 0.5]);
end
    
carydografu(range,cas,els)
%set(gca,'YTickLabel',cellstr(int2str([1:30]'))');

end


function []=carydografu(range,cas,els)
    %cara casu 0
    line([0 0],[min(range) max(range)], 'Color',[0.5 0.5 0.5]);
    %cary hranic elektrod
    for e = els
        line(cas,[e e], 'Color',[0.5 0.5 0.5]);
    end
end
