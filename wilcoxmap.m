function [ W ] = wilcoxmap( ALLEEG, pacient,blok,W,range)
% promenne podle soucasnych datasetu v eeglabu
%blok: 1=allo vs baseline, 2=ego vs baseline, 3=control vs baseline, 1 2, 1 3, 2 3

if  pacient == 73 %p73
    els = [ 11 21 31 43 53 64 75 92 101]; %hranice elektrod pro p71 - konce
    %cas = [-200 2800]; %cas = [-500 1500];
    cas = [-100 900]; %PPA experiment
elseif  pacient == 71 %p71
    els = [ 16 32 48  64 77 89 105 115 125]; %hranice elektrod pro p71 - konce
    cas = [-200 2800]; %cas = [-500 1500];
elseif pacient == 69 %p69 - AEdist-p69 1Hz Ep01-1 Allo.set, AEdist-p69 1Hz Ep01-1 Control.set, AEdist-p69 1Hz Ep01-1 Ego.set
    els = [12 20 27 37 47 56 64 75 83 90 99 109 119 125]; %p69
    cas = [-100 1000];
elseif pacient == 68 %p68
    els = [ 10 18 28  41 54 64 72 82 91 100 114 125]; %hranice elektrod pro p71 - konce
    cas = [-200 1000];
else
    disp('nezname cislo pacienta');
    exit;
end


if numel(blok)==1
    CondA = ALLEEG(blok).data;    
else
    CondA =  ALLEEG(blok(1)).data;   
    CondB =  ALLEEG(blok(2)).data;
end

% AlloI = 5; EgoI = 4; ControlI = 2;
% if(sum(blok==AlloI)>0),     Allo = ALLEEG(AlloI).data;          end %v PPA to je Scene
% if(sum(blok==EgoI)>0),      Ego = ALLEEG(EgoI).data;            end %v PPA to je Face
% if(sum(blok==ControlI)>0),  Control = ALLEEG(ControlI).data;    end % v PPA to je Object

stimul = ceil(cas(1)/1000*-512); %vzork frekvence 512 Hz

if nargin < 5
    range = 1:ALLEEG(blok(1)).nbchan; %rozsah elektrod na vykresleni
end

if numel(blok)==1
    %obrazek hrubych dat - 3.11.2014
    figure('Name','prumerne odpovedi')
    M = mean(CondA,3);
    imagesc(cas,range,M(range,:));
    colorbar;
    caxis([-30,30]);
    carydografu(range,cas,els);
end

%bipolarni reference
reference = 0; %0=bipolarni
disp(['reference ' iff(reference==1,'prumerna El',iff(reference==2,'prumerna celkove','bipolarni'))]);

CondA = bipolarRef(CondA,els,reference);  
if exist('CondB','var') , CondB = bipolarRef(CondB,els,reference);   end %12.11.2014


% if(sum(blok==AlloI)>0),     Allo = bipolarRef(Allo,els,reference);          end
% if(sum(blok==EgoI)>0),      Ego = bipolarRef(Ego,els,reference);            end
% if(sum(blok==ControlI)>0),  Control = bipolarRef(Control,els,reference);    end

if numel(blok)==1
    %obrazek hrubych dat - 3.11.2014
    figure('Name','prumerne odpovedi bipolarni')
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
    fdr = 0; %pro tento kontrast nebudu pouzivat fdr korekci
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

WK = klouzaveokno(W,25); %25 nebo 50
%WK = W;

%nakreslim 2D obrazek signifikanci
figure('Name','W map');
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
function [W2] = klouzaveokno(W,oknosirka)
%oknosirka je v poctu bodu
    if oknosirka == 0 
        W2 = W;
    else
        %oknosirka = 10; %sirka klouzaveho okna
        pulsirka = ceil(oknosirka/2); %ktery sloupec se povazuje za pulku okna - tam se hodnota ulozi
        W2 = zeros(size(W,1),size(W,2)); %musim udelat kopii, jinak si prepisuju hodnoty ze kterych pak pocitam
        for sloupec = 1:size(W,2); 
            iW = max([1 sloupec-pulsirka+1]) : min([size(W,2) sloupec-pulsirka+oknosirka]); 
            W2(:,sloupec)=max(W(:,iW),[],2);
        end
    end
end

function []=carydografu(range,cas,els)
    %cara casu 0
    line([0 0],[min(range) max(range)], 'Color',[0.5 0.5 0.5]);
    %cary hranic elektrod
    for e = els
        line(cas,[e e], 'Color',[0.5 0.5 0.5]);
    end
end
