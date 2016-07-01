classdef CStat
    %CSTAT Class for custom statistical functions
    %   Kamil Vlcek, FGU AVCR, since 2016 06
    
    properties
    end
    
    methods (Static,Access = public)
        
        function W = Wilcox2D(A,B,print,fdr)
            %srovna dve 3D matice proti sobe, ohledne hodnot v poslednim rozmeru
            %A musi mit oba prvni rozmery > rozmery B, 
            %B muze mit jeden nebo oba prvni rozmer = 1 - pak se porovnava se vsemi hodnotami v A
            %pokud fdr=1 nebo prazne, provadi fdr korekci
            %pokud print = 0 nebo prazne, netiskne nic
            if ~exist('print','var'), print = 0; end
            if ~exist('fdr','var'), fdr = 1; end            
            
            W = zeros(size(A,1),size(A,2));
            if print, fprintf('Wilcox Test 2D, 1.dim :'); end
            for j = 1:size(A,1)
                if print && mod(j,50)==0, fprintf('%d ', j); end %tisknu jen cele padesatky
                for k = 1:size(A,2)                   
                   aa = squeeze (A(j,k,:)); 
                   bb = squeeze (B( min(j,size(B,1)) , min(k,size(B,2)) , :));
                   W(j,k)= ranksum(aa,bb);
                end
            end
            if fdr
                if fdr == 2, method = 'dep'; else method = 'pdep'; end
                [~, ~, adj_p]=fdr_bh(W,0.05,method,'no'); %dep je striktnejsi nez pdep
                W = adj_p; %prepisu puvodni hodnoty korigovanymi podle FDR
            end
            if print, fprintf('%d .. done\n',j); end
        end
        
        function [W] = Klouzaveokno(W,oknosirka, funkce,dimension)
            % oknosirka je v poctu bodu, funkce muze byt min, max, mean
            % 27.4.2015 - vynato z wilcoxmap, 21.6.2016 presunuto do CStat      
                        
            if oknosirka >= 1 
                if dimension == 1 %pokud chci pocitat klouzave okno v prvnim rozmeru, transponuju na zacatku i na konci
                    W = W';
                end                
                pulsirka = ceil(oknosirka/2); %ktery sloupec se povazuje za pulku okna - tam se hodnota ulozi
                W2 = zeros(size(W,1),size(W,2)); %musim udelat kopii, jinak si prepisuju hodnoty ze kterych pak pocitam
                for sloupec = 1:size(W,2); 
                    iW = max([1 sloupec-pulsirka+1]) : min([size(W,2) sloupec-pulsirka+oknosirka]); 
                    switch funkce
                        case 'min'
                            W2(:,sloupec)=min(W(:,iW),[],2);
                        case 'max'
                            W2(:,sloupec)=max(W(:,iW),[],2);
                        case 'mean'
                            W2(:,sloupec)=mean(W(:,iW),2);
                        otherwise
                            W2(:,sloupec)=W(:,sloupec); %pokud jina funkce, nemenim matici
                    end 
                end
                if dimension == 1 
                    W = W2';
                else
                    W = W2;
                end
            end
        end
        
        function N = round(n,dig)
            %zakrouhuje na dany pocet desetinnych mist
            N = round(n * 10^dig) / 10 ^ dig;
        end
            
    end
    
end

