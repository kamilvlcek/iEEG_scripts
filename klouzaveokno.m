function [W2] = klouzaveokno(W,oknosirka, funkce)
% oknosirka je v poctu bodu, funkce muze byt min, max, mean
% 27.4.2015 - vynato z wilcoxmap

    funkce_handle = str2func(funkce);
    if oknosirka == 0 
        W2 = W;
    else
        %oknosirka = 10; %sirka klouzaveho okna
        pulsirka = ceil(oknosirka/2); %ktery sloupec se povazuje za pulku okna - tam se hodnota ulozi
        W2 = zeros(size(W,1),size(W,2)); %musim udelat kopii, jinak si prepisuju hodnoty ze kterych pak pocitam
        for sloupec = 1:size(W,2); 
            iW = max([1 sloupec-pulsirka+1]) : min([size(W,2) sloupec-pulsirka+oknosirka]); 
            switch funkce
                case 'min'
                case 'max'
                   W2(:,sloupec)=funkce_handle(W(:,iW),[],2);
                case 'mean'
                   W2(:,sloupec)=funkce_handle(W(:,iW),2);
                otherwise
                   W2(:,sloupec)=W(:,sloupec);
            end 
        end
    end
end