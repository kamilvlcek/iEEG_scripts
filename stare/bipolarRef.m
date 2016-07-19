function [ dataB ] = bipolarRef( data,els,prumery )
%BIPOLARREF vrati data s bipolarni nebo prumernou (za elektrody) referenci 
%  pokud prumery=0, bipolarni reference
%  pokud prumery=1, udela prumernou referenci
%  pokud prumery=2, udela celkove prumernou referenci,za vsechny elektrody
fprintf('bipolarRef: '); 
dataB = zeros(size(data,1),size(data,2),size(data,3)); %dimenze:elektrody,time,events
e0=1;               %nejnizsi cislo kontaktu aktualni elektrody
for e=els           %cyklus pres jednotlive elektrody
                    %e je nejvyssi cislo kontaktu pro tutu elektrodu
    fprintf('El %i, ',e);                
    if nargin < 3 || prumery==0
        %bipolarni reference
        for j=e0:e-1    %cyklus pres kontakty v elektrody, od nejnizsiho 
            dataB(j,:,:)=data(j,:,:)-data(j+1,:,:); %odecitam od aktualni elektrody vyssi cislo el
            %protoze vyssi cislo znamena vice k povrchu mozku
            %dataB(e,:,:) se vynechavam tam uz je 0
        end
        e0 = e+1;       %nejnizsi kontakt nasledujici elektrody
    elseif prumery==1
        %prumerna reference za kazdou elektrodu zvlast - 16.9.2014
        eavg = mean(data(e0:e,:,:));     %prumer za jednu elektrodu (pro kazdy cas a event zvlast)          
        for j=e0:e    %cyklus pres kontakty v elektrody, od nejnizsiho 
            dataB(j,:,:)=data(j,:,:)-eavg(1,:,:); %odecitam od aktualni elektrody vyssi cislo el
            %protoze vyssi cislo znamena vice k povrchu mozku
            %dataB(e,:,:) se vynechavam tam uz je 0
        end
        e0 = e+1;       %nejnizsi kontakt nasledujici elektrody
    elseif prumery==2
        %celkove prumerna reference
        eavg = mean(data(:,:,:));
        for j=e0:e    %cyklus pres kontakty v elektrody, od nejnizsiho 
            dataB(j,:,:)=data(j,:,:)-eavg(1,:,:); %odecitam od aktualni elektrody vyssi cislo el
            %protoze vyssi cislo znamena vice k povrchu mozku
            %dataB(e,:,:) se vynechavam tam uz je 0
        end
        e0 = e+1;   
    end
end
disp('OK');
end

