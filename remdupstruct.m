function [evts2] =  remdupstruct(evts)
% REMDUPSTRUCT vrati strukturu bez duplikatu 
% pouzivam na odstraneni duplikatu z evts
    evts2 = struct(evts(1));
    for a = 2:numel(evts) %jedu po polozkach od druhe do konce
        duplikat = 0;
        for b = 1:a - 1 %porovnavam s predchozimi polozkami
            if isequal(evts(a),evts(b))
                duplikat = 1; %je nejaka predchozi stejna polozka, nechci ji
                break;
            end
        end
        if ~duplikat %nenasel jsem zadnou predchozi stejnou polozku kopirujue
            evts2(end+1) = evts(a); %#ok<AGROW>
        end
    end              
end

