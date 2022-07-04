function [ val ] = cellval( kat, n )
%CELLVAL vrati hodnotu n pole kat, at uz to je cellarray nebo numeric array
    if exist('n','var')
        if iscell(kat)
            val = kat{n};
        else
            val = kat(n);
        end
    else
        if iscell(kat)
            kat = cell2mat(kat);
        end
        val = kat;
    end
        

end

