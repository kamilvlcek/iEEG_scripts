function [ val ] = cellval( kat, n )
%CELLVAL vrati hodnotu n pole kat, at uz to je cellarray nebo numeric array

    if iscell(kat)
        val = kat{n};
    else
        val = kat(n);
    end

end

