function [b] = cell2double(a)
%CELL2DOUBLE converts cell array with intermixed double values and char arrays to double array
%   also if cell contains multiple double vectors, it explodes them to the main vector, similarly to cell2mat
%   if a is not a cell, it returns the same value
if iscell(a)
    numa = 0;
    for ia = 1:numel(a)
        if isnumeric(a{ia})
            numa = numa + numel(a{ia});
        else
            numa = numa + 1;
        end
    end
    b = iff(isrow(a),zeros(1,numa),zeros(numa,1));
    ib = 1;
    for ia = 1:numel(a)
        if ischar(a{ia}) || isstring(a{ia})
            b(ib) = str2double(a{ia});
            ib = ib + 1;
        elseif isnumeric(a{ia})
            b(ib:ib+numel(a{ia})-1) = a{ia};
            ib = ib + numel(a{ia});
        else
            b(ib) = a{ia};
            ib = ib + 1;
        end
    end
else
    b = a;
end
end

