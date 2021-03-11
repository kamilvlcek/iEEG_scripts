function [b] = cell2double(a)
%CELL2DOUBLE converts cell array with intermixed double values and char arrays to double array
%   Detailed explanation goes here
if iscell(a)
    b = zeros(size(a));
    for ia = 1:numel(a)
        if ischar(a{ia}) || isstring(a{ia})
            b(ia) = str2double(a{ia});
        else
            b(ia) = a{ia};
        end
    end
end
end

