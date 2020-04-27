function [str] = double2str(n,decimals, delimiter)
%DOUBLE2STR converts double to character array
%   enables to customize the delimiter, uses sprintf + strrep
if ~exist('decimals','var') || isempty(decimals) , decimals = 1; end
if ~exist('delimiter','var') || isempty(delimiter) , delimiter = ','; end
str = sprintf(['%.' num2str(decimals) 'f'],n);
str = strrep(str,'.',delimiter);
end

