function [SEM] = stderr(data)
%SEM - standard error of the mean
% ignors NaN values using nanstd function

SEM = nanstd(data,1)/sqrt(size(data,1));