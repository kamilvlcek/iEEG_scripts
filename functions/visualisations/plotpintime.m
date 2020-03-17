function [] = plotpintime(results, x, y)
%PLOTPINTIME plots a graph of p values with X axis of time
%   Allows plotting a size of a floating effect in timedomain comparisons.
%   Original author: Lukáš Hejtmánek
%   
%   results: 2D Matrix with time in the X axis and values corresponsing to 
%   x: dimensions for the x axis. default 1:size(results, 1)
%   y: dimensions for the y axis. default 1:size(results, 2)
%   example:
%       wilcox = CStat.Wilcox2D(dataA, dataB, 1, [], 'mean vs baseline');
%       plotpintime(wilcox, [0 500], [1 20])z
if ~exist('x', 'var'), x = 1:size(results, 1); end
if ~exist('y', 'var'), y = 1:size(results, 2); end

figure('Name', 'P value in time');
imagesc(x, y, 1 - results, [0.95 1]);
axis ij;
xlabel('time');
colorbar;
end

