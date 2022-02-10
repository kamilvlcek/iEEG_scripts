function [] = plotpintime(results, x, y)
%PLOTPINTIME plots a graph of p values with columns representing timeseries
%   of different variables. Allows plotting a size of a floating effect 
%   in timedomain comparisons.
%   Original author: Lukáš Hejtmánek
%   
%   results: 2D Matrix with each column corresponding to a single variable
%       timeseries along the x axis. E.g. if plotting p values for various
%       frequencies, each column in the results represents a single
%       frequency and each row is a point in time
%   x: dimensions for the x axis. default 1:size(results, 1)
%   y: dimensions for the y axis. default 1:size(results, 2)
%   example:
%       wilcox = CStat.Wilcox2D(dataA, dataB, 1, [], 'mean vs baseline');
%       plotpintime(wilcox, [0 500], [1 20])z
if ~exist('x', 'var'), x = 1:size(results, 1); end
if ~exist('y', 'var'), y = 1:size(results, 2); end

figure('Name', 'P value in time');
% Transposes results - from each column is timeseries for a particular
% variable to each column is a point in time fro multiple variables
plt = imagesc(x, y, (1 - results)', [0.95 1]);
set(plt, 'AlphaData', ~isnan(results'));
axis xy;
xlabel('time');
colorbar;
end

