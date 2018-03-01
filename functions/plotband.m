function hh = plotband(x,y,ey,bcolor,transparency,varargin)
% function plotband - plot graph with errorband
%
% syntax: h = plotband(x,y,ey,bcolor,line-properties)
%
% draws a graph centered to a band of deviations (e.g. standard error)
% x: vector points on the x-Axis, where y is provided
% y: vector, the same length as x
% ey: deviation from y, that will determine the witdh of the band. The
% lower boundary of the band will be at y-ey and the upper one at y+ey.
% Should be a vector, the same size as x and y, or a scalar.
% bcolor: (optional) color of the band, either as a 3-element color def. or
% one of the usual color strings (like 'r' for red, 'g' for green).
% For [], the default color, a nice shade of lavender, is used.
% line-properties: any line propteries that can be given with the usual
% plot() command, which will be applied to the central 
% line of y-data.
% poslano od JIRKY HAMMERA 13.3.2017

if nargin<5
    transparency = 0.5;
end
if nargin<4 || isempty(bcolor)
    bcolor= [0.7 0.5 1];
end

x= x(:);
y= y(:);
ey= ey(:);

ype= y+ey;
yme= y-ey;

washold = ishold;

h(1,1)= fill([x;x(end:-1:1)],[ype;yme(end:-1:1)], bcolor, 'edgecolor','none');
alpha(h(1,1), transparency);
hold on
h(2,1)= plot(x,y,'color',bcolor, varargin{:}); % co je tohle?

%if nargout>0, varargout= {h}; end %puvodni vystupni parametr byl varargout
hh = h(1,1); %kamil 1.3.2018
if ~washold, hold off, end 