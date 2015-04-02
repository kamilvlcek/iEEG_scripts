function h=circle(r,x,y,format,filled)
% circle(r,x,y,format), eg. circle(140,0,0,'k')
% jako druhy argument se muze vlozit i complexni cislo
% v tom pripade se pouzivaji jen 3 argumenty
if(~isreal(x)) %jako druhy argument se muze vlozit i complexni cislo
    % v tom pripade se pouzivaji jen 3 argumenty
    format = y;
    y = imag(x);
    x = real(x);
end
if ~exist('filled','var') , filled = 0; end %7.11.2014
if filled
    N = 126;
    t = 2*pi/N*(1:N); 
    fill(x+r*cos(t), y+r*sin(t), format);
    
else
    plot(r*exp(i*[0:0.05:2*pi+0.05])+x+y*i,format)
end