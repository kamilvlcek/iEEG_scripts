function h=circle(r,x,y,format,filled)
% circle(r,x,y,format), eg. circle(140,0,0,'k')
% jako druhy argument se muze vlozit i complexni cislo
% v tom pripade se pouzivaji jen 3 argumenty
% jako druhy argument muze byt i vektor 3 cisel -> plot in 3D
if ~isreal(x)   %jako druhy argument se muze vlozit i complexni cislo 
    % v tom pripade se pouzivaji jen 3 argumenty    
    format = y;
    y = imag(x);
    x = real(x);
end
if ~exist('filled','var') , filled = 0; end %7.11.2014
if filled %plot filled circle in 2D
    N = 126;
    t = 2*pi/N*(1:N); 
    h=fill(x+r*cos(t), y+r*sin(t), format);
elseif numel(x) == 3 % if the second argument is 3 numbers -> plot the circle in 3D
    % the x contains the center, the y contains code of the plane 
    % circle in 3D cant be filled
    teta=-pi:0.05:pi+0.05; %angles
    vals = [r*cos(teta); r*sin(teta)]; %values of the circle
    if y == 1 %circle in x plane 
        h=plot3(x(1)+zeros(1,numel(teta)), x(2)+vals(1,:), x(3)+vals(2,:),format);   
    elseif y == 2 %circle in y plane
        h=plot3(x(1)+vals(1,:), x(2)+zeros(1,numel(teta)) , x(3)+vals(2,:),format);
    elseif y == 3 %circle in z plane
        h=plot3(x(1)+vals(1,:), x(2)+vals(2,:) , x(3)+zeros(1,numel(teta)),format);        
    end
else %plot empty circle in 2D 
    h=plot(r*exp(1i* 0:0.05:2*pi+0.05 )+x+y*1i,format);
end