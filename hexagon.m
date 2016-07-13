function [xunit,yunit] = hexagon(x,y,r)
th = 0:pi/3:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
end