function [x,y] = uniformMS(xa,xb,ya,yb,n)
% Generate n random points in the 2 coordinate system
% with coordinate x ranging from xa to xb
% and coordinate y ranging from ya to yb

x = xa + (xb - xa) .* rand(n, 1);
y = ya + (yb - ya) .* rand(n, 1);

end
