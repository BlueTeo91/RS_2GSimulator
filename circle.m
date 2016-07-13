function circle(x,y,r)
%Plot a circle with radius 'r' and center at the coordinates 'x' and 'y'

hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit,'k');
axis equal

end

