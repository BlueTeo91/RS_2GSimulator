function p = random(a,b,n)
% Plot n random points with coordinates (x,y) ranging from a to b
hold on
% Generate n random points in the 2 coordinate system 
x = a + (b - a) * rand(1, n); 
y = a + (b - a) * rand(1, n); 

% plot
p = plot(x, y, 'o');
axis equal
end

