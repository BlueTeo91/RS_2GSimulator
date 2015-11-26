% Generate hexagonal grid
Rad3Over2 = sqrt(3) / 2;
[X Y] = meshgrid(0:1:41);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);

for i=1:42
    for j=1:42
       
        hold on
    end
end

% Plot the hexagonal mesh, including cell borders
[XV YV] = voronoi(X(:),Y(:)); 
plot(X,Y,'.',XV,YV,'b-')
axis equal, axis([10 100 10 100]), zoom on