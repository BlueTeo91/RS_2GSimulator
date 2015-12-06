% Generate hexagonal grid
N = 19;
Rad3Over2 = sqrt(3) / 2;
[X, Y] = meshgrid(0:1:N);
n = size(X,1);
X = (Rad3Over2 * X)/N;
Y = (Y + repmat([0 0.5],[n,n/2]))/N;

% Plot the hexagonal mesh, including cell borders
[XV, YV] = voronoi(X(:),Y(:)); 
hold on
plot(X,Y,'.',XV,YV,'b-')
axis equal, axis([0.1 1 0.1 1]), zoom on

radius = (1/sqrt(3))/N;

% Plot N*N circles centered in X,Y coordinates
for i=1:N
    for j=1:N
        circle(X(i,j),Y(i,j),radius)
    end
end
