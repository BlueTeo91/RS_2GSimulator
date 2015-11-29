% Generate hexagonal grid
N = 21;
Rad3Over2 = sqrt(3) / 2;
[X, Y] = meshgrid(0:1:N);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);

% Plot the hexagonal mesh, including cell borders
[XV, YV] = voronoi(X(:),Y(:)); 
hold on
plot(X,Y,'.',XV,YV,'b-')
axis equal, axis([10 50 10 50]), zoom on

radius = 1/sqrt(3);

%Plot N*N circles centered in X,Y coordinates
for i=1:N
    for j=1:N
        circle(X(i,j),Y(i,j),radius)
    end
end

