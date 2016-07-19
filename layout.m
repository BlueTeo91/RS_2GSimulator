%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%      HEXAGONAL LAYOUT GENERATOR       %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all    % Close all open figures
clear        % Reset variables
clc          % Clear the command window


%% Generate hexagonal grid

N = 19;                             % NxN grid
Rad3Over2 = sqrt(3)/2;
[X, Y] = meshgrid(0:1:N);
n = size(X,1);
X = (Rad3Over2 * X)/N;
Y = (Y + repmat([0 0.5],[n,n/2]))/N;
radius = (1/sqrt(3))/N;  

hold on
plot(X,Y,'b.');                     % BS center
for i = 1:N;
    for j = 1:N
    xTemp = X(i,j);
    yTemp = Y(i,j);
    [xunit,yunit] = hexagon(xTemp,yTemp,radius);
    plot(xunit, yunit,'k');         % Hexagonal Cell borders
    end
end
axis equal, zoom on
