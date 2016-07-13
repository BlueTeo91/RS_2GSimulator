function d = computeDistance(coordsA, coordsB)
% Computes distance between MS and BS
x1 = coordsA(:,1);
y1 = coordsA(:,2);
x2 = coordsB(:,1);
y2 = coordsB(:,2);

d = sqrt((x2-x1).^2+(y2-y1).^2);
end