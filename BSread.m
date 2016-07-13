function [x,y] = BSread(file)
% Read BS coordinates from file
x = dlmread(file,'\t','A1..A:');
y = dlmread(file,'\t','B1..B:');
end