function x = UnitCircleGenerate(nTheta)
% x = UnitCircleGenerate(nTheta)
%
% Generate a set of points on the unit circle in two dimensions.
%   nTheta - number of samples around the azimuthal theta (0 to 2pi)
%
% Coordinates are returned in an 2 by (nTheta) matrix, with the rows
% being the x, y coordinates of the points.
% 
% 7/2/16  dhb  Wrote it.

% Generate a unit sphere in 3D
theta=linspace(0,2*pi,nTheta);
rho=1;
xCoords=rho*cos(theta(:));
yCoords=rho*sin(theta(:));

% Stuff the coordinates into a single nTheta by 2 matrix
x = [xCoords yCoords]';
