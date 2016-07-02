function x = UnitSphereGenerate(nTheta,nPhi)
% x = UnitSphereGenerate(nTheta,nPhi)
%
% Generate a set of points on the unit sphere in three
% dimensions.
%   nTheta - number of samples around the azimuthal theta (0 to 2pi)
%   nPhi   - number of points around the elevation direction (0 to pi);
%
% Coordinates are returned in an 3 by(nTheta*nPhi) matrix, with the rows
% being the x, y, z coordinates of the points.
%
% Follows code posted at
%   https://www.mathworks.com/matlabcentral/answers/44481-how-to-pick-n-random-points-on-a-sphere
%
% You could also use the Matlab fuction sphere() to do this.
% 
% 6/27/16  dhb  Wrote it.

% Generate a unit sphere in 3D
thetaBase=linspace(0,2*pi,nTheta);
phiBase=linspace(0,pi,nPhi);
[theta,phi]=meshgrid(thetaBase,phiBase);
rho=1;
xCoords=rho*sin(phi(:)).*cos(theta(:));
yCoords=rho*sin(phi(:)).*sin(theta(:));
zCoords=rho*cos(phi(:));

% Stuff the coordinates into a single nTheta*nPhi by 3 matrix
x = [xCoords yCoords zCoords]';
