function xEllipsoid = PointsOnEllipsoidFind(Q,xDir)
% xEllipsoid = PointsOnEllipsoidFind(Q,xDir)
%
% Given the quadratic matrix Q that defines an ellipsoid and a set of 
% vector directions, find the point in each direction on the ellipsoid.
%
% x, y, and z coordinates are in the rows of xDir and xEllipsoid.
%
% 7/2/16  dhb  Wrote it.

%% Find lengths in each direction
theLengths = diag(xDir'*Q*xDir);

%% Scale
xEllipsoid = xDir*diag(1./sqrt(theLengths));

