function xEllipsoid = PointsOnEllipsoidFind(Q,xDir,xCenter)
% xEllipsoid = PointsOnEllipsoidFind(Q,xDir)
%
% Given the quadratic matrix Q that defines an ellipsoid and a set of
% vector directions, find the point in each direction on the ellipsoid.
%
% If 3-dimensional column vector xCenter is passed, add that back into the
% computed xEllipsoid.  Note, however, that xDir should be directions
% around the origin.
%
% x, y, and z coordinates are in the rows of xDir and xEllipsoid.
%
% 7/2/16  dhb  Wrote it.

if (nargin < 3 | isempty(xCenter))
    xCenter = [0 0 0]';
end

%% Find lengths in each direction
theLengths = diag(xDir'*Q*xDir);

%% Scale
xEllipsoid = xDir*diag(1./sqrt(theLengths));

%% Add back in center
for ii = 1:3
    xEllipsoid(ii,:) = xEllipsoid(ii,:) + xCenter(ii);
end

