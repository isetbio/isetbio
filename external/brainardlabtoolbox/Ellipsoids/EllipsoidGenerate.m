function xEllipsoid = EllipsoidGenerate(ellParams,nTheta,nPhi)
% xEllipsoid = EllipsoidGenerate(ellParams,nTheta,nPhi)
%
% Generate a set of points on the 3D ellipsoid by stretching and rotating
% points on the 3D unit sphere.
%
% See EllipsoidMatricesGenerate for specification of ellParams vector.
%
% 6/27/16  dhb  Wrote it.

xSphere = UnitSphereGenerate(nTheta,nPhi);
[A,Ainv,Q] = EllipsoidMatricesGenerate(ellParams);
xEllipsoid = Ainv*xSphere;
