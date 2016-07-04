function EllipsoidTest
% EllipsoidTest
%
% Test the ellipsoid code
%
% 6/27/16  dhb  Wrote it.

%% Clear
clear; close all;

%% Generate points on unit sphere for simulating data
nTheta = 10;
nPhi = 10;
xSphere = UnitSphereGenerate(nTheta,nPhi);
figure; clf; hold on
plot3(xSphere(1,:),xSphere(2,:),xSphere(3,:),'ro','MarkerSize',8,'MarkerFaceColor','r');
axis('square');
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Unit Sphere');

%% Generate the ellipsoid matrices from a vector
% of ellipsoid parameters
ellParams = [3 0.8 2 pi/4 pi/8 pi]';
[A,Ainv,Q] = EllipsoidMatricesGenerate(ellParams);

%% Map the unit sphere into the ellipsoid
xEllipsoid = Ainv*xSphere;
figure; clf; hold on
plot3(xEllipsoid(1,:),xEllipsoid(2,:),xEllipsoid(3,:),'ro','MarkerSize',8,'MarkerFaceColor','r');
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Ellipsoid');

%% Make sure individual ellipsoid points satisfy xEllipsoid'*Q*xEllipsoid = 1
check = diag(xEllipsoid'*Q*xEllipsoid);
tolerance = 1e-8;
if (max(abs(check(:)-1)) > tolerance)
    error('Failure of ellipsoid points to have transformed length of 1');
end

%% Map ellipsoid back into unit sphere space and make sure it is a sphere
xSphereCheck = A*xEllipsoid;
figure; clf; hold on
plot3(xSphereCheck(1,:),xSphereCheck(2,:),xSphereCheck(3,:),'ro','MarkerSize',8,'MarkerFaceColor','r');
axis('square');
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Unit Sphere Check');

%% Add some noise to the ellipsoid, so we can try to fit it
noiseSd = 0.2;
xEllipsoid = EllipsoidGenerate(ellParams,nTheta,nPhi);
xNoisyEllipsoid = xEllipsoid + normrnd(0,noiseSd,size(xEllipsoid));
figure; clf; hold on
plot3(xNoisyEllipsoid(1,:),xNoisyEllipsoid(2,:),xNoisyEllipsoid(3,:),'ro','MarkerSize',8,'MarkerFaceColor','r');
axis('square');
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Noisy Ellipsoid With Fit');

%% Fit that sucker
%
% Options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');

% Have a go at reasonable initial values
ellRanges = max(xNoisyEllipsoid,[],2)-min(xNoisyEllipsoid,[],2);
ellParams0 = [ellRanges' 0 0 0]';

% Set reasonable bounds on parameters
vlb = [0.1 0.1 0.1 0 0 0]';
vub = [1e3 1e3 1e3 2*pi 2*pi 2*pi]';

% Fit
nThetaFit = 50;
nPhiFit = 50;
ellParamsFit = fmincon(@(x)FitEllipseFunction(x,xNoisyEllipsoid),ellParams0,[],[],[],[],vlb,vub,[],options);
xFitEllipsoid = EllipsoidGenerate(ellParamsFit,nThetaFit,nPhiFit);

% Plot the fit as a nice surface
xCoords = squeeze(xFitEllipsoid(1,:));
yCoords = squeeze(xFitEllipsoid(2,:));
zCoords = squeeze(xFitEllipsoid(3,:));
tri = delaunay(xCoords, yCoords, zCoords);
h = trisurf(tri, xCoords, yCoords, zCoords);
set(h,'FaceAlpha',0.25)
set(h,'EdgeColor',[0.5 0.5 0.5])
set(h,'FaceColor',[0.6 0.6 0.6]);
lighting phong;


end

function f = FitEllipseFunction(ellParams,x)

[A,Ainv,Q] = EllipsoidMatricesGenerate(ellParams);
vectorLengths = diag(x'*Q*x);
f = sqrt(mean((vectorLengths-1).^2));
end