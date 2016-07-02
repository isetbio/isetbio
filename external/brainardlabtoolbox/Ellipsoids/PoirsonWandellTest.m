% PoirsonWandellTest
%
% Illustrate routines that generate data according to the pattern-color
% separable model given in Poirson-Wandell
%    Poirson AB, Wandell BA. 1996. Pattern-color separable pathways predict
%    sensitivity to simple colored patterns Vision Res 36: 515-26.
%
% 6/29/16  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Get the pattern-color model parameters from the tables in the paper
conditionStr = 'HT,cc';
theSf = 2;
[A,Ainv,Q,theBgLMS] = PoirsonWandellEllipsoidParameters(conditionStr,theSf);

%% Generate and plot the ellipsoid in 3D
%
% The ellipsoid comes back in cone difference coordinates,
% not contrast.  Divide by the background to get cone contrast.
nThetaEllipsoid = 20;
nPhiEllipsoid = 20;
xSphere = UnitSphereGenerate(nThetaEllipsoid,nPhiEllipsoid);
xEllipsoid = Ainv*xSphere;
xEllipsoid = bsxfun(@times,xEllipsoid,1./theBgLMS);

% Plot the fit as a nice surface
xCoords = squeeze(xEllipsoid(1,:));
yCoords = squeeze(xEllipsoid(2,:));
zCoords = squeeze(xEllipsoid(3,:));
tri = delaunay(xCoords, yCoords, zCoords);
figure; clf; hold on
h = trisurf(tri, xCoords, yCoords, zCoords);
set(h,'FaceAlpha',0.25)
set(h,'EdgeColor',[0.5 0.5 0.5])
set(h,'FaceColor',[0.6 0.6 0.6]);
lighting phong;
xlabel('L contrast'); ylabel('M contrast'); zlabel('S contrast'); title('Poirson Wandell Ellipsoid');
%xlim([-0.15 0.15]); ylim([-0.15 0.15]); zlim([-0.25 0.25]);
axis('square');

%% Generate 2D slices 
%
% This plot reproduces Figure 3 in the Poirson & Wandell paper, and
% gives us reason to believe we understand the model and how its parameters
% are specified in the paper.
nThetaEllipse = 200;
xCircle = UnitCircleGenerate(nThetaEllipse);
xCirclePlane = [xCircle(1,:) ; xCircle(2,:) ; zeros(size(xCircle(1,:)))];
xEllipsoidPlane = PointsOnEllipsoidFind(Q,xCirclePlane);
xEllipsoidPlane = bsxfun(@times,xEllipsoidPlane,1./theBgLMS);
figure; clf;
subplot(1,3,1); hold on
plot(xEllipsoidPlane(1,:),xEllipsoidPlane(2,:),'r','LineWidth',3);
plot([-0.015 0.015],[0 0],'k:','LineWidth',2);
plot([0 0],[-0.015 0.015],'k:','LineWidth',2);
xlabel('L contrast'); ylabel('M contrast'); title('Poirson Wandell Ellipsoid');
xlim([-0.015 0.015]); ylim([-0.015 0.015]);
set(gca,'XTick',[-0.015 -0.010 -0.005 0 0.005 0.010 0.015]);
set(gca,'XTickLabel',{'-0.015' '' '' '' '' '' '0.015'});
set(gca,'YTick',[-0.015 -0.010 -0.005 0 0.005 0.010 0.015]);
set(gca,'YTickLabel',{'-0.015' '' '' '' '' '' '0.015'});
axis('square');

xCircle = UnitCircleGenerate(nThetaEllipse);
xCirclePlane = [xCircle(1,:) ; zeros(size(xCircle(1,:))) ; xCircle(2,:) ];
xEllipsoidPlane = PointsOnEllipsoidFind(Q,xCirclePlane);
xEllipsoidPlane = bsxfun(@times,xEllipsoidPlane,1./theBgLMS);
subplot(1,3,2); hold on
plot(xEllipsoidPlane(1,:),xEllipsoidPlane(3,:),'r','LineWidth',3);
plot([-0.02 0.02],[0 0],'k:','LineWidth',2);
plot([0 0],[-0.04 0.04],'k:','LineWidth',2);
xlabel('L contrast'); ylabel('S contrast'); title('Poirson Wandell Ellipsoid');
xlim([-0.02 0.02]); ylim([-0.04 0.04]);
axis('square');

xCirclePlane = [zeros(size(xCircle(1,:))) ; xCircle(1,:) ; xCircle(2,:)];
xEllipsoidPlane = PointsOnEllipsoidFind(Q,xCirclePlane);
xEllipsoidPlane = bsxfun(@times,xEllipsoidPlane,1./theBgLMS);
subplot(1,3,3); hold on
plot(xEllipsoidPlane(2,:),xEllipsoidPlane(3,:),'r','LineWidth',3);
plot([-0.02 0.02],[0 0],'k:','LineWidth',2);
plot([0 0],[-0.04 0.04],'k:','LineWidth',2);
xlabel('M contrast'); ylabel('S contrast'); title('Poirson Wandell Ellipsoid');
xlim([-0.02 0.02]); ylim([-0.04 0.04]);
axis('square');

