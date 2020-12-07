% EllipseTest
%
% Put an ellipse through a set of sample points.
%
% This illustrates use of FitEllipseQ, EllipsoidMatricesGenerate, and
% PointsOnEllipseQ.
%
% We can also try to fit ellipses using the EllipsoidFit routine, hacking
% in some ways to try force it to behave reasonably in the 2D plane.  This
% is coded here conditionally.  It doesn't currently work well, but might
% be fixed up with some work.

% History
%  11/25/20  dhb  Rewrote to illustrate a method that actually works

%% Clear
clear; close all;

%% Show 3D?
DoEllipsoidFit = false;

%% Generate some elliptical data
%
% Parameter format is reciprocol major axis length, reciprocol minor axis
% length, major axis angle (clockwise from x axis in degrees, in range -90
% to 90);
ellParamsTrue = [0.5 2 45];
nDataPoints = 10;
noiseSd = 0.1;
theDirsFit = UnitCircleGenerate(nDataPoints);
[~,~,QTrue] = EllipsoidMatricesGenerate(ellParamsTrue,'dimension',2);
theDataToFit = PointsOnEllipseQ(QTrue,theDirsFit);
theDataToFit = theDataToFit + normrnd(0,noiseSd,2,nDataPoints);

%% Fit using general ellipse fit routine
[ellParamsFit,fitA,fitAinv,fitQ,fitErr] = FitEllipseQ(theDataToFit);
[~,~,QFit] = EllipsoidMatricesGenerate(ellParamsFit,'dimension',2);
nPlotPoints = 200;
theDirsPlot = UnitCircleGenerate(nPlotPoints);
theEllTrue = PointsOnEllipseQ(QTrue,theDirsPlot);
theEllFit = PointsOnEllipseQ(QFit,theDirsPlot);

%% Fit using ellipsoid fit routines
if (DoEllipsoidFit)
    fitCenter = zeros(3,1);
    [~,~,Q3DFit,fitEllParams] = EllipsoidFit([theDataToFit ; zeros(1,nDataPoints)],[],false,true);
    theEllFit3D = PointsOnEllipsoidFind(Q3DFit,[theDirsPlot ; zeros(1,nPlotPoints)],fitCenter);
end

%% Plot
theColors = ['r' 'k' 'b' 'b' 'y' 'c'];
figure; clf; hold on;
plot(theDataToFit(1,:),theDataToFit(2,:),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
theLim = 2;
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
plot(theEllTrue(1,:),theEllTrue(2,:),'k.','MarkerSize',4,'MarkerFaceColor','k');
plot(theEllFit(1,:),theEllFit(2,:),'r.','MarkerSize',8,'MarkerFaceColor','r');
if (DoEllipsoidFit)
    plot(theEllFit3D(1,:),theEllFit3D(2,:),'g','MarkerSize',8,'MarkerFaceColor','g');
end
