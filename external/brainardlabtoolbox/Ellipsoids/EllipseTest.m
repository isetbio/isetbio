% EllipseTest
%
% Put an ellipse through a set of sample points.
%
% Then map the data back to space where it is on a unit sphere
% and plot the xy plane in that space.
%
% Finally, make a histogram of the data vector lengths in the spherical
% space.

%% Clear
clear; close all;

%% Fit center offset?
fitCenterOffset = false;

%% Load in the data
whichDataSet = 4;
theRawData = load('dataForEllipseTest.mat');
theData0 = theRawData.dataForEllipse(theRawData.ind==whichDataSet,:)';

%% Embed in 3D so we can use ellipsoid code
theData0(3,:) = 0*ones(1,size(theData0,2));
theData1 = 0.5*theData0;
theData1(3,:) = 0.5;
theData2 = 0.5*theData0;
theData2(3,:) = -0.5;
theDataToFit = [theData0 theData1 theData2];

%% Deal with offset
%
% Our ellipse fitting routine tested here does not allow fitting of the
% center offset. Although we should fix this, for now deal with by
% subtracting the mean of the data from the data
if (~fitCenterOffset)
    theDataToFit = theDataToFit-mean(theDataToFit,1);
end

%% Fit using general routine
[fitA,fitAinv,fitQ,fitEllParams] = EllipsoidFit(theDataToFit,[],fitCenterOffset,true);

%% Grab ellipsoid center from params.  
if (fitCenterOffset)
    fitXCenter = fitEllParams(7:9);
else
    fitXCenter = zeros(3,1);
end

%% Get the LM plane ellipse from the fit
nThetaEllipse = 200;
circleIn2D = UnitCircleGenerate(nThetaEllipse);
circleInLMPlane = [circleIn2D(1,:) ; circleIn2D(2,:) ; zeros(size(circleIn2D(1,:)))];
fitEllipse = PointsOnEllipsoidFind(fitQ,circleInLMPlane,fitXCenter);

%% Plot
figure; clf; hold on
plot(theData0(1,:),theData0(2,:),'b.','MarkerSize',8);
plot(fitEllipse(1,:),fitEllipse(2,:),'r');
title('Ellipse fit to data');
xlabel('Coordinate 1');
ylabel('Coordinate 2');

%% Subtract off ellipsoid center location from the data
for ii = 1:3
    theData0Centered(ii,:) = theData0(ii,:)-fitXCenter(ii);
    ellipseFitCentered(ii,:) = fitEllipse(ii,:)-fitXCenter(ii);
end

%% Map points to unit sphere, and plot the xy plane 
%
% It is not guaranteed that fit data in the xy plane map back into
% the xy plane in the spherized space, but in practice it appears close to 
% true;
theDataSphere = fitA*theData0Centered;
theFitSphere = fitA*ellipseFitCentered;
figure; clf; hold on
plot(theDataSphere(1,:),theDataSphere(2,:),'b.','MarkerSize',8);
plot(theFitSphere(1,:),theFitSphere(2,:),'r');
title('Ellipse fit to data');
xlabel('Coordinate 1');
ylabel('Coordinate 2');
axis('square'); axis([-1.5 1.5 -1.5 1.5]);

%% Make a map of the vector lengths of the spherized data
figure; clf;
fitLengths = sqrt(diag(theDataSphere'*theDataSphere));
hist(fitLengths); xlim([0 2]);
xlabel('Data Length in Spherical Space');
ylabel('Count');

