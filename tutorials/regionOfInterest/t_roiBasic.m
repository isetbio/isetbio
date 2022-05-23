% Demo basic usage of the @regionOfInterest object
%
% Description:
%    Shows basic usage of the @regionOfInterest object
%    Here, we generate an ellipticalROI and a rectangular ROI. Then we
%    visualize them, and finally compute the indices of a set of random
%    points that either lie inside these ROIs or outside of them
%
% See Also:

% History:
%    11/16/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

regionOfInterest('help');
regionOfInterest('params');

% Instantiate an elliptical ROI
theROI = regionOfInterest(...
    'geometryStruct', struct(...
        'units', 'degs', ...
        'shape', 'ellipse', ...
        'center', [15.5, 1.5], ...
        'minorAxisDiameter', 10, ...
        'majorAxisDiameter', 5, ...
        'rotation', 13.0...
    ));
% Visualize it
theROI.visualize()

% Change the shape from 'ellipse' to 'rect', and also change the center, width & height
theROI.set('shape', 'rect', 'center', [-3 4], 'width', 10, 'height', 2);
theROI.visualize()

% Change the shape to 'line'
theROI.shape = 'line';

% Set the 'from', and 'to' params for the 'line' ROI
theROI.set('from', [-4 3], 'to', [2 2]);
% Visualize it
theROI.visualize();


theThickLineROI = regionOfInterest('shape', 'line', 'from', [-4 3], 'to', [2 2], 'thickness', 0.1);
theThickLineROI.visualize();

% Generate random points
randomPoints = bsxfun(@plus, [0 2], randn(2000,2)*1);

% Compute the indices of the random points that lie within theROI
indicesOfPointsInside = theROI.indicesOfPointsInside(randomPoints);

% Compute the indices of the random points that lie outside of theROI
indicesOfPointsOutside  = theROI.indicesOfPointsOutside(randomPoints);

% Compute the indices of the random points that lie around theROI
samplingPoints = 1000; % sample the perimeter using 1000 points
pointsPerSample = 30;  % up to 30 points for each sample along the perimeter
maxDistance = 0.5;     % up to 0.5 units aray from the closest point on the perimeter
indicesOfPointsNearROI = theROI.indicesOfPointsAround(randomPoints, pointsPerSample, samplingPoints, maxDistance);
indicesOfClosestPointsNearROI = theROI.indicesOfPointsAround(randomPoints, 1, samplingPoints, maxDistance);

% Compute the indices of the random points that lie around theThickLineROI
indicesOfPointsWithinThickLineROI = theThickLineROI.indicesOfPointsInside(randomPoints);

% Visualize
hFig = figure();

ax = subplot(3,1,1);
% Visualize the lineROI
theROI.visualize('figureHandle', hFig, 'axesHandle', ax);
% Plot all points
scatter(ax,randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);
% Label points  with a small distance from the lineROI
scatter(ax,randomPoints(indicesOfPointsNearROI,1), randomPoints(indicesOfPointsNearROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [1 0.3 0.0], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.1);
legend({'roi', 'all points', 'lots of points near ROI outline'}); 

ax = subplot(3,1,2);
% Visualize the lineROI
theROI.visualize('figureHandle', hFig, 'axesHandle', ax);
% Plot all points
scatter(ax,randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);
% Label points  with a large distance from the lineROI
scatter(ax,randomPoints(indicesOfClosestPointsNearROI,1), randomPoints(indicesOfClosestPointsNearROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [1 0.3 0.0], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.3);
legend({'roi', 'all points', 'closest points near ROI outline'}); 

ax = subplot(3,1,3);
% Visualize theThickLineROI
theThickLineROI.visualize('figureHandle', hFig, 'axesHandle', ax);
% Plot all points
scatter(ax,randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);
% Label points  with a large distance from the lineROI
scatter(ax,randomPoints(indicesOfPointsWithinThickLineROI,1), randomPoints(indicesOfPointsWithinThickLineROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [1 0.3 0.0], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.3);
legend({'roi', 'all points', 'points inside thick line ROI'}); 


