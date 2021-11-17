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

% Instantiate an elliptical ROI
opticDiskROI = regionOfInterest(...
    struct(...
        'units', 'degs', ...
        'shape', 'ellipse', ...
        'center', [15.5, 1.5], ...
        'minorAxisDiameter', 6.3, ...
        'majorAxisDiameter', 6.7, ...
        'rotation', 13.0...
    ));

% Instantiate a rectangular ROI
stimulusROI = regionOfInterest(...
    struct(...
    'units', 'degs', ...
    'shape', 'rect', ...
    'center', [14, 2], ...
    'width', 10, ...
    'height', 5, ...
    'rotation', 30.0...
));

% Generate random points
randomPoints = bsxfun(@plus, [14 2], rand(100,2)*5);

% Compute the indices of the random points that lie within the stimulusROI
indicesOfPointsInsideStimulusROI = stimulusROI.indicesOfPointsInside(randomPoints);

% Compute the indices of the random points that lie within the opticDiskROI
indicesOfPointsInsideOpticDiskROI = opticDiskROI.indicesOfPointsInside(randomPoints);

% Compute the indices of the random points that lie outside of the stimulusROI
indicesOfPointsOutsideStimulusROI = stimulusROI.indicesOfPointsOutside(randomPoints);

% Compute the indices of the random points that lie outside of the opticDiskROI
indicesOfPointsOutsideOpticDiskROI = opticDiskROI.indicesOfPointsOutside(randomPoints);

% Visualize everything
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1100 1000]);

% Render the opticDiskROI and label the random points that lie within it
ax = subplot(2,3,1);
opticDiskROI.visualize('figureHandle', hFig, 'axesHandle', ax);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0]);

% Plot points inside the opticDiskROI in red
scatter(ax, randomPoints(indicesOfPointsInsideOpticDiskROI,1), randomPoints(indicesOfPointsInsideOpticDiskROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.3 0.3], 'MarkerEdgeColor', [1 0 0]);
title(ax, 'points inside ellipse ROI');

% Render the opticDiskROI and label the random points that lie outside
ax = subplot(2,3,4);
opticDiskROI.visualize('figureHandle', hFig, 'axesHandle', ax);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0]);

% Plot points outside the opticDiskROI in red
scatter(ax, randomPoints(indicesOfPointsOutsideOpticDiskROI,1), randomPoints(indicesOfPointsOutsideOpticDiskROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.3 0.3], 'MarkerEdgeColor', [1 0 0]);
title(ax, 'points outside ellipse ROI');



% Render the stimulusROI and label the random points that lie within it
ax = subplot(2,3,2);
stimulusROI.visualize('figureHandle', hFig, 'axesHandle', ax);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0]);

% Plot points inside the stimulusROI in blue
scatter(ax, randomPoints(indicesOfPointsInsideStimulusROI,1), randomPoints(indicesOfPointsInsideStimulusROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.9], 'MarkerEdgeColor', [0 0 1]);
title(ax, 'points inside rect ROI');


% Render the stimulusROI and label the random points that lie outside
ax = subplot(2,3,5);
stimulusROI.visualize('figureHandle', hFig, 'axesHandle', ax);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0]);

% Plot points inside the stimulusROI in blue
scatter(ax, randomPoints(indicesOfPointsOutsideStimulusROI,1), randomPoints(indicesOfPointsOutsideStimulusROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.9], 'MarkerEdgeColor', [0 0 1]);
title(ax, 'points outside rect ROI');



% Render the 2 ROIs superimposed
ax = subplot(2,3,3);
opticDiskROI.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [1 1 0]);
stimulusROI.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [0 1 1]);

