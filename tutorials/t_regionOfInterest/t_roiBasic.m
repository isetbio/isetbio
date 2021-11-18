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

% BW Notes
%
%   1. Typing regionOfInterest should return a proper ROI (the default)
%   2. Typing regionOfInterest('help') should return a useful list or
%       some documentation.
%   3. I don’t think we should require the input to be a struct. I prefer the
%      key/val pairs. It seems to me that the ability to write like this is useful:
%
%            thisROI = regionOfInterest('shape','ellipse');
%            thisROI = regionOfInterest('shape','ellipse','minor axis diameter',1);
%
%      For a tutorial, this is an easier way to introduce the basic idea, and then 
%     buil uild up from  there.  Then we can adjust the properties
%  
%            thisROI.set('minor axis diameter',1);
%  
%      This method allows you to submit as a struct, too.  But requiring the 
%      struct() means you have to know the upper/lower case and 
%      the words precisely in order to make the struct.
%  
%


% Instantiate an elliptical ROI
opticDiskROI = regionOfInterest(...
    'geometryStruct', struct(...
        'units', 'degs', ...
        'shape', 'ellipse', ...
        'center', [15.5, 1.5], ...
        'minorAxisDiameter', 10, ...
        'majorAxisDiameter', 5, ...
        'rotation', 13.0...
    ));

% Instantiate a rectangular ROI
stimulusROI = regionOfInterest(...
    'geometryStruct', struct(...
        'units', 'degs', ...
        'shape', 'rect', ...
        'center', [14, 1], ...
        'width', 10, ...
        'height', 5, ...
        'rotation', 30.0...
));


% Instantiate a line ROI
lineROI1 = regionOfInterest(...
    'geometryStruct', struct(...
    'units', 'degs', ...
    'shape', 'line', ...
    'from', [10 -1], ...
    'to', [20 2] ...
));

% And another lineROI
lineROI2 = regionOfInterest(...
    'geometryStruct', struct(...
    'units', 'degs', ...
    'shape', 'line', ...
    'from', [6 9], ...
    'to', [17 -2] ...
));


% Generate random points
randomPoints = bsxfun(@plus, [14 2], randn(600,2)*3);


% Compute the indices of the random points that lie within the opticDiskROI
indicesOfPointsInsideOpticDiskROI = opticDiskROI.indicesOfPointsInside(randomPoints);

% Compute the indices of the random points that lie outside of the opticDiskROI
indicesOfPointsOutsideOpticDiskROI = opticDiskROI.indicesOfPointsOutside(randomPoints);

% Compute the indices of the random points that lie around the opticDiskROI outline
maxDistance = 0.5;
samplingPoints = 1000;
indicesOfPointsAroundOpticDiskROI = opticDiskROI.indicesOfPointsAround(randomPoints, maxDistance, samplingPoints);


% Compute the indices of the random points that lie within the stimulusROI
indicesOfPointsInsideStimulusROI = stimulusROI.indicesOfPointsInside(randomPoints);

% Compute the indices of the random points that lie outside of the stimulusROI
indicesOfPointsOutsideStimulusROI = stimulusROI.indicesOfPointsOutside(randomPoints);

% Compute the indices of the random points that lie outside of the stimulusROI
maxDistance = 0.1;
indicesOfPointsAroundStimulusROI = stimulusROI.indicesOfPointsAround(randomPoints, maxDistance, samplingPoints);


% Compute the indices of the random points that lie within the lineROI
maxDistance = 0.5;
indicesOfPointsAlongLineROI1 = lineROI1.indicesOfPointsAround(randomPoints, maxDistance, samplingPoints);
indicesOfPointsAlongLineROI2 = lineROI2.indicesOfPointsAround(randomPoints, maxDistance, samplingPoints);



% Visualize everything
sv = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 3, ...
       'heightMargin',  0.07, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.03);

hFig = figure(1); clf;
set(hFig, 'Position', [10 10 2000 1100]);

% Render the opticDiskROI and label the random points that lie within it
ax = subplot('Position', sv(1,1).v);
opticDiskROI.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [0.9 0.9 0]);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

% Plot points inside the opticDiskROI in red
scatter(ax, randomPoints(indicesOfPointsInsideOpticDiskROI,1), randomPoints(indicesOfPointsInsideOpticDiskROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.3 0.3], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.5);
title(ax, 'points inside ellipse ROI');

% Render the opticDiskROI and label the random points that lie outside
ax = subplot('Position', sv(1,2).v);
opticDiskROI.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [0.9 0.9 0]);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

% Plot points outside the opticDiskROI in red
scatter(ax, randomPoints(indicesOfPointsOutsideOpticDiskROI,1), randomPoints(indicesOfPointsOutsideOpticDiskROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.3 0.3], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.5);
title(ax, 'points outside ellipse ROI');


% Render the opticDiskROI and label the random points that lie outside
ax = subplot('Position', sv(1,3).v);
opticDiskROI.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [0.9 0.9 0]);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

% Plot points around the opticDiskROI in red
scatter(ax, randomPoints(indicesOfPointsAroundOpticDiskROI,1), randomPoints(indicesOfPointsAroundOpticDiskROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.3 0.3], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.5);
title(ax, 'points along ellipse ROI outline');





% Render the stimulusROI and label the random points that lie within it
ax = subplot('Position', sv(2,1).v);
stimulusROI.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [0.9 0.9 0]);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

% Plot points inside the stimulusROI in red
scatter(ax, randomPoints(indicesOfPointsInsideStimulusROI,1), randomPoints(indicesOfPointsInsideStimulusROI,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.3 0.3], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.5);
title(ax, 'points inside rect ROI');



% Render the 2 ROIs superimposed
ax = subplot('Position', sv(2,2).v);
opticDiskROI.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [1 1 0]);
stimulusROI.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [0 1 1]);
title(ax, 'rect ROI & ellipse ROI superimposed');

% Plot points inside both the stimulusROI and the opticDiskROI
idx = intersect(indicesOfPointsInsideStimulusROI,indicesOfPointsInsideOpticDiskROI);
scatter(ax, randomPoints(idx,1), randomPoints(idx,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.3 0.9], 'MarkerEdgeColor', [1 0 1]);
title(ax, 'points inside both rect ROI & ellipse ROI');

% Render the 2 lineROIs
ax = subplot('Position', sv(2,3).v);
lineROI1.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [1 0 0]);
lineROI2.visualize('figureHandle', hFig, 'axesHandle', ax, 'fillColor', [0 1 1]);

% Plot all points
scatter(ax, randomPoints(:,1), randomPoints(:,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

% Plot points along the lineROI in red
scatter(ax, randomPoints(indicesOfPointsAlongLineROI1,1), randomPoints(indicesOfPointsAlongLineROI1,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.9 0.3 0.3], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

scatter(ax, randomPoints(indicesOfPointsAlongLineROI2,1), randomPoints(indicesOfPointsAlongLineROI2,2), 64, ...
    'o', 'filled', 'MarkerFaceColor', [0.3 0.9 0.9], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

title(ax, 'points along line ROI1 and line ROI2');


