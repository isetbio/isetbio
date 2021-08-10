% Demo basic usage of the @mRGCMosaic object
%
% Description:
%    Shows basic usage of the new @mRGCmosaic.
%    Here, we generate a midget RGC mosaic object and visualize it along
%    with its input @cMosaic
%
% See Also:
%
%

% History:
%    04/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

%% Initialize
ieInit;
clear;
close all;

%% Generate the 10 deg foveal mosaic
onMRGCmosaic = mRGCMosaic(...
    'whichEye', 'left eye', ...
    'sizeDegs', 54.0*[1 1], ...     
    'eccentricityDegs', [0 0] ...
    );

save('ONMRGCmosaic54degs.mat', 'onMRGCmosaic', '-v7.3');

return;

rfPos = onMRGCmosaic.inputConeMosaic.coneRFpositionsDegs;
rfSpacings = onMRGCmosaic.inputConeMosaic.coneRFspacingsDegs;
sampledPositions{1} = -5:0.2:5;
sampledPositions{2} = -5:0.2:5;
densityMapCones = cMosaic.densityMap(rfPos, rfSpacings, sampledPositions);


rfPos = onMRGCmosaic.rgcRFpositionsDegs;
rfSpacings = onMRGCmosaic.rgcRFspacingsDegs;
densityMapRGCs = cMosaic.densityMap(rfPos, rfSpacings, sampledPositions);


figure(10);
subplot(1,3,1);
contourLabelSpacing = 4000;
densityContourLevels = 10;
[cH, hH] = contour(sampledPositions{1} , sampledPositions{2} , ...
                densityMapCones, densityContourLevels, 'LineColor', 'k', 'LineWidth', 2.0, ...
                'ShowText', 'on', 'LabelSpacing', contourLabelSpacing);
clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, ...
                'Color', [0 0 0], 'BackgroundColor', 'none');
axis 'equal'      
            
subplot(1,3,2);

[cH, hH] = contour(sampledPositions{1} , sampledPositions{2} , ...
                densityMapRGCs, densityContourLevels, 'LineColor', 'k', 'LineWidth', 2.0, ...
                'ShowText', 'on', 'LabelSpacing', contourLabelSpacing);
clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, ...
                'Color', [0 0 0], 'BackgroundColor', 'none');
axis 'equal'            

subplot(1,3,3);
densityContourLevels = 0.5:0.1:1.0;
[cH, hH] = contour(sampledPositions{1} , sampledPositions{2} , ...
                densityMapCones./densityMapRGCs, densityContourLevels, 'LineColor', 'k', 'LineWidth', 2.0, ...
                'ShowText', 'on', 'LabelSpacing', contourLabelSpacing);
clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, ...
                'Color', [0 0 0], 'BackgroundColor', 'none');
axis 'equal'                  
pause
xyLims = [-1 1 -1 1];
onMRGCmosaic.visualize(...
    'domainVisualizationLimits', xyLims, ...
    'covisualizeInputConeMosaic', true);

