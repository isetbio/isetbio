function plotHexMosaic(obj, varargin)
% Visualize different aspects of the hex grid
%
% Name-Value options
%   showCorrespondingRectangularMosaicInstead - False
%   overlayNullSensorsPerfectHexMesh          - False
%   overlayPerfectHexMesh       - False
%   overlayConeDensityContour   - 'none'
%   coneDensityContourLevelStep - 5000
%
% Includes functions, some of which might get moved out
%   computeDensityMap
%   renderPatchArray
%   renderHexMesh
%
% NPC, ISETBIO TEAM, 2015

%% parse input
p = inputParser;
p.addParameter('showCorrespondingRectangularMosaicInstead', false, @islogical);
p.addParameter('overlayNullSensors', false, @islogical);
p.addParameter('overlayPerfectHexMesh', false, @islogical);
p.addParameter('overlayConeDensityContour', 'none', @ischar);
p.addParameter('coneDensityContourLevelStep', 5000, @isnumeric);
p.parse(varargin{:});

showCorrespondingRectangularMosaicInstead = p.Results.showCorrespondingRectangularMosaicInstead;
showNullSensors = p.Results.overlayNullSensors;
showPerfectHexMesh = p.Results.overlayPerfectHexMesh;
showConeDensityContour = p.Results.overlayConeDensityContour;
coneDensityContourLevelStep = p.Results.coneDensityContourLevelStep;

sampledHexMosaicXaxis = obj.patternSupport(1,:,1) + obj.center(1);
sampledHexMosaicYaxis = obj.patternSupport(:,1,2) + obj.center(2);

% Choose the radius of the aperture obj.pigment.pdWidth vs obj.pigment.width
radiusComesFrom = 'ConeAperture';
if (strcmp(radiusComesFrom, 'ConeAperture'))
    dx = obj.pigment.pdWidth;
else
    dx = obj.pigment.width;
end

pixelOutline.x = [-1 -1 1 1 -1]*dx/2;
pixelOutline.y = [-1 1 1 -1 -1]*dx/2;
originalPixelOutline.x = [-1 -1 1 1 -1]*dx/2.0;
originalPixelOutline.y = [-1 1 1 -1 -1]*dx/2.0;

iTheta = (0:15:360)/180*pi;
apertureOutline.x = dx/2.0 * cos(iTheta);
apertureOutline.y = dx/2.0 * sin(iTheta);

rectCoords = obj.coneLocsOriginatingRectGrid;
hexCoords = obj.coneLocsHexGrid;


% Clear axes
cla(obj.hdl.CurrentAxes, 'reset');

%% Do the display
switch showConeDensityContour
    case 'measured'
        [densityMap, densityMapSupportX, densityMapSupportY] = computeDensityMap(obj, 'from mosaic');
    case 'theoretical'
        [densityMap, densityMapSupportX, densityMapSupportY] = computeDensityMap(obj, 'from model');
    case 'none'
    otherwise
        error('coneMosaicHex.visualizeGrid: ''coneDensityContourOverlay'' must be set to one of the following: ''measured'', ''theoretical'', ''none''. ');
end

if (~showCorrespondingRectangularMosaicInstead)
    lineStyle = '-';
    if (showNullSensors)
        idx = find(obj.pattern==1);
        [iRows,iCols] = ind2sub(size(obj.pattern), idx);
        edgeColor = [0.4 0.4 0.4]; faceColor = 'none';
        renderPatchArray(pixelOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    end
    
    % L-cones
    idx = find(obj.pattern == 2);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = [1 0 0]; faceColor = [1.0 0.7 0.7];
    renderPatchArray(apertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    
    % M-cones
    idx = find(obj.pattern == 3);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = [0 0.7 0]; faceColor = [0.7 1.0 0.7];
    renderPatchArray(apertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    
    % S-cones
    idx = find(obj.pattern == 4);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = [0 0 1]; faceColor = [0.7 0.7 1.0];
    renderPatchArray(apertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    
    if (showPerfectHexMesh)
        % Superimpose hex mesh showing the locations of the perfect hex grid
        meshFaceColor = [0.8 0.8 0.8]; meshEdgeColor = [0.5 0.5 0.5]; meshFaceAlpha = 0.0; meshEdgeAlpha = 0.5; lineStyle = '-';
        renderHexMesh(hexCoords(:,1), hexCoords(:,2), meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle);
    end
else
    % Show the corresponding rectangular mosaic
    
    % The original rect sensors
    idx = find(obj.patternOriginatingRectGrid==2);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [1.0 0.7 0.7]; lineStyle = '-';
    renderPatchArray(originalPixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
    
    idx = find(obj.patternOriginatingRectGrid==3);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [0.7 1.0 0.7];
    renderPatchArray(originalPixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
    
    idx = find(obj.patternOriginatingRectGrid==4);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [0.7 0.7 1.0];
    renderPatchArray(originalPixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
end

if (~strcmp(showConeDensityContour, 'none'))
    contourLevels = coneDensityContourLevelStep: coneDensityContourLevelStep: 250000;
    [cH, hH] = contour(densityMapSupportX, densityMapSupportY, densityMap, contourLevels, 'LineColor', 'k', 'LineWidth', 3.0, 'ShowText', 'on', 'LabelSpacing', 500);
    clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, 'Color', [0 0 0])
    set(gca, 'CLim', [10000 250000]);
end

%% Arrange axis and fonts

hold off
axis 'equal'; axis 'xy'
xTicks = [sampledHexMosaicXaxis(1) obj.center(1) sampledHexMosaicXaxis(end)];
yTicks = [sampledHexMosaicYaxis(1) obj.center(2) sampledHexMosaicYaxis(end)];
xTickLabels = sprintf('%2.0f um\n', xTicks*1e6);
yTickLabels = sprintf('%2.0f um\n', yTicks*1e6);
set(gca, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
set(gca, 'FontSize', 16, 'XColor', [0.1 0.2 0.9], 'YColor', [0.1 0.2 0.9], 'LineWidth', 1.0);
box on; grid off;
set(gca, 'XLim', [sampledHexMosaicXaxis(1)-dx sampledHexMosaicXaxis(end)+dx]);
set(gca, 'YLim', [sampledHexMosaicYaxis(1)-dx sampledHexMosaicYaxis(end)+dx]);

drawnow;

end

%% Separate function?  Utility?
function [densityMap, densityMapSupportX, densityMapSupportY] = computeDensityMap(obj, computeConeDensityMap)

deltaX = 3*obj.lambdaMin*1e-6;
areaHalfWidth = deltaX*4;

mosaicRangeX = obj.center(1) + obj.width/2*[-1 1]  + [-obj.lambdaMin obj.lambdaMin]*1e-6;
mosaicRangeY = obj.center(2) + obj.height/2*[-1 1] + [-obj.lambdaMin obj.lambdaMin]*1e-6;

gridXPos = (mosaicRangeX(1)+areaHalfWidth):deltaX:(mosaicRangeX(2)-areaHalfWidth);
gridYPos = (mosaicRangeY(1)+areaHalfWidth):deltaX:(mosaicRangeY(2)-areaHalfWidth);
measurementAreaInMM2 = (2*areaHalfWidth*1e3)^2;
densityMap = zeros(numel(gridYPos), numel(gridXPos));

if (strcmp(computeConeDensityMap, 'from mosaic'))
    for iYpos = 1:numel(gridYPos)
        for iXpos = 1:numel(gridXPos)
            xo = gridXPos(iXpos);
            yo = gridYPos(iYpos);
            conesWithin = numel(find( ...
                obj.coneLocsHexGrid(:,1) >= xo-areaHalfWidth  & ...
                obj.coneLocsHexGrid(:,1) <= xo+areaHalfWidth  & ...
                obj.coneLocsHexGrid(:,2) >= yo-areaHalfWidth  & ...
                obj.coneLocsHexGrid(:,2) <= yo+areaHalfWidth ));
            densityMap(iYpos,iXpos) = conesWithin / measurementAreaInMM2;
        end
    end
else
    [X,Y] = meshgrid(gridXPos, gridYPos);
    eccInMeters = sqrt(X.^2 + Y.^2);
    ang = atan2(Y, X)/pi*180;
    [~, ~, densities] = coneSize(eccInMeters(:),ang(:));
    densityMap = reshape(densities, size(X));
end

gridXPos = mosaicRangeX(1) + (gridXPos - min(gridXPos))/(max(gridXPos) - min(gridXPos))*(mosaicRangeX(2)-mosaicRangeX(1));
gridYPos = mosaicRangeY(1) + (gridYPos - min(gridYPos))/(max(gridYPos) - min(gridYPos))*(mosaicRangeY(2)-mosaicRangeY(1));
[densityMapSupportX, densityMapSupportY] = meshgrid(gridXPos, gridYPos);
end

%% Maybe put in utility directory
function renderPatchArray(pixelOutline, xCoords, yCoords, edgeColor, faceColor, lineStyle)

verticesNum = numel(pixelOutline.x);
x = zeros(verticesNum, numel(xCoords));
y = zeros(verticesNum, numel(xCoords));

for vertexIndex = 1:verticesNum
    x(vertexIndex, :) = pixelOutline.x(vertexIndex) + xCoords;
    y(vertexIndex, :) = pixelOutline.y(vertexIndex) + yCoords;
end
patch(x, y, [0 0 0], 'EdgeColor', edgeColor, 'FaceColor', faceColor, 'LineWidth', 1.0, 'LineStyle', lineStyle);
end

%% Separate function??
function renderHexMesh(xHex, yHex, meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle)
x = []; y = [];
triangleConeIndices = delaunayn([xHex(:), yHex(:)]);
for triangleIndex = 1:size(triangleConeIndices,1)
    coneIndices = triangleConeIndices(triangleIndex, :);
    xCoords = xHex(coneIndices);
    yCoords = yHex(coneIndices);
    for k = 1:numel(coneIndices)
        x = cat(2, x, xCoords);
        y = cat(2, y, yCoords);
    end
end
patch(x, y, [0 0 0], 'EdgeColor', meshEdgeColor, 'EdgeAlpha', meshEdgeAlpha, 'FaceAlpha', meshFaceAlpha, 'FaceColor', meshFaceColor, 'LineWidth', 1.5, 'LineStyle', lineStyle);
end

