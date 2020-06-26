function plotHexMosaic(obj, varargin)
% Visualize the hex grid
%
% Syntax:
%   plotHexMosaic(obj, [varargin])
%
% Description:
%   Using the key/value pairs you can visualize different aspects of
%   the hexagonal cone mosaic.
%
% Inputs:
%    obj                                       - The cone mosaic hex object
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%	 showCorrespondingRectangularMosaicInstead - Boolean, Default False
%    overlayNullSensorsPerfectHexMesh          - Boolean, Default False
%    overlayPerfectHexMesh                     - Boolean, Default False
%    overlayConeDensityContour                 - String, Default 'none'
%    coneDensityContourLevelStep               - Integer, Default 5000
%
% Note:
%    * TODO - Assign someone to fix the CurrentAxes problem that DHB
%      describes below.
%

% History:
%    xx/xx/15  NPC  ISETBIO TEAM, 2015
%    02/20/18  jnm  Formatting
%    03/08/19  NPC  Fixed subroutine calling issues
%    5/26/20   NPC  Fixed iRow issue, which was causing the mosaic plotting Y-coord flip 
%                   (rows grow top -> bottom, whereas Y-coords grow bottom -> top)

%% parse input
p = inputParser;
p.KeepUnmatched = true;   % This allows unrecognized parameters

% Defaulting this to true until we solve the rendering speed problem
p.addParameter('showCorrespondingRectangularMosaicInstead', ...
    false, @islogical);
p.addParameter('overlayNullSensors', false, @islogical);
p.addParameter('overlayPerfectHexMesh', false, @islogical);
p.addParameter('overlayConeDensityContour', 'none', @ischar);
p.addParameter('coneDensityContourLevelStep', 5000, @isnumeric);

% The commented out line fails because referencing obj.hdl.CurrentAxes
% generates a "Struct contents reference from a non-struct array object"
% error. I commented out and made the default empty. This, at least, 
% gets the tutorials to run. DHB
%p.addParameter('hf', obj.hdl.CurrentAxes, @isgraphics);
p.addParameter('hf', [], @isgraphics);

p.parse(varargin{:});

% From the coneMosaicWindow, this is normally the rectangular version.
% From the Plot | Mosaic | Cone Mosaic we show the hex data
showCorrespondingRectangularMosaicInstead = ...
    p.Results.showCorrespondingRectangularMosaicInstead;
showNullSensors = p.Results.overlayNullSensors;
showPerfectHexMesh = p.Results.overlayPerfectHexMesh;
showConeDensityContour = p.Results.overlayConeDensityContour;
coneDensityContourLevelStep = p.Results.coneDensityContourLevelStep;

% Sometimes this is the main window and sometimes this is a new window
% passed in.
thisAxes = p.Results.hf;

sampledHexMosaicXaxis = obj.patternSupport(1, :, 1) + obj.center(1);
sampledHexMosaicYaxis = obj.patternSupport(:, 1, 2) + obj.center(2);

% Choose radius of the aperture obj.pigment.pdWidth vs obj.pigment.width
radiusComesFrom = 'ConeAperture';
if (strcmp(radiusComesFrom, 'ConeAperture'))
    dx = obj.pigment.pdWidth;
else
    dx = obj.pigment.width;
end

pixelOutline.x = [-1 -1 1 1 -1] * dx / 2;
pixelOutline.y = [-1 1 1 -1 -1] * dx / 2;
originalPixelOutline.x = [-1 -1 1 1 -1] * dx / 2.0;
originalPixelOutline.y = [-1 1 1 -1 -1] * dx / 2.0;

dAngle = 30;
iTheta = (0:dAngle:360 - dAngle) / 180 * pi;
apertureOutline.x = dx / 2.0 * cos(iTheta);
apertureOutline.y = dx / 2.0 * sin(iTheta);

rectCoords = obj.coneLocsOriginatingRectGrid;
hexCoords = obj.coneLocsHexGrid;

% Clear axes, which sometimes is a figure handle that has one axis
cla(thisAxes, 'reset');

%% Do the display
switch showConeDensityContour
    case 'measured'
        [densityMap, densityMapSupportX, densityMapSupportY] = ...
            obj.computeDensityMap('from mosaic');
    case 'theoretical'
        [densityMap, densityMapSupportX, densityMapSupportY] = ...
            obj.computeDensityMap('from model');
    case 'none'
    otherwise
        error(['coneMosaicHex.visualizeGrid: '...
            'coneDensityContourOverlay'...
            ' must be set to one of the following: ''measured'', '...
            'theoretical'', ''none''. ']);
end

if (~showCorrespondingRectangularMosaicInstead)
    % Show Hex mosaic (don't show the rectangular mosaic)
    % We are not using this code right now. We are defaulting to rect
    % because of speed of the rendering. We are going to make a pull down
    % that renders the hex mosaic in a separate window for the mean time.
    % When we the solve the problem we will use that code in here.
    lineStyle = '-';
    lineWidth = 1.0;
    if (showNullSensors)
        idx = find(obj.pattern == 1);
        [iRows, iCols] = ind2sub(size(obj.pattern), idx);
        edgeColor = [0.4 0.4 0.4];
        faceColor = 'none';
        renderPatchArray(gca, pixelOutline, sampledHexMosaicXaxis(iCols), ...
            sampledHexMosaicYaxis(end-iRows+1), edgeColor, faceColor, lineStyle, lineWidth);
    end

    % L-cones - The 'finds' take a long time, 3-times. Let's see if we
    % can't speed it up.
    idx = find(obj.pattern == 2);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = [1 0 0];
    faceColor = [1.0 0.7 0.7];
    obj.renderPatchArray(gca, apertureOutline, sampledHexMosaicXaxis(iCols), ...
        sampledHexMosaicYaxis(end-iRows+1), edgeColor, faceColor, lineStyle, lineWidth);

    % M-cones
    idx = find(obj.pattern == 3);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = [0 0.7 0];
    faceColor = [0.7 1.0 0.7];
    obj.renderPatchArray(gca, apertureOutline, sampledHexMosaicXaxis(iCols), ...
        sampledHexMosaicYaxis(end-iRows+1), edgeColor, faceColor, lineStyle, lineWidth);

    % S-cones
    idx = find(obj.pattern == 4);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = [0 0 1];
    faceColor = [0.7 0.7 1.0];
    obj.renderPatchArray(gca, apertureOutline, sampledHexMosaicXaxis(iCols), ...
        sampledHexMosaicYaxis(end-iRows+1), edgeColor, faceColor, lineStyle, lineWidth);

    if (showPerfectHexMesh)
        % Superimpose hex mesh showing locations of the perfect hex grid
        meshFaceColor = [0.8 0.8 0.8];
        meshEdgeColor = [0.5 0.5 0.5];
        meshFaceAlpha = 0.0;
        meshEdgeAlpha = 0.5;
        lineStyle = '-';
        obj.renderHexMesh(gca, hexCoords(:, 1), hexCoords(:, 2), meshEdgeColor, ...
            meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle);
    end
else
    % Show the corresponding rectangular mosaic
    % This is the code we use for now

    % The original rect sensors
    idx = find(obj.patternOriginatingRectGrid==2);
    % [iRows, iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3];
    faceColor = [1.0 0.7 0.7];
    lineStyle = '-';
    lineWidth = 1.0;
    renderPatchArray(gca, originalPixelOutline, rectCoords(idx, 1), ...
        rectCoords(idx, 2), edgeColor, faceColor, lineStyle, lineWidth);

    idx = find(obj.patternOriginatingRectGrid==3);
    %[iRows, iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3];
    faceColor = [0.7 1.0 0.7];
    renderPatchArray(gca, originalPixelOutline, rectCoords(idx, 1), ...
        rectCoords(idx, 2), edgeColor, faceColor, lineStyle, lineWidth);

    idx = find(obj.patternOriginatingRectGrid==4);
    %[iRows, iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3];
    faceColor = [0.7 0.7 1.0];
    renderPatchArray(gca, originalPixelOutline, rectCoords(idx, 1), ...
        rectCoords(idx, 2), edgeColor, faceColor, lineStyle, lineWidth);
end

if (~strcmp(showConeDensityContour, 'none'))
    contourLevels = coneDensityContourLevelStep: ...
        coneDensityContourLevelStep:250000;
    [cH, hH] = contour(densityMapSupportX, densityMapSupportY, ...
        densityMap, contourLevels, 'LineColor', 'k', 'LineWidth', 3.0, ...
        'ShowText', 'on', 'LabelSpacing', 500);
    clabel(cH, hH, 'FontWeight', 'bold', 'FontSize', 16, 'Color', [0 0 0])
    set(gca, 'CLim', [10000 250000]);
end

%% Arrange axis and fonts
hold off;
axis 'equal';
axis 'xy'
xTicks = [sampledHexMosaicXaxis(1), obj.center(1), ...
    sampledHexMosaicXaxis(end)];
yTicks = [sampledHexMosaicYaxis(1), obj.center(2), ...
    sampledHexMosaicYaxis(end)];
xTickLabels = sprintf('%2.0f um\n', xTicks * 1e6);
yTickLabels = sprintf('%2.0f um\n', yTicks * 1e6);
set(gca, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', xTickLabels, ...
    'YTickLabel', yTickLabels);
set(gca, 'FontSize', 16, ...
    'LineWidth', 1.0); 
    % 'XColor', [0.1 0.2 0.9], ...
    % 'YColor', [0.1 0.2 0.9]);
box on;
grid off;
set(gca, 'XLim', [sampledHexMosaicXaxis(1) - dx, ...
    sampledHexMosaicXaxis(end) + dx]);
set(gca, 'YLim', [sampledHexMosaicYaxis(1) - dx, ...
    sampledHexMosaicYaxis(end) + dx]);

drawnow;

end

%% Separate function?  Utility?
function [densityMap, densityMapSupportX, densityMapSupportY] = ...
    computeDensityMap(obj, computeConeDensityMap)
% Compute the cone mosaic hex density map
%
% Syntax:
%   [densityMap, densityMapSupportX, densityMapSupportY] = ...
%       computeDensityMap(obj, computeConeDensityMap)
%
% Description:
%    Compute the density map from the cone mosaic hex object and a provided
%    cone density map computation
%
% Inputs:
%    obj                   - The cone mosaic hex object
%    computeConeDensityMap - The computed cone density map
%
% Outputs:
%    densityMap            - The calculated density map
%    densityMapSupportX    - The X-axis support for the density map
%    densityMapSupportY    - The Y-axis support for the density map
%
% Optional key/value pairs:
%    None.
%

deltaX = 3 * obj.lambdaMin * 1e-6;
areaHalfWidth = deltaX * 4;

mosaicRangeX = obj.center(1) + obj.width / 2 * [-1 1]  + ...
    [-obj.lambdaMin obj.lambdaMin] * 1e-6;
mosaicRangeY = obj.center(2) + obj.height / 2 * [-1 1] + ...
    [-obj.lambdaMin obj.lambdaMin] * 1e-6;

gridXPos = (mosaicRangeX(1) + areaHalfWidth):deltaX:(mosaicRangeX(2) - ...
    areaHalfWidth);
gridYPos = (mosaicRangeY(1) + areaHalfWidth):deltaX:(mosaicRangeY(2) - ...
    areaHalfWidth);
measurementAreaInMM2 = (2 * areaHalfWidth * 1e3) ^ 2;
densityMap = zeros(numel(gridYPos), numel(gridXPos));

if (strcmp(computeConeDensityMap, 'from mosaic'))
    for iYpos = 1:numel(gridYPos)
        for iXpos = 1:numel(gridXPos)
            xo = gridXPos(iXpos);
            yo = gridYPos(iYpos);
            conesWithin = numel(find( ...
                obj.coneLocsHexGrid(:, 1) >= xo-areaHalfWidth  & ...
                obj.coneLocsHexGrid(:, 1) <= xo+areaHalfWidth  & ...
                obj.coneLocsHexGrid(:, 2) >= yo-areaHalfWidth  & ...
                obj.coneLocsHexGrid(:, 2) <= yo+areaHalfWidth ));
            densityMap(iYpos, iXpos) = conesWithin / measurementAreaInMM2;
        end
    end
else
    [X, Y] = meshgrid(gridXPos, gridYPos);
    eccInMeters = sqrt(X .^ 2 + Y .^ 2);
    ang = atan2(Y, X) / pi * 180;
    [~, ~, densities] = coneSizeReadData(...
        'eccentricity', eccInMeters(:), 'angle', ang(:));
    densityMap = reshape(densities, size(X));
end

gridXPos = mosaicRangeX(1) + ...
    (gridXPos - min(gridXPos)) / (max(gridXPos) - min(gridXPos)) * ...
    (mosaicRangeX(2) - mosaicRangeX(1));
gridYPos = mosaicRangeY(1) + ...
    (gridYPos - min(gridYPos)) / (max(gridYPos) - min(gridYPos)) * ...
    (mosaicRangeY(2) - mosaicRangeY(1));
[densityMapSupportX, densityMapSupportY] = meshgrid(gridXPos, gridYPos);
end
