function hFig = visualizeGrid(obj, varargin)
% Visualize different aspects of the hex grid
%
% Name-Value options
%   axesHandle - Empty - Axes handle to draw on
%   generateNewFigure  - False
%   panelPosition      - [1 1]
%   showCorrespondingRectangularMosaicInstead - False
%   overlayNullSensorsPerfectHexMesh          - False
%   overlayPerfectHexMesh       - False
%   overlayConeDensityContour   - 'none'
%   coneDensityContourLevelStep - 5000
%
% Includes functions, some of which might get moved out
%   renderPatchArray
%   renderHexMesh
%
% NPC, ISETBIO TEAM, 2015

%% parse input
p = inputParser;
p.addParameter('generateNewFigure', false, @islogical);
p.addParameter('panelPosition', [1 1]);
p.addParameter('axesHandle', []);
p.addParameter('showCorrespondingRectangularMosaicInstead', false, @islogical);
p.addParameter('visualizedConeAperture', 'lightCollectingArea', @(x)ismember(x, {'lightCollectingArea', 'geometricArea'}));
p.addParameter('overlayNullSensors', false, @islogical);
p.addParameter('overlayPerfectHexMesh', false, @islogical);
p.addParameter('overlayConeDensityContour', 'none', @ischar);
p.addParameter('coneDensityContourLevelStep', 5000, @isnumeric);
p.parse(varargin{:});

showCorrespondingRectangularMosaicInstead = p.Results.showCorrespondingRectangularMosaicInstead;
showNullSensors = p.Results.overlayNullSensors;
showPerfectHexMesh = p.Results.overlayPerfectHexMesh;
showConeDensityContour = p.Results.overlayConeDensityContour;
generateNewFigure = p.Results.generateNewFigure;
panelPosition = p.Results.panelPosition;
coneDensityContourLevelStep = p.Results.coneDensityContourLevelStep;
visualizedConeAperture = p.Results.visualizedConeAperture;

%% Set up cone coordinates and outline
sampledHexMosaicXaxis = obj.patternSupport(1,:,1) + obj.center(1);
sampledHexMosaicYaxis = obj.patternSupport(:,1,2) + obj.center(2);

% Choose the radius of the aperture obj.pigment.pdWidth or obj.pigment.width
if (strcmp(visualizedConeAperture, 'lightCollectingArea'))
    % Note that pigment.pdWidth defines the size of a square collective
    % aperture. Here we compute the equivalent circular aperture
    dx = sqrt((obj.pigment.pdWidth^2)/pi)*2;
elseif (strcmp(visualizedConeAperture, 'geometricArea'))
    dx = obj.pigment.width;
end

if (showCorrespondingRectangularMosaicInstead)
    titleString = sprintf('<RECT grid> cones: %d x %d (%d total)', ...
        size(obj.patternOriginatingRectGrid,2), size(obj.patternOriginatingRectGrid,1), numel(obj.patternOriginatingRectGrid));
else
    titleString = sprintf('<RECT grid> cones: %d x %d (%d total), <HEX grid> cones: %d (active), %d (total), resampling factor: %d, visualized aperture: %s', ...
        size(obj.patternOriginatingRectGrid,2), size(obj.patternOriginatingRectGrid,1), numel(obj.patternOriginatingRectGrid), ...
        numel(find(obj.pattern > 1)), numel(obj.pattern), ...
        obj.resamplingFactor, visualizedConeAperture);
end

pixelOutline.x = [-1 -1 1 1 -1]*dx/2;
pixelOutline.y = [-1 1 1 -1 -1]*dx/2;
originalPixelOutline.x = [-1 -1 1 1 -1]*dx/2.0;
originalPixelOutline.y = [-1 1 1 -1 -1]*dx/2.0;

iTheta = (0:20:360)/180*pi;
apertureOutline.x = dx/2.0 * cos(iTheta);
apertureOutline.y = dx/2.0 * sin(iTheta);

rectCoords = obj.coneLocsOriginatingRectGrid;
hexCoords = obj.coneLocsHexGrid;

%% Set up figure
axesHandle = p.Results.axesHandle;
if (isempty(axesHandle))
    if (generateNewFigure)
        hFig = figure(round(rand()*100000));
        if (isempty(panelPosition))
            figPosition = [rand()*2000 rand()*1000 1300 1300];
        else
            figPosition = [(panelPosition(1)-1)*980 (panelPosition(2)-1)*700 1300 1300];
        end
    else
        % We want to use the coneMosaic window 
        if (isempty(panelPosition))
            hFig = figure(1);
            figPosition = [rand()*2000 rand()*1000 1300 1300];
        else
            hFig = figure(panelPosition(1)*10+panelPosition(2));
            figPosition = [(panelPosition(1)-1)*980 (panelPosition(2)-1)*700  1300 1300];
        end
    end
    cla;
    set(hFig, 'Position', figPosition, 'Color', [1 1 1]); % , 'MenuBar', 'none', 'NumberTitle', 'off');
    set(hFig, 'Name', titleString);
    subplot('Position', [0.02 0.01 0.97 0.97]);
    axesHandle = gca;
end

hold(axesHandle, 'on');

%% Do the display
switch showConeDensityContour
    case 'measured'
        [densityMap, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from mosaic');
    case 'theoretical'
        [densityMap, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from model');
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
        renderPatchArray(axesHandle, pixelOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    end
    
    % L-cones
    idx = find(obj.pattern == 2);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = 'none'; % [1 0 0]; 
    faceColor = [1.0 0. 0.];
    renderPatchArray(axesHandle, apertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    
    % M-cones
    idx = find(obj.pattern == 3);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = 'none';% = [0 0.7 0]; 
    faceColor = [0. 0.7 0.];
    renderPatchArray(axesHandle, apertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    
    % S-cones
    idx = find(obj.pattern == 4);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = 'none';% = [0 0 1]; 
    faceColor = [0. 0. 1.0];
    renderPatchArray(axesHandle, apertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    
    if (showPerfectHexMesh)
        % Superimpose hex mesh showing the locations of the perfect hex grid
        meshFaceColor = [0.8 0.8 0.8]; meshEdgeColor = [0.5 0.5 0.5]; meshFaceAlpha = 0.0; meshEdgeAlpha = 0.5; lineStyle = '-';
        renderHexMesh(axesHandle, hexCoords(:,1), hexCoords(:,2), meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle);
    end
else
    % Show the corresponding rectangular mosaic
    
    % The original rect sensors
    idx = find(obj.patternOriginatingRectGrid==2);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [1.0 0.7 0.7]; lineStyle = '-';
    renderPatchArray(axesHandle, originalPixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
    
    idx = find(obj.patternOriginatingRectGrid==3);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [0.7 1.0 0.7];
    renderPatchArray(axesHandle, originalPixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
    
    idx = find(obj.patternOriginatingRectGrid==4);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [0.7 0.7 1.0];
    renderPatchArray(axesHandle, originalPixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
end

if (~strcmp(showConeDensityContour, 'none'))
    contourLevels = coneDensityContourLevelStep: coneDensityContourLevelStep: 250000;
    plotContoursOverHalfField = true;
    if (plotContoursOverHalfField)
        idx = find(~((densityMapSupportX >= 0) & (densityMapSupportY >= 0)));
        densityMap(idx) = NaN;
    end
    [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, densityMap, contourLevels, 'LineColor', 'k', 'LineWidth', 1.5, 'ShowText', 'on', 'LabelSpacing', 2000);
    clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, 'Color', [0 0 0]);
    set(gca, 'CLim', [10000 250000]);
end

%% Arrange axis and fonts

hold(axesHandle, 'off')
axis(axesHandle, 'equal'); axis(axesHandle, 'xy')
xTicks = [sampledHexMosaicXaxis(1) obj.center(1) sampledHexMosaicXaxis(end)];
yTicks = [sampledHexMosaicYaxis(1) obj.center(2) sampledHexMosaicYaxis(end)];
xTickLabels = sprintf('%2.0f um\n', xTicks*1e6);
yTickLabels = sprintf('%2.0f um\n', yTicks*1e6);
set(axesHandle, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', {}, 'YTickLabel', {});
set(axesHandle, 'FontSize', 16, 'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 1.0);
box(axesHandle, 'on'); grid(axesHandle, 'off');
title(axesHandle, sprintf('%2.0f microns', obj.width*1e6), 'FontSize', 16);
set(axesHandle, 'XLim', [sampledHexMosaicXaxis(1)-dx sampledHexMosaicXaxis(end)+dx]);
set(axesHandle, 'YLim', [sampledHexMosaicYaxis(1)-dx sampledHexMosaicYaxis(end)+dx]);

drawnow;

end


%% Maybe put in utility directory
function renderPatchArray(axesHandle, pixelOutline, xCoords, yCoords, edgeColor, faceColor, lineStyle)

verticesNum = numel(pixelOutline.x);
x = zeros(verticesNum, numel(xCoords));
y = zeros(verticesNum, numel(xCoords));

for vertexIndex = 1:verticesNum
    x(vertexIndex, :) = pixelOutline.x(vertexIndex) + xCoords;
    y(vertexIndex, :) = pixelOutline.y(vertexIndex) + yCoords;
end
patch(x, y, [0 0 0], 'EdgeColor', edgeColor, 'FaceColor', faceColor, 'LineWidth', 0.2, 'LineStyle', lineStyle, 'Parent', axesHandle);
end

%% Separate function??
function renderHexMesh(axesHandle, xHex, yHex, meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle)
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
patch(x, y, [0 0 0], 'EdgeColor', meshEdgeColor, 'EdgeAlpha', meshEdgeAlpha, 'FaceAlpha', meshFaceAlpha, 'FaceColor', meshFaceColor, 'LineWidth', 1.5, 'LineStyle', lineStyle, 'Parent', axesHandle);
end

function Cout = getContourStruct(C)
    K = 0; n0 = 1;
    while n0<=size(C,2)
       K = K + 1;
       n0 = n0 + C(2,n0) + 1;
    end

    % initialize output struct
    el = cell(K,1);
    Cout = struct('level',el,'length',el,'x',el,'y',el);

    % fill the output struct
    n0 = 1;
    for k = 1:K
       Cout(k).level = C(1,n0);
       idx = (n0+1):(n0+C(2,n0));
       Cout(k).length = C(2,n0);
       Cout(k).x = C(1,idx);
       Cout(k).y = C(2,idx);
       n0 = idx(end) + 1; % next starting index
    end
end

