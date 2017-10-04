function hFig = visualizeGrid(obj, varargin)
% Visualize different aspects of the hex grid
%
% Name-Value options
%   axesHandle - Empty - Axes handle to draw on
%   generateNewFigure  - False
%   panelPosition      - [1 1]
%   showCorrespondingRectangularMosaicInstead - False
%   overlayNullSensors         - False
%   overlayEMpath Microns  - a single EM path (specified in microns) that can be overlayed on top of the mosaic
%   overlayHexMesh       - False
%   overlayConeDensityContour   - 'none'
%   coneDensityContourLevelStep - 5000
%
% NPC, ISETBIO TEAM, 2015

%% parse input
p = inputParser;
p.addParameter('generateNewFigure', false, @islogical);
p.addParameter('panelPosition', [1 1]);
p.addParameter('axesHandle', []);
p.addParameter('labelConeTypes', true, @islogical);
p.addParameter('showCorrespondingRectangularMosaicInstead', false, @islogical);
p.addParameter('visualizedConeAperture', 'lightCollectingArea', @(x)ismember(x, {'lightCollectingArea', 'geometricArea', 'both'}));
p.addParameter('apertureShape', 'hexagons', @(x)ismember(x, {'hexagons', 'disks'}));
p.addParameter('overlayNullSensors', false, @islogical);
p.addParameter('overlayEMpathMicrons', [], @(x)(isnumeric(x) && ((isempty(x)) || (ndims(x)==2))));
p.addParameter('overlayHexMesh', false, @islogical);
p.addParameter('overlayConeDensityContour', 'none', @(x)ismember(x, {'none', 'theoretical', 'measured', 'theoretical_and_measured'}));
p.addParameter('coneDensityContourLevels', [100:20:250]*1000, @isnumeric);
p.parse(varargin{:});

showCorrespondingRectangularMosaicInstead = p.Results.showCorrespondingRectangularMosaicInstead;
showNullSensors = p.Results.overlayNullSensors;
overlaidEMpathMicrons = p.Results.overlayEMpathMicrons;
showHexMesh = p.Results.overlayHexMesh;
showConeDensityContour = p.Results.overlayConeDensityContour;
generateNewFigure = p.Results.generateNewFigure;
panelPosition = p.Results.panelPosition;
coneDensityContourLevels = p.Results.coneDensityContourLevels;
visualizedConeAperture = p.Results.visualizedConeAperture;
apertureShape = p.Results.apertureShape;
labelConeTypes = p.Results.labelConeTypes;

%% Set up cone coordinates and outline
sampledHexMosaicXaxis = obj.patternSupport(1,:,1) + obj.center(1);
sampledHexMosaicYaxis = obj.patternSupport(:,1,2) + obj.center(2);

% Choose the radius of the aperture obj.pigment.pdWidth or obj.pigment.width
if (strcmp(visualizedConeAperture, 'lightCollectingArea'))
    % Note that pigment.pdWidth defines the size of a square collective
    % aperture. Here we compute the equivalent circular aperture
    dxInner = diameterForCircularApertureFromWidthForSquareAperture(obj.pigment.pdWidth);
    dxOuter = [];
elseif (strcmp(visualizedConeAperture, 'geometricArea'))
    dxOuter = diameterForCircularApertureFromWidthForSquareAperture(obj.pigment.width);
    dxInner = [];
elseif (strcmp(visualizedConeAperture, 'both'))
    dxInner = diameterForCircularApertureFromWidthForSquareAperture(obj.pigment.pdWidth);
    dxOuter = diameterForCircularApertureFromWidthForSquareAperture(obj.pigment.width);
end

if (showCorrespondingRectangularMosaicInstead)
    titleString = sprintf('<RECT grid> cones: %d x %d (%d total)', ...
        size(obj.patternOriginatingRectGrid,2), size(obj.patternOriginatingRectGrid,1), numel(obj.patternOriginatingRectGrid));
else
    titleString = sprintf('cones: %d (LMS), %d (LMSK), resampleF: %d, aperture: %s', ...
        numel(find(obj.pattern > 1)), numel(obj.pattern), obj.resamplingFactor, visualizedConeAperture);
end

% The outline of pixels in the original rect grid
pixelOutline.x = [0 0 1 1 0]*obj.patternSampleSize(1);
pixelOutline.y = [0 1 1 0 0]*obj.patternSampleSize(1);
    

if strcmp(apertureShape, 'hexagons')
    iTheta = (0:60:360)/180*pi;
else
    iTheta = (0:10:360)/180*pi;
end
if (~isempty(dxOuter))  
    outerApertureOutline.x = dxOuter/2.0 * cos(iTheta);
    outerApertureOutline.y = dxOuter/2.0 * sin(iTheta);
else
    outerApertureOutline = [];
end
if (~isempty(dxInner))  
    innerApertureOutline.x = dxInner/2.0 * cos(iTheta);
    innerApertureOutline.y = dxInner/2.0 * sin(iTheta);
else
    innerApertureOutline = [];
end

rectCoords = obj.coneLocsOriginatingRectGrid;
hexCoords = obj.coneLocsHexGrid;

%% Set up figure
axesHandle = p.Results.axesHandle;
if (isempty(axesHandle))
    if (generateNewFigure)
        hFig = figure(round(rand()*100000));
        if (isempty(panelPosition))
            figPosition = [rand()*2000 rand()*1000 750 750];
        else
            figPosition = [(panelPosition(1)-1)*980 (panelPosition(2)-1)*700 750 750];
        end
    else
        % We want to use the coneMosaic window 
        if (isempty(panelPosition))
            hFig = figure(1);
            figPosition = [rand()*2000 rand()*1000 750 750];
        else
            hFig = figure(panelPosition(1)*10+panelPosition(2));
            figPosition = [(panelPosition(1)-1)*980 (panelPosition(2)-1)*700  750 750];
        end
    end
    cla;
    set(hFig, 'Position', figPosition, 'Color', [1 1 1]); % , 'MenuBar', 'none', 'NumberTitle', 'off');
    set(hFig, 'Name', titleString);
    subplot('Position', [0.1 0.04 0.89 0.92]);
    axesHandle = gca;
end

hold(axesHandle, 'on');

%% Do the display
switch showConeDensityContour
    case 'measured'
        [densityMapMeasured, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from mosaic');
    case 'theoretical'
        [densityMapTheoretical, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from model');
    case 'theoretical_and_measured'
        [densityMapMeasured, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from mosaic');
        [densityMapTheoretical, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from model');
    case 'none'
    otherwise
        error('coneMosaicHex.visualizeGrid: ''coneDensityContourOverlay'' must be set to one of the following: ''measured'', ''theoretical'', ''none''. ');
end

if (showHexMesh)
        % Superimpose hex mesh showing the locations of the perfect hex grid
        meshFaceColor = [0.8 0.8 0.8]; meshEdgeColor = [0.5 0.5 0.5]; meshFaceAlpha = 0.0; meshEdgeAlpha = 0.5; lineStyle = '-';
        coneMosaicHex.renderHexMesh(axesHandle, hexCoords(:,1), hexCoords(:,2), meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle);
end
    
if (~showCorrespondingRectangularMosaicInstead)
    lineStyle = '-';
    if (showNullSensors)
        idx = find(obj.pattern==1);
        [iRows,iCols] = ind2sub(size(obj.pattern), idx);
        edgeColor = [0.4 0.4 0.4]; faceColor = 'none';
        coneMosaicHex.renderPatchArray(axesHandle, pixelOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColor, lineStyle);
    end
    
    % L-cones
    idx = find(obj.pattern == 2);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = 'none'; % [1 0 0]; 
    if (labelConeTypes)
        faceColorInner = [1 0 0];
        faceColorOuter = [1 0.5 0.5];
    else
        edgeColor = [0 0 0]; 
        faceColorInner = 0.8*[1 1 1];
        faceColorOuter = 0.9*[1 1 1];
    end
    if (~isempty(outerApertureOutline))
        coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColorOuter, lineStyle);
    end
    if (~isempty(innerApertureOutline))
        coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColorInner, lineStyle);
    end
    
    % M-cones
    idx = find(obj.pattern == 3);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = 'none';% = [0 0.7 0]; 
    if (labelConeTypes)
        faceColorInner = [0 1 0];
        faceColorOuter = [0.5 1 0.5];
    else
        edgeColor = [0 0 0]; 
        faceColorInner = 0.8*[1 1 1];
        faceColorOuter = 0.9*[1 1 1];
    end
    
    if (~isempty(outerApertureOutline))
        coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColorOuter, lineStyle);
    end
    if (~isempty(innerApertureOutline))
        coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColorInner, lineStyle);
    end
    
    % S-cones
    idx = find(obj.pattern == 4);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);
    edgeColor = 'none';% = [0 0 1]; 
    if (labelConeTypes)
        faceColorInner = [0 0 1];
        faceColorOuter = [0.5 0.5 1];
    else
        edgeColor = [0 0 0]; 
        faceColorInner = 0.8*[1 1 1];
        faceColorOuter = 0.9*[1 1 1];
    end
    
    if (~isempty(outerApertureOutline))
        coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColorOuter, lineStyle);
    end
    if (~isempty(innerApertureOutline))
        coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutline, sampledHexMosaicXaxis(iCols), sampledHexMosaicYaxis(iRows), edgeColor, faceColorInner, lineStyle);
    end
    
else
    % Show the corresponding rectangular mosaic
    % The original rect sensors
    idx = find(obj.patternOriginatingRectGrid==2);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [1.0 0.7 0.7]; lineStyle = '-';
    coneMosaicHex.renderPatchArray(axesHandle, pixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
    
    idx = find(obj.patternOriginatingRectGrid==3);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [0.7 1.0 0.7];
    coneMosaicHex.renderPatchArray(axesHandle, pixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
    
    idx = find(obj.patternOriginatingRectGrid==4);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [0.7 0.7 1.0];
    coneMosaicHex.renderPatchArray(axesHandle, pixelOutline, rectCoords(idx,1), rectCoords(idx,2), edgeColor, faceColor, lineStyle);
end


contourLevels = coneDensityContourLevels;
plotContoursOverHalfField = false;
    
switch showConeDensityContour
    case 'measured'
        [densityMapMeasured, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from mosaic');
        if (plotContoursOverHalfField)
            idx = find(~((densityMapSupportX >= 0) & (densityMapSupportY >= 0)));
            densityMapMeasured(idx) = NaN;
        end
        [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, densityMapMeasured, contourLevels, 'LineColor', 'r', 'LineWidth', 2.0, 'ShowText', 'on', 'LabelSpacing', 2000);
        clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, 'Color', [1 0 0], 'BackgroundColor', [1 1 1]);
        set(gca, 'CLim', [10000 250000]);
    
    case 'theoretical'
        [densityMapTheoretical, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from model');
        if (plotContoursOverHalfField)
            idx = find(~((densityMapSupportX >= 0) & (densityMapSupportY >= 0)));
            densityMapTheoretical(idx) = NaN;
        end
        [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, densityMapTheoretical, contourLevels, 'LineColor', 'b', 'LineWidth', 2.0, 'ShowText', 'on', 'LabelSpacing', 2000);
        clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, 'Color', [0 0 1], 'BackgroundColor', [1 1 1]);
        set(gca, 'CLim', [10000 250000]);
        
    case 'theoretical_and_measured'
        [densityMapMeasured, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from mosaic');
        [densityMapTheoretical, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from model');
        if (plotContoursOverHalfField)
            idx = find(~((densityMapSupportX >= 0) & (densityMapSupportY >= 0)));
            densityMapMeasured(idx) = NaN;
        end
        [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, densityMapMeasured, contourLevels, 'LineColor', 'r', 'LineWidth', 2.0, 'ShowText', 'on', 'LabelSpacing', 2000);
        clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, 'Color', [1 0 0], 'BackgroundColor', [1 1 1]);
        if (plotContoursOverHalfField)
            idx = find(~((densityMapSupportX >= 0) & (densityMapSupportY >= 0)));
            densityMapTheoretical(idx) = NaN;
        end
        [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, densityMapTheoretical, contourLevels, 'LineColor', 'b', 'LineWidth', 2.0, 'ShowText', 'on', 'LabelSpacing', 3000);
        clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, 'Color', [0 0 1], 'BackgroundColor', [1 1 1]);
        set(gca, 'CLim', [10000 250000]);
end

if (~isempty(overlaidEMpathMicrons))
    color = 'k';
    if (~labelConeTypes)
        color = 'r';
    end
    plot(overlaidEMpathMicrons(:,1)*1e-6, overlaidEMpathMicrons(:,2)*1e-6, 'ks-', 'Color', color, 'LineWidth', 1.5);
end

%% Arrange axis and fonts

hold(axesHandle, 'off')
axis(axesHandle, 'equal'); axis(axesHandle, 'xy')

if (isempty(p.Results.axesHandle))
    xTicks = [sampledHexMosaicXaxis(1) obj.center(1) sampledHexMosaicXaxis(end)];
    yTicks = [sampledHexMosaicYaxis(1) obj.center(2) sampledHexMosaicYaxis(end)];
    xTickLabels = sprintf('%2.0f um\n', xTicks*1e6);
    yTickLabels = sprintf('%2.0f um\n', yTicks*1e6);
    set(axesHandle, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
    set(axesHandle, 'FontSize', 16, 'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 1.0);
    box(axesHandle, 'on'); grid(axesHandle, 'off');
    title(axesHandle, sprintf('%2.0f microns', obj.width*1e6), 'FontSize', 16);
    set(axesHandle, 'XLim', [sampledHexMosaicXaxis(1)-1.5*1e-6 sampledHexMosaicXaxis(end)+1.5*1e-6]);
    set(axesHandle, 'YLim', [sampledHexMosaicYaxis(1)-1.5*1e-6 sampledHexMosaicYaxis(end)+1.5*1e-6]);
    drawnow;
end

end