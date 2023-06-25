function hFig = visualizeGrid(obj, varargin)
% Visualize different aspects of the hexagonal cone mosaic
%
% Syntax:
%   hFig = visualizeGrid(obj, [varargin])
%
% Description:
%    Use this function to visualize different aspects of the hex grid.
%
% Inputs:
%    obj                       - The cone mosaic hex object
%
% Outputs:
%    hFig                      - The figure handle
%
% Optional key/value pairs:
%    axesHandle                - Handle, Axes handle to draw on. Default
%                                is Empty.
%    generateNewFigure         - Boolean, Whether or not to use a new
%                                figure? Default is False.
%    panelPosition             - 2 x 1. Panel position. Default is [1 1]
%    showCorrespondingRectangularMosaicInstead
%                              - Boolean, Whether or not to show the
%                                rectangular mosaic in place of the hex.
%                                Default is False.
%    visualizedConeAperture    - Char. Which element to visualize.
%                                Choose b/n: {'lightCollectingArea',
%                                             'geometricArea', 'both'}
%                                Default is 'lightCollectingArea'
%    apertureShape             - Char. What shape to use for cones
%                                Choose b/n: {'hexagons', 'disks'}
%                                Default is 'hexagons'
%    overlayNullSensors        - Boolean. Whether or not to overlay the
%                                null sensors. Default is False.
%    overlayEMpathMicrons      - 2 x 1. a single EM path (specified in
%                                microns) that can be overlayed on top of
%                                the mosaic. Default is Empty.
%    overlayHexMesh            - Boolean. Whether or not to overlay the
%                                hex mesh. Default is False.
%    overlayConeDensityContour - Char. Options are {'theoretical',
%                                'measured', 'theoretical_and_measured',
%                                'none'} Default is 'none'.
%    coneDensityContourLevels  - Array. Contour levels of cone density.
%                                Default is [100:20:250] * 1000.
%    overlayContourLabels      - Whether to label the contrours. Default: false
%    backgroundColor           - Background color. Default: [0.75 0.75 0.75]
%    foregroundColor           - Foreground (axes) color. Default; [0 0 0]
%

% History:
%    xx/xx/15  NPC  ISETBIO TEAM, 2015
%    02/21/18  jnm  Formatting
%    04/23/18  npc  Added coverage, innerSegmentCoverage properties
%    5/26/20   NPC  Fixed iRow issue, which was causing the mosaic plotting Y-coord flip 
%                   (rows grow top -> bottom, whereas Y-coords grow bottom -> top)

%% parse input
p = inputParser;
p.addParameter('generateNewFigure', false, @islogical);
p.addParameter('panelPosition', [1 1]);
p.addParameter('axesHandle', []);
p.addParameter('labelConeTypes', true, @islogical);
p.addParameter('showCorrespondingRectangularMosaicInstead', ...
    false, @islogical);
p.addParameter('visualizedConeAperture', 'geometricArea', ...
    @(x)ismember(x, {'lightCollectingArea', 'geometricArea', 'both'}));
p.addParameter('apertureShape', 'disks', @(x)ismember(x, ...
    {'hexagons', 'disks'}));
p.addParameter('overlayNullSensors', false, @islogical);
p.addParameter('overlayEMpathMicrons', [], @(x)(isnumeric(x) && ...
    ((isempty(x)) || (ndims(x) == 2))));
p.addParameter('overlayHexMesh', false, @islogical);
p.addParameter('overlayConeDensityContour', 'none', @(x)ismember(x, ...
    {'none', 'theoretical', 'measured', 'theoretical_and_measured'}));
p.addParameter('coneDensityContourLevels', (100:20:250) * 1000, ...
    @isnumeric);
p.addParameter('overlayContourLabels', false, @islogical);
p.addParameter('backgroundColor', [0.75 0.75 0.75]);
p.addParameter('foregroundColor', [0 0 0]);
p.addParameter('visualizedFOV', [], @(x)(isnumeric(x)||(isstruct(x))));
p.addParameter('ticksInVisualDegs', false, @islogical);
p.addParameter('ticksInMicrons', false, @islogical);
p.addParameter('tickInc', [], @isnumeric);
p.addParameter('noXaxisLabel', false, @islogical);
p.addParameter('noYaxisLabel', false, @islogical);
p.addParameter('scaleBarLengthMicrons', [], @(x)(isnumeric(x)));
p.addParameter('plotTitle', '', @ischar);
p.parse(varargin{:});

showCorrespondingRectangularMosaicInstead = ...
    p.Results.showCorrespondingRectangularMosaicInstead;
showNullSensors = p.Results.overlayNullSensors;
overlaidEMpathMicrons = p.Results.overlayEMpathMicrons;
overlayHexMesh = p.Results.overlayHexMesh;
overlayConeDensityContour = p.Results.overlayConeDensityContour;
generateNewFigure = p.Results.generateNewFigure;
panelPosition = p.Results.panelPosition;
coneDensityContourLevels = p.Results.coneDensityContourLevels;
visualizedConeAperture = p.Results.visualizedConeAperture;
visualizedFOV = p.Results.visualizedFOV;
apertureShape = p.Results.apertureShape;
labelConeTypes = p.Results.labelConeTypes;
backgroundColor = p.Results.backgroundColor;
foregroundColor = p.Results.foregroundColor;
scaleBarLengthMicrons = p.Results.scaleBarLengthMicrons;
ticksInMicrons = p.Results.ticksInMicrons;
ticksInVisualDegs = p.Results.ticksInVisualDegs;
tickInc = p.Results.tickInc;
noXaxisLabel = p.Results.noXaxisLabel;
noYaxisLabel = p.Results.noYaxisLabel;
if (p.Results.overlayContourLabels)
    overlayContourLabels = 'on';
else
    overlayContourLabels = 'off';
end
plotTitle = p.Results.plotTitle;

%% Set up cone coordinates and outline
sampledHexMosaicXaxis = obj.patternSupport(1, :, 1) + obj.center(1);
sampledHexMosaicYaxis = obj.patternSupport(:, 1, 2) + obj.center(2);

% Choose the radius of aperture obj.pigment.pdWidth or obj.pigment.width
if (strcmp(visualizedConeAperture, 'lightCollectingArea'))
    % Note that pigment.pdWidth defines the size of a square collective
    % aperture. Here we compute the equivalent circular aperture
    dxInner = diameterForCircularApertureFromWidthForSquareAperture(...
        obj.pigment.pdWidth);
    dxOuter = [];
elseif (strcmp(visualizedConeAperture, 'geometricArea'))
    dxOuter = diameterForCircularApertureFromWidthForSquareAperture(...
        obj.pigment.width);
    dxInner = [];
elseif (strcmp(visualizedConeAperture, 'both'))
    dxInner = diameterForCircularApertureFromWidthForSquareAperture(...
        obj.pigment.pdWidth);
    dxOuter = diameterForCircularApertureFromWidthForSquareAperture(...
        obj.pigment.width);
end


if (showCorrespondingRectangularMosaicInstead)
    titleString = sprintf('<RECT grid> cones: %d x %d (%d total)', ...
        size(obj.patternOriginatingRectGrid, 2), ...
        size(obj.patternOriginatingRectGrid, 1), ...
        numel(obj.patternOriginatingRectGrid));
else
    titleString = sprintf(['cones: %d (LMS), %d (LMSK), resampleF: %d,' ...
        ' aperture: %s'], numel(find(obj.pattern > 1)), ...
        numel(obj.pattern), obj.resamplingFactor, visualizedConeAperture);
end


% The outline of pixels in the original rect grid
pixelOutline.x = [0 0 1 1 0] * obj.patternSampleSize(1);
pixelOutline.y = [0 1 1 0 0] * obj.patternSampleSize(1);

if strcmp(apertureShape, 'hexagons')
    iTheta = ((0:60:360)+obj.rotationDegs)/180*pi;
else
    iTheta = (0:10:360) / 180 * pi;
end
if (~isempty(dxOuter))
    outerApertureOutline.x = dxOuter / 2.0 * cos(iTheta);
    outerApertureOutline.y = dxOuter / 2.0 * sin(iTheta);
else
    outerApertureOutline = [];
end
if (~isempty(dxInner))
    innerApertureOutline.x = dxInner / 2.0 * cos(iTheta);
    innerApertureOutline.y = dxInner / 2.0 * sin(iTheta);
else
    innerApertureOutline = [];
end

rectCoords = obj.coneLocsOriginatingRectGrid;
hexCoords = obj.coneLocsHexGrid;


if (isstruct(visualizedFOV))
    clippingRect = visualizedFOV;
    
    % Check that the clippingRect is valid
    coneMosaic.validateClippingRect(clippingRect);
    rr = abs(bsxfun(@minus, hexCoords, [clippingRect.xo clippingRect.yo]*obj.micronsPerDegree*1e-6));
    % 5 micron margin
    marginMeters = 5*1e-6;
    % Find cone indices within the ROI
    coneIndicesWithinROI = find((rr(:,1) <= clippingRect.width/2*obj.micronsPerDegree*1e-6+marginMeters) & (rr(:,2) <= clippingRect.height/2*obj.micronsPerDegree*1e-6+marginMeters));
    hexCoords = hexCoords(coneIndicesWithinROI,:); 
    xxx = obj.patternSupport(:,:,1);
    yyy = obj.patternSupport(:,:,2);
    xxx = abs(xxx(:) - clippingRect.xo*obj.micronsPerDegree*1e-6);
    yyy = abs(yyy(:) - clippingRect.yo*obj.micronsPerDegree*1e-6);
    pixelIndicesWithinROI = (xxx <= clippingRect.width/2*obj.micronsPerDegree*1e-6+marginMeters) &  (yyy <= clippingRect.height/2*obj.micronsPerDegree*1e-6+marginMeters);
else
    pixelIndicesWithinROI = [];
end


%% Set up figure
axesHandle = p.Results.axesHandle;
if (isempty(axesHandle))
    if (generateNewFigure)
        hFig = figure(round(rand() * 100000));
        if (isempty(panelPosition))
            figPosition = [rand() * 2000, rand() * 1000, 750, 750];
        else
            figPosition = [(panelPosition(1) - 1) * 980, ...
                (panelPosition(2) - 1) * 700, 750, 750];
        end
    else
        % We want to use the coneMosaic window
        if (isempty(panelPosition))
            hFig = figure(1);
            figPosition = [rand() * 2000, rand() * 1000, 750, 750];
        else
            hFig = figure(panelPosition(1) * 10 + panelPosition(2));
            figPosition = [(panelPosition(1) - 1) * 980, ...
                (panelPosition(2) - 1) * 700, 750, 750];
        end
    end
    cla;
    
    set(hFig, 'Position', figPosition, 'Color', backgroundColor);
    set(hFig, 'Name', titleString);
    subplot('Position', [0.1 0.04 0.89 0.92]);
    axesHandle = gca;
else
    hFig = get(gca,'Parent');
end


hold(axesHandle, 'on');

%% Do the display
if (overlayHexMesh)
    % Superimpose hex mesh showing the locations of the perfect hex grid
    meshFaceColor = [0.8 0.8 0.8]; meshEdgeColor = [0.5 0.5 0.5];
    meshFaceAlpha = 0.0; meshEdgeAlpha = 0.5; lineStyle = '-';
    coneMosaicHex.renderHexMesh(axesHandle, hexCoords(:,1), hexCoords(:,2), ...
        meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle);
    
end

if (~showCorrespondingRectangularMosaicInstead)
    lineStyle = '-';
    lineWidth = 0.2;
    if (showNullSensors)
        idx = find(obj.pattern == 1);
        [iRows, iCols] = ind2sub(size(obj.pattern), idx);
        edgeColor = [0.4 0.4 0.4];
        faceColor = 'none';
        coneMosaicHex.renderPatchArray(axesHandle, pixelOutline, ...
            sampledHexMosaicXaxis(iCols),  sampledHexMosaicYaxis(end-iRows+1), ...
            edgeColor, faceColor, lineStyle, lineWidth);
    end
    
    % L-cones
    if (isempty(pixelIndicesWithinROI))
        idx = find(obj.pattern == 2);
    else
        idx = find((obj.pattern(:) == 2) & (pixelIndicesWithinROI));
    end

    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords =  sampledHexMosaicYaxis(end-iRows+1);
    
    edgeColor = 'none'; % [1 0 0];
    if (labelConeTypes)
        faceColorInner = [1 0 0];
        faceColorOuter = [1 0 0];
    else
        edgeColor = [0 0 0];
        faceColorInner = 0.3*[1 1 1];
        faceColorOuter = 0.3*[1 1 1];
    end
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [innerApertureOutlineVarying, outerApertureOutlineVarying] = coneMosaicHex.computeApertureSizes(...
            dxInner, dxOuter, ...
            innerApertureOutline, outerApertureOutline, ...
            coneXcoords, coneYcoords ...
        );
        if (~isempty(outerApertureOutlineVarying))
            coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutlineVarying, ...
            coneXcoords, coneYcoords, ...
            edgeColor, faceColorOuter, lineStyle, lineWidth);
        end
        if (~isempty(innerApertureOutlineVarying))
            coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutlineVarying, ...
                coneXcoords, coneYcoords, ...
                edgeColor, faceColorInner, lineStyle, lineWidth);
        end
    else
        if (~isempty(outerApertureOutline))
            coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutline, ...
            coneXcoords, coneYcoords, ...
            edgeColor, faceColorOuter, lineStyle, lineWidth);
        end
        if (~isempty(innerApertureOutline))
            coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutline, ...
                coneXcoords, coneYcoords, ...
                edgeColor, faceColorInner, lineStyle, lineWidth);
        end
    end

    
    
    % M-cones
    if (isempty(pixelIndicesWithinROI))
        idx = find(obj.pattern == 3);
    else
        idx = find((obj.pattern(:) == 3) & (pixelIndicesWithinROI));
    end
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords =  sampledHexMosaicYaxis(end-iRows+1);
    
    edgeColor = 'none';  % = [0 0.7 0];
    if (labelConeTypes)
        if (mean(backgroundColor) < 0.5)
            faceColorInner = [0 1 0];
            faceColorOuter = [0.2 1 0.2];
        else
            faceColorInner = [0 1 0];
            faceColorOuter = [0 1 0];
        end
    else
        edgeColor = [0 0 0];
        faceColorInner = 0.3*[1 1 1];
        faceColorOuter = 0.3*[1 1 1];
    end
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [innerApertureOutlineVarying, outerApertureOutlineVarying] = coneMosaicHex.computeApertureSizes(...
            dxInner, dxOuter, ...
            innerApertureOutline, outerApertureOutline, ...
            coneXcoords, coneYcoords ...
        );
        if (~isempty(outerApertureOutlineVarying))
            coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutlineVarying, ...
            coneXcoords, coneYcoords, ...
            edgeColor, faceColorOuter, lineStyle, lineWidth);
        end
        if (~isempty(innerApertureOutlineVarying))
            coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutlineVarying, ...
                coneXcoords, coneYcoords, ...
                edgeColor, faceColorInner, lineStyle, lineWidth);
        end
    else
        if (~isempty(outerApertureOutline))
            coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutline, ...
                coneXcoords, coneYcoords, ...
                edgeColor, faceColorOuter, lineStyle, lineWidth);
        end
        if (~isempty(innerApertureOutline))
            coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutline, ...
                coneXcoords, coneYcoords, ...
                edgeColor, faceColorInner, lineStyle, lineWidth);
        end
    end
    
    % S-cones
    if (isempty(pixelIndicesWithinROI))
        idx = find(obj.pattern == 4);
    else
        idx = find((obj.pattern(:) == 4) & (pixelIndicesWithinROI));
    end
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords =  sampledHexMosaicYaxis(end-iRows+1);
    
    edgeColor = 'none';  % = [0 0 1];
    if (labelConeTypes)
        if (mean(backgroundColor) < 0.5)
            faceColorInner = [0 .4 1];
            faceColorOuter = [0.1 0.4 1];
        else
            faceColorInner = [0 0 1];
            faceColorOuter = [0 0 1];
        end
    else
        edgeColor = [0 0 0];
        faceColorInner = 0.3*[1 1 1];
        faceColorOuter = 0.3*[1 1 1];
    end
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [innerApertureOutlineVarying, outerApertureOutlineVarying] = coneMosaicHex.computeApertureSizes(...
            dxInner, dxOuter, ...
            innerApertureOutline, outerApertureOutline, ...
            coneXcoords, coneYcoords ...
        );
        if (~isempty(outerApertureOutlineVarying))
            coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutlineVarying, ...
            coneXcoords, coneYcoords, ...
            edgeColor, faceColorOuter, lineStyle, lineWidth);
        end
        if (~isempty(innerApertureOutlineVarying))
            coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutlineVarying, ...
                coneXcoords, coneYcoords, ...
                edgeColor, faceColorInner, lineStyle, lineWidth);
        end
    else
        if (~isempty(outerApertureOutline))
            coneMosaicHex.renderPatchArray(axesHandle, outerApertureOutline, ...
                coneXcoords, coneYcoords, ...
                edgeColor, faceColorOuter, lineStyle, lineWidth);
        end
        if (~isempty(innerApertureOutline))
            coneMosaicHex.renderPatchArray(axesHandle, innerApertureOutline, ...
                coneXcoords, coneYcoords, ...
                edgeColor, faceColorInner, lineStyle, lineWidth);
        end
    end
    
else
    lineWidth = 0.5;
    % Show the corresponding rectangular mosaic
    % The original rect sensors
    idx = find(obj.patternOriginatingRectGrid==2);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [1.0 0.7 0.7]; lineStyle = '-';
    coneMosaicHex.renderPatchArray(axesHandle, pixelOutline, ...
        rectCoords(idx,1), rectCoords(idx,2), edgeColor, ...
        faceColor, lineStyle, lineWidth);
    
    idx = find(obj.patternOriginatingRectGrid==3);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [0.7 1.0 0.7];
    coneMosaicHex.renderPatchArray(axesHandle, pixelOutline, ...
        rectCoords(idx,1), rectCoords(idx,2), edgeColor, ...
        faceColor, lineStyle, lineWidth);
    
    idx = find(obj.patternOriginatingRectGrid==4);
    %[iRows,iCols] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
    edgeColor = [0.3 0.3 0.3]; faceColor = [0.7 0.7 1.0];
    coneMosaicHex.renderPatchArray(axesHandle, pixelOutline, ...
        rectCoords(idx,1), rectCoords(idx,2), edgeColor, ...
        faceColor, lineStyle, lineWidth);
end

contourLevels = coneDensityContourLevels;
contourLabelSpacing = 4000;

plotContoursOverHalfField = false;

switch overlayConeDensityContour
    case 'measured'
        [densityMapMeasured, densityMapSupportX, densityMapSupportY] = ...
            obj.computeDensityMap('from mosaic');
        if (plotContoursOverHalfField)
            idx = find(~((densityMapSupportX >= 0) & ...
                (densityMapSupportY >= 0)));
            densityMapMeasured(idx) = NaN;
        end
        [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, ...
            densityMapMeasured, contourLevels, 'LineColor', 'k', 'LineWidth', 2.0, ...
            'ShowText', overlayContourLabels, 'LabelSpacing', contourLabelSpacing);
        clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, ...
            'Color', [1 0 0], 'BackgroundColor', [1 1 1]);
        set(gca, 'CLim', [10000 250000]);
        
    case 'theoretical'
        [densityMapTheoretical, densityMapSupportX, ...
            densityMapSupportY] = obj.computeDensityMap('from model');
        if (plotContoursOverHalfField)
            idx = find(~((densityMapSupportX >= 0) & ...
                (densityMapSupportY >= 0)));
            densityMapTheoretical(idx) = NaN;
        end
        
        if (p.Results.overlayContourLabels)
            [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, ...
                densityMapTheoretical, contourLevels, 'LineColor', [0.0 0.4 1.0], ...
                'LineWidth', 3.0, 'ShowText', overlayContourLabels, ...
                'LabelSpacing', contourLabelSpacing);
            clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, ...
                'Color', [0 0 1], 'BackgroundColor', [1 1 1]);
        else
            contour(axesHandle, densityMapSupportX, densityMapSupportY, ...
                densityMapTheoretical, contourLevels, 'LineColor', [0.0 1.0 0.3], ...
                'LineWidth', 3.0);
            %clabel(cH,hH,'FontWeight','bold', 'FontSize', 1, 'Color', 'none', 'BackgroundColor', 'none');
        end
        set(gca, 'CLim', [10000 250000]);
        
    case 'theoretical_and_measured'
        [densityMapMeasured, densityMapSupportX, densityMapSupportY] = ...
            obj.computeDensityMap('from mosaic');
        [densityMapTheoretical, densityMapSupportX, ...
            densityMapSupportY] = obj.computeDensityMap('from model');
        if (plotContoursOverHalfField)
            idx = find(~((densityMapSupportX >= 0) & ...
                (densityMapSupportY >= 0)));
            densityMapMeasured(idx) = NaN;
        end
        
        if (p.Results.overlayContourLabels)
            [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, ...
                densityMapMeasured, contourLevels, 'LineColor', [1 0.0 0.0], ...
                'LineWidth', 3.0, 'ShowText', overlayContourLabels, ...
                'LabelSpacing', contourLabelSpacing);
            clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, 'Color', [1 0 0], ...
                'BackgroundColor', [1 1 1]);
        else
            contour(axesHandle, densityMapSupportX, densityMapSupportY, densityMapMeasured, ...
                contourLevels, 'LineColor', [1 0.0 0.0], 'LineWidth', 3.0);
        end
        
        if (plotContoursOverHalfField)
            idx = find(~((densityMapSupportX >= 0) & ...
                (densityMapSupportY >= 0)));
            densityMapTheoretical(idx) = NaN;
        end
        
        [cH, hH] = contour(axesHandle, densityMapSupportX, densityMapSupportY, ...
            densityMapTheoretical, contourLevels, 'LineColor', [0.0 1.0 0.3], ...
            'LineWidth', 3.0, 'ShowText', 'on', 'LabelSpacing', contourLabelSpacing);
        clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, 'Color', [0 0 1], 'BackgroundColor', [1 1 1]);
        set(gca, 'CLim', [10000 250000]);
end

set(gca, 'Color', backgroundColor, 'XColor', foregroundColor, 'YColor', foregroundColor);

if (~isempty(overlaidEMpathMicrons))
    color = 'k';
    if (~labelConeTypes), color = 'r'; end
    plot(overlaidEMpathMicrons(:, 1) * 1e-6, ...
        overlaidEMpathMicrons(:, 2) * 1e-6, 'k.-', 'Color', color, ...
        'LineWidth', 3);
end

if (~isempty(scaleBarLengthMicrons))
   scaleBarLengthMeters = scaleBarLengthMicrons*1e-6;
   scaleBarX = sampledHexMosaicXaxis(end) * 0.99 + [-scaleBarLengthMeters 0];
   scaleBarY = sampledHexMosaicYaxis(1) * 0.95 * [1 1];
   plot(scaleBarX, scaleBarY, 'w-', 'LineWidth', 10.0);
   plot(scaleBarX, scaleBarY, 'k-', 'LineWidth', 5.0);
end

%% Arrange axis and fonts
hold(axesHandle, 'off')
axis(axesHandle, 'xy');
axis(axesHandle, 'square');

if (ticksInVisualDegs)
    rangeDegs = max(obj.fov);
    if (~isempty(p.Results.tickInc))
        tickInc = p.Results.tickInc;
    else
        if (rangeDegs < 0.25)
            tickInc = 0.01;
        elseif (rangeDegs < 0.5)
            tickInc = 0.05;
        elseif (rangeDegs < 1.0)
            tickInc = 0.1;
        elseif (rangeDegs < 4.0)
            tickInc = 0.25;
        elseif (rangeDegs < 8.0)
            tickInc = 0.5;
        else
            tickInc = 1;
        end
    end
    
    xMax = round(rangeDegs/tickInc)*tickInc;
    ticksDegs = (-xMax:tickInc:xMax);
    ticksMeters = ticksDegs * obj.micronsPerDegree * 1e-6;
    if (tickInc < 0.1)
        tickLabels = sprintf('%1.2f\n', ticksDegs);
    else
        tickLabels = sprintf('%1.1f\n', ticksDegs);
    end
    set(axesHandle, 'XTick', ticksMeters, 'YTick', ticksMeters, ...
        'XTickLabel', tickLabels, 'YTickLabel', tickLabels);
    if (~noXaxisLabel)
        xlabel(axesHandle,'space (degs)');
    end
    if (~noYaxisLabel)
        ylabel(axesHandle,'space (degs)');
    end
elseif (ticksInMicrons)
    rangeMicrons = max(obj.fov) * obj.micronsPerDegree;
    if (~isempty(p.Results.tickInc))
        tickInc = p.Results.tickInc;
    else
        if (rangeMicrons < 20)
            tickInc = 5;
        elseif (rangeMicrons < 50)
            tickInc = 10;
        elseif (rangeMicrons < 100)
            tickInc = 20;
        elseif (rangeMicrons < 250)
            tickInc = 50;
        elseif (rangeMicrons < 400)
            tickInc = 100;
        elseif (rangeMicrons < 800)
            tickInc = 200;
        elseif (rangeMicrons < 1000)
            tickInc = 250;
        else
            tickInc = 500;
        end
    end
    xMax = round(rangeMicrons/tickInc)*tickInc;
    ticksMicrons = (-xMax:tickInc:xMax);
    tickLabels = sprintf('%02.0f\n', ticksMicrons);
    ticksMeters = ticksMicrons * 1e-6;
    set(axesHandle, 'XTick', ticksMeters, 'YTick', ticksMeters, ...
        'XTickLabel', tickLabels, 'YTickLabel', tickLabels);
    if (~noXaxisLabel)
        xlabel(axesHandle,'space (microns)');
    end
    if (~noYaxisLabel)
        ylabel(axesHandle,'space (microns)');
    end
end

axis(axesHandle,'equal')

if (isempty(visualizedFOV))
    set(axesHandle, 'XLim', [sampledHexMosaicXaxis(1)-1.5*1e-6 sampledHexMosaicXaxis(end)+1.5*1e-6]);
    set(axesHandle, 'YLim', [sampledHexMosaicYaxis(1)-1.5*1e-6 sampledHexMosaicYaxis(end)+1.5*1e-6]);
else
    if (isstruct(visualizedFOV))
        clippingRect = visualizedFOV;
        % Check that the clippingRect is valid
        coneMosaic.validateClippingRect(clippingRect);
        spatialExtentMetersX = (clippingRect.xo+clippingRect.width/2*[-1 1])*obj.micronsPerDegree*1e-6;
        spatialExtentMetersY = (clippingRect.yo+clippingRect.height/2*[-1 1])*obj.micronsPerDegree*1e-6;
        set(axesHandle, 'XLim', spatialExtentMetersX, 'YLim', spatialExtentMetersY);
    else
        spatialExtentMeters = 0.5 * visualizedFOV * [-1 1] * obj.micronsPerDegree * 1e-6;
        set(axesHandle, 'XLim', spatialExtentMeters, 'YLim', spatialExtentMeters);
    end
end

box(axesHandle, 'on'); grid(axesHandle, 'off');
set(axesHandle, 'FontSize', 18, 'LineWidth', 1.0);

if (~isempty(plotTitle))
    title(axesHandle, plotTitle);
end

end
