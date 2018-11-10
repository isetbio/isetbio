function renderActivationMap(obj, axesHandle, activation, varargin)
% Render (in the passed axesHandle) an activation map for the hex mosaic
%
% Syntax:
%   renderActivationMap(obj, axesHandle, activation, [varargin])
%
% Description:
%    Render (draw) an activation map for the hex mosaic on the passed
%    axesHandle figure.
%
% Inputs:
%    obj        - The cone mosaic hex object
%    axesHandle - The axes figure handle
%    activation - The activation to map
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/17  NPC  ISETBIO TEAM, 2017
%    02/20/18  jnm  Formatting

    p = inputParser;
    p.addRequired('axesHandle', @ishandle);
    p.addRequired('activation', @isnumeric);
    p.addParameter('signalName', ' ', @ischar);
    p.addParameter('visualizedConeAperture', 'lightCollectingArea', ...
        @(x)ismember(x, {'lightCollectingArea', 'geometricArea'}));
    p.addParameter('mapType', 'modulated hexagons', ...
        @(x)ismember(x, {'modulated hexagons', 'modulated disks'}));
    p.addParameter('colorMap', gray(1024), @isnumeric);
    p.addParameter('titleForColorBar', '', @ischar);
    p.addParameter('showColorBar', false, @islogical);
    p.addParameter('labelColorBarTicks', false, @islogical);
    p.addParameter('showXLabel', true, @islogical);
    p.addParameter('showYLabel', true, @islogical);
    p.addParameter('outlineConesAlongHorizontalMeridian', false, @islogical);
    p.addParameter('crossHairPosition', [], @isnumeric);
    p.addParameter('visualizedFOV', [], @isnumeric);
    p.addParameter('signalRange', [], @isnumeric);
    p.addParameter('xRange', [], @isnumeric);
    p.addParameter('yRange', [], @isnumeric);
    p.addParameter('activationTime', [], @isnumeric);
    p.addParameter('backgroundColor', [0 0 0], @isnumeric);
    p.parse(axesHandle, activation, varargin{:});

    axesHandle = p.Results.axesHandle;
    activation = p.Results.activation;
    mapType = p.Results.mapType;
    visualizedConeAperture = p.Results.visualizedConeAperture;
    colorMap = p.Results.colorMap;
    signalRange = p.Results.signalRange;
    xRange =  p.Results.xRange;
    yRange = p.Results.yRange;
    showColorBar = p.Results.showColorBar;
    labelColorBarTicks = p.Results.labelColorBarTicks;
    backgroundColor = p.Results.backgroundColor;
    crossHairPosition = p.Results.crossHairPosition;
    visualizedFOV = p.Results.visualizedFOV;
    titleForColorBar = p.Results.titleForColorBar;
    showXLabel = p.Results.showXLabel;
    showYLabel = p.Results.showYLabel;
    outlineConesAlongHorizontalMeridian = p.Results.outlineConesAlongHorizontalMeridian;
    
    if (any(size(activation) ~= size(obj.pattern)))    
       activation = obj.reshapeHex2DmapToHex3Dmap(activation);
    end

    sampledHexMosaicXaxis = squeeze(obj.patternSupport(1, :, 1)) + ...
        obj.center(1);
    sampledHexMosaicYaxis = squeeze(obj.patternSupport(:, 1, 2)) + ...
        obj.center(2);

    % Choose aperture radius from obj.pigment.pdWidth or obj.pigment.width
    if (strcmp(visualizedConeAperture, 'lightCollectingArea'))
        % Note that pigment.pdWidth defines the size of a square collective
        % aperture. Here we compute the equivalent circular aperture
        %dx = sqrt((obj.pigment.pdWidth ^ 2) / pi) * 2;
        dx = diameterForCircularApertureFromWidthForSquareAperture(...
            obj.pigment.pdWidth);
    
    elseif (strcmp(visualizedConeAperture, 'geometricArea'))
        %dx = obj.pigment.width;
        dx = diameterForCircularApertureFromWidthForSquareAperture(...
            obj.pigment.width);
    end

    if strcmp(mapType, 'modulated disks') 
        iTheta = (0:15:360) / 180 * pi;
        apertureOutline.x = dx / 2.0 * cos(iTheta);
        apertureOutline.y = dx / 2.0 * sin(iTheta);
    elseif strcmp(mapType, 'modulated hexagons')
        iTheta = (0:60:360) / 180 * pi;
        apertureOutline.x = dx / 2.0 * cos(iTheta);
        apertureOutline.y = dx / 2.0 * sin(iTheta);
    end

    activeConesActivations = activation(obj.pattern > 1);
    if (isempty(signalRange))
        activationRange = [min(activeConesActivations(:)), ...
            max(activeConesActivations(:))];
    else
        activationRange = signalRange;
    end

    if (isempty(xRange))
        xRange = ...
            [sampledHexMosaicXaxis(1) sampledHexMosaicXaxis(end)] + ...
            1.5*1e-6*[-1 1];
    else
        xRange = xRange * 1e-6;
    end
    if (isempty(yRange))
        yRange = [sampledHexMosaicYaxis(1) sampledHexMosaicYaxis(end)] + ...
            1.5*1e-6*[-1 1];
     else
        yRange = yRange * 1e-6;
    end

    
    cMapLevels = size(colorMap, 1);
    idx = find(obj.pattern > 1);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);  
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(iRows);
    activationsNlevels = round((activation(idx) - activationRange(1)) / ...
        (activationRange(2) - activationRange(1)) * cMapLevels);
    faceColorsNormalizedValues = activationsNlevels / cMapLevels;
    faceColorsNormalizedValues(faceColorsNormalizedValues < 0) = 0;
    faceColorsNormalizedValues(faceColorsNormalizedValues > 1) = 1;

    edgeColor = [0 0 0];
    lineWidth = 0.1;

    if (obj.shouldCorrectAbsorptionsWithEccentricity())
    % Compute ecc-varying apertures
        [apertureOutline, ~, maxAperture] = coneMosaicHex.computeApertureSizes(...
            dx, [], apertureOutline, [], coneXcoords, coneYcoords);
    end

    hold(axesHandle, 'on');
    x = [xRange(1) xRange(2) xRange(2) xRange(1)];
    y = [yRange(1) yRange(1) yRange(2) yRange(2)];
    patch(x, y, [0 0 0], 'FaceColor', 'none')
    renderModulatedColorPatchArray(axesHandle, apertureOutline, ...
        coneXcoords, coneYcoords, ...
        faceColorsNormalizedValues, edgeColor, lineWidth);
    if (outlineConesAlongHorizontalMeridian)
        coneXcoordsDegs = coneXcoords * 1e6 / obj.micronsPerDegree;
        coneYcoordsDegs = coneYcoords * 1e6 / obj.micronsPerDegree;
        % Find cones lying near the y=0 axis
        indicesOfConesAlongXaxis = find(abs(coneYcoordsDegs) < dx * 1e6 / obj.micronsPerDegree);
        coneXcoordsDegs = coneXcoordsDegs(indicesOfConesAlongXaxis);
        coneYcoordsDegs = coneYcoordsDegs(indicesOfConesAlongXaxis);
        identitiesOfConesAlongXaxis = obj.pattern(idx(indicesOfConesAlongXaxis));
        apertureOutline.x = apertureOutline.x(indicesOfConesAlongXaxis,:);
        apertureOutline.y = apertureOutline.y(indicesOfConesAlongXaxis,:);
        for kkk = 1:numel(identitiesOfConesAlongXaxis)
            switch (identitiesOfConesAlongXaxis(kkk))
                case 2
                    color = [1 0.2 0.1];
                case 3
                    color = [0 1 0];
                case 4
                    color = [0 0.6 1];    
            end
            plot(apertureOutline.x(kkk,:)+coneXcoords(indicesOfConesAlongXaxis(kkk)), ...
                apertureOutline.y(kkk,:)+coneYcoords(indicesOfConesAlongXaxis(kkk)), ...
                'r-', 'Color', color, 'LineWidth', 2.0);
        end
        %plot(coneXcoords(indicesOfConesAlongXaxis), coneYcoords(indicesOfConesAlongXaxis), 'ro');
    end
    
    xlim(axesHandle, xRange);
    ylim(axesHandle, yRange);
    
    if (~isempty(crossHairPosition))
        plot(axesHandle, [xRange(1) xRange(end)], -crossHairPosition(2)*[1 1], 'g-', 'LineWidth', 1.5);
        plot(axesHandle, crossHairPosition(1)*[1 1],[yRange(1) yRange(end)], 'g-', 'LineWidth', 1.5);
    end
    
    hold(axesHandle, 'off');
    set(axesHandle, 'CLim', [0 1], 'XTick', [], 'YTick', [], ...
        'Color', backgroundColor);

    colormap(axesHandle, colorMap);     
    if (showColorBar)
        hC = colorbar();
        if (labelColorBarTicks)
            ticks = [0 0.5 1.0];
            hC.Ticks = ticks; 
            tickLabels = (ticks * (activationRange(2) - activationRange(1)) + activationRange(1));
            
            if (max(abs(tickLabels) < 0.01))
                if (min(tickLabels) < 0)
                    hC.TickLabels = sprintf('%+2.2f\n',tickLabels); 
                else
                    hC.TickLabels = sprintf('%2.2f\n',tickLabels); 
                end
            elseif (max(abs(tickLabels) < 0.4))
                if (min(tickLabels) < 0)
                    hC.TickLabels = sprintf('%+2.1f\n',tickLabels); 
                else
                    hC.TickLabels = sprintf('%2.1f\n',tickLabels); 
                end
            else
                if (min(tickLabels) < 0)
                    hC.TickLabels = sprintf('%+2.0f\n',round(tickLabels)); 
                else
                    hC.TickLabels = sprintf('%2.0f\n',round(tickLabels)); 
                end
            end
        else
            hC.TickLabels = {}; 
        end
        if (~isempty(titleForColorBar))
            hC.Label.String = titleForColorBar;
        end
    end
    
    axis(axesHandle, 'image'); 
    axis(axesHandle, 'xy');
    box(axesHandle, 'on');
    
    if (isempty(visualizedFOV))
        visualizedFOV = max(obj.fov);
    end
    ticksDegs = round(100*0.5 * visualizedFOV * 1.0*[-1 -0.5 0 0.5 1])/100;
    ticksMeters = ticksDegs * obj.micronsPerDegree * 1e-6;
    spatialExtentMeters = 0.5 * visualizedFOV * [-1 1] * obj.micronsPerDegree * 1e-6 + maxAperture/2*[-1 1];
    
    if (showXLabel)
        xlabel(axesHandle, 'space (degs)', 'FontWeight', 'bold');
    end
    if (showYLabel)
        ylabel(axesHandle, 'space (degs)', 'FontWeight', 'bold');
    end
    tickLabels = {sprintf('%2.2f', ticksDegs(1)), sprintf('%2.2f', ticksDegs(2)), sprintf('%2.2f', ticksDegs(3))};
    
    set(axesHandle, 'XTick', ticksMeters, 'YTick', ticksMeters, ...
        'XLim', spatialExtentMeters,  'YLim', spatialExtentMeters, ...
        'XTick', ticksMeters, 'YTick', ticksMeters, ...
        'XTickLabels', sprintf('%2.2f\n', ticksDegs), 'YTickLabels', sprintf('%2.2f\n', ticksDegs));
    xticks(axesHandle, ticksMeters); yticks(axesHandle, ticksMeters); 
    set(axesHandle, 'FontSize', 18, 'LineWidth', 1.0);
    yticklabels(axesHandle, tickLabels);
    
end

function renderModulatedColorPatchArray(axesHandle, pixelOutline, xCoords, yCoords, ...
    faceColorsNormalizedValues, edgeColor, lineWidth)
% Render the patch array
%
% Syntax:
%   renderModulatedColorPatchArray(axesHandle, pixelOutline, xCoords, yCoords, ...
%       faceColorNormalizedValues, edgeColor, lineStyle)
%
% Description:
%    Render the patch array
%
% Inputs:
%    axesHandle   - The handle to the patch axes
%    pixelOutline - The outline of the pixels
%    xCoords      - The X-axis coordinates
%    yCoords      - The Y-axis coordinates
%    faceColorNormalizedValues
%                 - Normalized values for the face colors (fill)
%    edgeColor    - The color of the edge of the shape (line)
%    lineStyle    - The style of the lines (dash, solid, etc...)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    %verticesPerCone = numel(pixelOutline.x)
    verticesPerCone = size(pixelOutline.x,2);
    
    verticesList = zeros(verticesPerCone * numel(xCoords), 2);
    facesList = [];
    colors = [];
    
    for coneIndex = 1:numel(xCoords)
        idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
        
        if (size(pixelOutline.x,1) == 1)
            verticesList(idx, 1) = pixelOutline.x(:) + xCoords(coneIndex);
            verticesList(idx, 2) = pixelOutline.y(:) + yCoords(coneIndex);
        else
            verticesList(idx, 1) = squeeze(pixelOutline.x(coneIndex,:)) + xCoords(coneIndex);
            verticesList(idx, 2) = squeeze(pixelOutline.y(coneIndex,:)) + yCoords(coneIndex);
        end
        facesList = cat(1, facesList, idx);
        colors = cat(1, colors, ...
            repmat(faceColorsNormalizedValues(coneIndex), ...
            [verticesPerCone 1]));
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end
