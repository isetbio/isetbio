function renderActivationMap(obj, axesHandle, activation, varargin)
% Render (in the passed axesHandle) an activation map for the hex mosaic
%
% NPC, ISETBIO TEAM, 2017

    p = inputParser;
    p.addRequired('axesHandle', @ishandle);
    p.addRequired('activation', @isnumeric);
    p.addParameter('signalName', ' ', @ischar);
    p.addParameter('visualizedConeAperture', 'lightCollectingArea', @(x)ismember(x, {'lightCollectingArea', 'geometricArea'}));
    p.addParameter('mapType', 'modulated hexagons', @(x)ismember(x, {'modulated hexagons', 'modulated disks'}));
    p.addParameter('colorMap', jet(1024), @isnumeric);
    p.addParameter('signalRange', [], @isnumeric);
    p.addParameter('xRange', [], @isnumeric);
    p.addParameter('yRange', [], @isnumeric);
    p.addParameter('activationTime', [], @isnumeric);
    p.parse(axesHandle, activation, varargin{:});
    
    axesHandle = p.Results.axesHandle;
    activation = p.Results.activation;
    mapType = p.Results.mapType;
    visualizedConeAperture = p.Results.visualizedConeAperture;
    signalName = p.Results.signalName;
    colorMap = p.Results.colorMap;
    signalRange = p.Results.signalRange;
    xRange =  p.Results.xRange;
    yRange = p.Results.yRange;
    activationTime = p.Results.activationTime;
    
    if (any(size(activation) ~= size(obj.pattern)))    
       activation = obj.reshapeHex2DmapToHex3Dmap(activation);
    end
    
    sampledHexMosaicXaxis = squeeze(obj.patternSupport(1,:,1)) + obj.center(1);
    sampledHexMosaicYaxis = squeeze(obj.patternSupport(:,1,2)) + obj.center(2);
    
    % Choose the radius of the aperture obj.pigment.pdWidth or obj.pigment.width
    if (strcmp(visualizedConeAperture, 'lightCollectingArea'))
        % Note that pigment.pdWidth defines the size of a square collective
        % aperture. Here we compute the equivalent circular aperture
        dx = sqrt((obj.pigment.pdWidth^2)/pi)*2;
    elseif (strcmp(visualizedConeAperture, 'geometricArea'))
        dx = obj.pigment.width;
    end
    
    if strcmp(mapType, 'modulated disks') 
        iTheta = (0:20:360)/180*pi;
        apertureOutline.x = dx/2.0 * cos(iTheta);
        apertureOutline.y = dx/2.0 * sin(iTheta);
    elseif strcmp(mapType, 'modulated hexagons')
        iTheta = (0:60:360)/180*pi;
        apertureOutline.x = dx/2.0 * cos(iTheta);
        apertureOutline.y = dx/2.0 * sin(iTheta);
    end
    
    activeConesActivations = activation(obj.pattern>1);
    if (isempty(signalRange))
        activationRange = [min(activeConesActivations(:)) max(activeConesActivations(:))];
    else
        activationRange = signalRange;
    end
    
    if (isempty(xRange))
        xRange = [sampledHexMosaicXaxis(1)-obj.pigment.width sampledHexMosaicXaxis(end)+obj.pigment.width];
    else
        xRange = xRange *1e-6;
    end
    if (isempty(yRange))
        yRange = [sampledHexMosaicYaxis(1)-obj.pigment.width sampledHexMosaicYaxis(end)+obj.pigment.width];
     else
        yRange = yRange *1e-6;
    end

    cMapLevels = size(colorMap,1);
    idx = find(obj.pattern > 1);
    [iRows,iCols] = ind2sub(size(obj.pattern), idx);  
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(iRows);
    
    activationsNlevels = round((activation(idx)-activationRange(1))/(activationRange(2)-activationRange(1))*cMapLevels);
    faceColorsNormalizedValues = activationsNlevels/cMapLevels;
    faceColorsNormalizedValues(faceColorsNormalizedValues<0) = 0;
    faceColorsNormalizedValues(faceColorsNormalizedValues>1) = 1;
            
    edgeColor = 'none';
    lineWidth = 0.1;

    hold(axesHandle, 'on');
    x = [xRange(1) xRange(2) xRange(2) xRange(1)];
    y = [yRange(1) yRange(1) yRange(2) yRange(2)];
    patch(x,y, [0 0 0], 'FaceColor', 'none')
    renderPatchArray(axesHandle, apertureOutline, coneXcoords, coneYcoords, faceColorsNormalizedValues,  edgeColor,  lineWidth);
    xlim(axesHandle, xRange);
    ylim(axesHandle, yRange);
    hold(axesHandle, 'off');
    set(axesHandle, 'CLim', [0 1], 'XTick', [], 'YTick', [], 'Color', 'none');
    colormap(axesHandle, colorMap);     
    axis(axesHandle, 'image'); 
    axis(axesHandle, 'xy');
    box(axesHandle, 'on');
end

function renderPatchArray(axesHandle, pixelOutline, xCoords, yCoords, faceColorsNormalizedValues,  edgeColor, lineWidth)
    verticesPerCone = numel(pixelOutline.x);
    verticesList = zeros(verticesPerCone * numel(xCoords), 2);
    facesList = [];
    colors = [];
    for coneIndex = 1:numel(xCoords)
        idx = (coneIndex-1)*verticesPerCone + (1:verticesPerCone);
        verticesList(idx,1) = pixelOutline.x(:) + xCoords(coneIndex);
        verticesList(idx,2) = pixelOutline.y(:) + yCoords(coneIndex);
        facesList = cat(1, facesList, idx);
        colors = cat(1, colors, repmat(faceColorsNormalizedValues(coneIndex), [verticesPerCone 1]));
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end

