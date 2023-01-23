function visualize(obj, varargin)
    
    % Parse optional input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('component', 'RF centers', @(x)ismember(x, {'RF centers'}));
    p.addParameter('activation', []);
    p.addParameter('activationRange', [],@(x)((isempty(x))||(numel(x)==2)));
    p.addParameter('activationColorMap', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.addParameter('backgroundColor', [], @(x)((ischar(x)&&(strcmp(x,'none')))||isempty(x)||((isvector(x))&&(numel(x) == 3))));
    p.addParameter('XLims', [], @isnumeric);
    p.addParameter('YLims', [], @isnumeric);
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    XLims = p.Results.XLims;
    YLims = p.Results.YLims;
    visualizedComponent = p.Results.component;
    activation = p.Results.activation;
    activationRange = p.Results.activationRange;
    activationColorMap = p.Results.activationColorMap;
    backgroundColor = p.Results.backgroundColor;



    switch (visualizedComponent)
        case 'RF centers'
            visualizeRFcenters(obj, hFig, ax, XLims, YLims, ...
                activation, activationRange, activationColorMap, backgroundColor);

        otherwise
            error('Uknown visualized component: ''%s''.', visualizedComponent);
    end
end


function visualizeRFcenters(obj,hFig, ax, XLims, YLims, ...
    activation, activationRange, activationColorMap, backgroundColor)

    coneCharacteristicRadiiDegs = ...
        obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
        obj.inputConeMosaic.coneApertureDiametersDegs;

    xSupport = [];
    ySupport = []; 
    spatialSupportSamples = 60;

    if (isempty(obj.rfCenterVisualizationOutlines))
        fprintf('\nPlease wait. Computing RF outlines for %d RF centers ... ', obj.rgcsNum);
        obj.rfCenterVisualizationOutlines = cell(1, obj.rgcsNum);
        for iRGC = 1:obj.rgcsNum
             % Retrieve the center cone indices & weights
             centerConnectivityVector = full(squeeze(obj.centerConePoolingMatrix(:, iRGC)));
             centerConeIndices = find(centerConnectivityVector > 0.0001);
    
             theConePoolingWeights = centerConnectivityVector(centerConeIndices);
             theConePositions = obj.inputConeMosaic.coneRFpositionsDegs(centerConeIndices,:);
             theConeCharacteristicRadii = coneCharacteristicRadiiDegs(centerConeIndices);
    
             obj.rfCenterVisualizationOutlines{iRGC} = subregionOutlineFromPooledCones(...
                 theConePositions, theConeCharacteristicRadii , theConePoolingWeights, ...
                 xSupport, ySupport, spatialSupportSamples);
    
        end

        fprintf('Done\n');
    end


    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1], 'Position', [10 10 850 800]);
        end
        ax = subplot('Position', [0.05 0.07 0.93 0.9]);
    end

    if (isempty(activation))
        % Visualize the RF centers
        faceAlpha = 0.5;
        lineStyle = '-';
        lineWidth = 1.0;
        cMap = [0.7 0.85 0.99; 0 0 0];
   
        for iRGC = 1:obj.rgcsNum
            rfVisualizationStruct = obj.rfCenterVisualizationOutlines{iRGC};
    
            % Plot contour of input source RFs
            transparentContourPlot(ax, ...
                rfVisualizationStruct.spatialSupportXY, ...
                rfVisualizationStruct.rfMap, ...
                rfVisualizationStruct.zLevels, ...
                faceAlpha, cMap, lineStyle, lineWidth);
        end
    else
        % Visualize the RF centers, color-coded with their activation level
        if (isempty(activationColorMap))
            cMap = gray(obj.rgcsNum);
        else
            cMap = activationColorMap;
        end
        
        if (isempty(activationRange))
            activationRange = [min(activation(:)) max(activation(:))];
        end

        activation = (activation - activationRange(1))/(activationRange(2)-activationRange(1));
        activation(activation<0) = 0;
        activation(activation>1) = 1;
        edgeColor = [0 1 0];
        faceAlpha = 1.0;
        lineWidth = 1.0;

        renderPatchArray(ax, obj.rfCenterVisualizationOutlines, activation, ...
            edgeColor, lineWidth, faceAlpha);
    end


    axis(ax, 'equal');

    if (~isempty(XLims))
        set(ax, 'XLim', XLims);
    end

    if (~isempty(YLims))
        set(ax, 'YLim', YLims);
    end

    set(ax, 'FontSize', 16);

    colormap(ax, cMap);
    if (isempty(backgroundColor))
        backgroundColor = [0.7 0.7 0.7];
    end
    set(ax, 'CLim', [0 1], 'Color', backgroundColor);

    box(ax, 'on')
    xlabel(ax, 'degrees'); 
    ylabel(ax, 'degrees');


end

function rfVisualizationStruct = subregionOutlineFromPooledCones(...
    conePos, coneRc, poolingWeights, ...
    xSupport, ySupport, spatialSupportSamples)

    % Compute spatial support
    xSep = max(coneRc)*3;
    if (isempty(xSupport))
        xx = conePos(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = conePos(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);

    RF = zeros(size(X));

    for iCone = 1:numel(poolingWeights)
        % Characteristic radius of the input RF
        rC = coneRc(iCone);
        % Compute aperture2D x weight
        XX = X-conePos(iCone,1);
        YY = Y-conePos(iCone,2);
        theAperture2D = poolingWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        % Accumulate 2D apertures
        RF = RF + theAperture2D;
    end

    zLevels(1) = 0.05*min(poolingWeights);
    zLevels(2) = max(poolingWeights);

    rfVisualizationStruct.spatialSupportXY = spatialSupportXY;
    rfVisualizationStruct.rfMap = RF;
    rfVisualizationStruct.zLevels = zLevels;
    rfVisualizationStruct.contourData = generateContourData(spatialSupportXY, RF, zLevels);

end

function transparentContourPlot(axesHandle, spatialSupportXY, zData, ...
                                zLevels, faceAlpha, cmap, lineStyle, lineWidth)

    xSupport = squeeze(spatialSupportXY(:,1));
    ySupport = squeeze(spatialSupportXY(:,2));
    C = contourc(xSupport, ySupport, zData, zLevels);
    dataPoints = size(C,2);
    startPoint = 1;
    hold(axesHandle, 'on');
    while (startPoint < dataPoints)
        theLevel = C(1,startPoint);
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);
        patch(axesHandle, 'Faces', f, 'Vertices', v, 'EdgeColor', [0 0 0], ...
            'FaceColor', cmap(max([1 round(theLevel*size(cmap,1))]),:), ...
            'FaceAlpha', faceAlpha, ... 
            'EdgeAlpha', 1, ...
            'LineStyle', lineStyle, 'LineWidth', lineWidth);
        startPoint = startPoint + theLevelVerticesNum+1;
    end

end

function cData = generateContourData(spatialSupportXY, zData, zLevels)

    xSupport = squeeze(spatialSupportXY(:,1));
    ySupport = squeeze(spatialSupportXY(:,2));
    C = contourc(xSupport, ySupport, zData, zLevels);
    dataPoints = size(C,2);
    startPoint = 1;

    iContour = 0;
    levelsNum = 0;
    while (startPoint < dataPoints)
        levelsNum = levelsNum + 1;
        if (levelsNum > 1)
        theLevel = C(1,startPoint)
        end
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);
        iContour = iContour+1;
        cData{iContour} = struct('faces', f, 'vertices', v);
        startPoint = startPoint + theLevelVerticesNum+1;
    end

end




function renderPatchArray(axesHandle, rfCenterVisualizationOutlines, faceColors, edgeColor, lineWidth, faceAlpha)


    rgcsNum = numel(rfCenterVisualizationOutlines);
    verticesNum = zeros(1, rgcsNum);
    for iRGC = 1:rgcsNum
        contourData = rfCenterVisualizationOutlines{iRGC}.contourData;
        s = contourData{1};
        verticesNum(iRGC) = size(s.vertices,1);
    end

    maxNoVertices = max(verticesNum);
    totalVerticesNum = sum(verticesNum);

    facesList = nan(rgcsNum, maxNoVertices);
    verticesList = zeros(totalVerticesNum,2);
    colors = zeros(totalVerticesNum,1);

    currentFacesNum = 0;
    for iRGC = 1:rgcsNum
        contourData = rfCenterVisualizationOutlines{iRGC}.contourData;
        s = contourData{1};
        newVerticesNum = size(s.vertices,1);
        idx = currentFacesNum+(1:newVerticesNum);
        verticesList(idx,:) =  s.vertices;
        colors(idx,:) = faceColors(iRGC);
        facesList(iRGC,1:newVerticesNum) = currentFacesNum + s.faces;
        currentFacesNum = currentFacesNum + newVerticesNum;
    end
    
    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.FaceAlpha = faceAlpha;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end
