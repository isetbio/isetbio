function visualize(obj, varargin)
    
    % Parse optional input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('component', 'RF centers', @(x)ismember(x, {'RF centers'}));
    p.addParameter('activation', []);
    p.addParameter('identifyInputCones', true, @islogical);
    p.addParameter('activationRange', [],@(x)((isempty(x))||(numel(x)==2)));
    p.addParameter('activationColorMap', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.addParameter('horizontalActivationColorBar', false, @islogical);
    p.addParameter('verticalActivationColorBar', false, @islogical);
    p.addParameter('horizontalActivationColorBarInside', false, @islogical);
    p.addParameter('verticalActivationColorBarInside', false, @islogical);
    p.addParameter('fontSize', 16, @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('labelRGCsWithIndices', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('colorbarFontSize', 16, @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('colorBarTickLabelPostFix', '', @ischar);
    p.addParameter('colorbarTickLabelColor',  [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 3))));
    p.addParameter('backgroundColor', [0 0 0], @(x)((ischar(x)&&(strcmp(x,'none')))||isempty(x)||((isvector(x))&&(numel(x) == 3))));
    p.addParameter('plotTitle', '', @(x)(isempty(x) || ischar(x) || islogical(x)));
    p.addParameter('plotTitleColor', [0 0 0], @isnumeric);
    p.addParameter('domainVisualizationLimits', [], @(x)((isempty(x))||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x)&&((isfield(x, 'x'))&&(isfield(x,'y'))))));
   
    
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    visualizedComponent = p.Results.component;
    identifyInputCones = p.Results.identifyInputCones;
    activation = p.Results.activation;
    activationRange = p.Results.activationRange;
    activationColorMap = p.Results.activationColorMap;
    colorBarTickLabelPostFix = p.Results.colorBarTickLabelPostFix;
    verticalColorBar = p.Results.verticalActivationColorBar;
    horizontalColorBar = p.Results.horizontalActivationColorBar;
    verticalColorBarInside = p.Results.verticalActivationColorBarInside;
    horizontalColorBarInside = p.Results.horizontalActivationColorBarInside;
    colorbarTickLabelColor = p.Results.colorbarTickLabelColor;
    backgroundColor = p.Results.backgroundColor;
    fontSize = p.Results.fontSize;
    labelRGCsWithIndices = p.Results.labelRGCsWithIndices;
    colorbarFontSize = p.Results.colorbarFontSize;
    plotTitle = p.Results.plotTitle;
    plotTitleColor = p.Results.plotTitleColor;


    % Generate the visualization cache
    xSupport = [];
    ySupport = []; 
    spatialSupportSamples = 60;
    generateVisualizationCache(obj, xSupport, ySupport, spatialSupportSamples);

    % Determine X,Y limits
    if (isempty(domainVisualizationLimits))
        xRange = obj.visualizationCache.surroundConesXrange;
        if (xRange(2) == xRange(1))
            xRange = xRange(1) + 0.02*[-1 1];
        end
        yRange = obj.visualizationCache.surroundConesYrange;
        if (yRange(2) == yRange(1))
            yRange = yRange(1) + 0.02*[-1 1];
        end
        xx = xRange(2)-xRange(1);
        yy = yRange(2)-yRange(1);
        XLims(1) = xRange(1)-xx*0.02;
        XLims(2) = xRange(2)+xx*0.02;
        YLims(1) = yRange(1)-yy*0.02;
        YLims(2) = yRange(2)+yy*0.02;
        domainVisualizationLimits(1) = XLims(1);
        domainVisualizationLimits(2) = XLims(2);
        domainVisualizationLimits(3) = YLims(1);
        domainVisualizationLimits(4) = YLims(2);
    else
        XLims(1) = domainVisualizationLimits(1);
        XLims(2) = domainVisualizationLimits(2);
        YLims(1) = domainVisualizationLimits(3);
        YLims(2) = domainVisualizationLimits(4);
    end
    

    if (isempty(domainVisualizationTicks))
        xo = (XLims(1)+XLims(2))/2;
        xx = XLims(2)-XLims(1);
        yo = (YLims(1)+YLims(2))/2;
        yy = YLims(2)-YLims(1);
        ticksX = xo + xx*0.5*[-0.75 0 0.75];
        ticksY = yo + yy*0.5*[-0.75 0 0.75];
        
        if (xx > 10)
            domainVisualizationTicks.x = round(ticksX);
        elseif (xx > 1)
            domainVisualizationTicks.x = round(ticksX*100)/100;
        else
            domainVisualizationTicks.x = round(ticksX*1000)/1000;
        end
        if (yy > 10)
            domainVisualizationTicks.y = round(ticksY);
        elseif (yy > 1)
            domainVisualizationTicks.y = round(ticksY*100)/100;
        else
            domainVisualizationTicks.y = round(ticksY*1000)/1000;
        end
        
    end

    switch (visualizedComponent)
        case 'RF centers'
            [hFig, ax] = visualizeRFcenters(obj, hFig, ax, ...
                XLims, YLims, domainVisualizationTicks, ...
                identifyInputCones, labelRGCsWithIndices, ...
                activation, activationRange, activationColorMap, ...
                colorBarTickLabelPostFix, colorbarTickLabelColor, ...
                verticalColorBar, horizontalColorBar, colorbarFontSize, ...
                verticalColorBarInside, horizontalColorBarInside, ...
                backgroundColor, fontSize, ...
                plotTitle, plotTitleColor);

        otherwise
            error('Uknown visualized component: ''%s''.', visualizedComponent);
    end

end


function [hFig, ax] = visualizeRFcenters(obj,hFig, ax, XLims, YLims, domainVisualizationTicks, ...
    identifyInputCones, labelRGCsWithIndices, ...
    activation, activationRange, activationColorMap, ...
    colorBarTickLabelPostFix, colorbarTickLabelColor, ...
    verticalColorBar, horizontalColorBar, colorbarFontSize, ...
    verticalColorBarInside, horizontalColorBarInside, ...
    backgroundColor, fontSize, plotTitle, plotTitleColor)

    
    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1], 'Position', [10 10 1120 1050], 'Name', obj.name);
        end
        ax = subplot('Position', [0.05 0.05 0.93 0.93]);
    end

    
    if (~isempty(activation))
        % Visualize RF centers, color-coded with their activation level
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

        currentFacesNum = 0;
        for iRGC = 1:obj.rgcsNum
            newVerticesNum = obj.visualizationCache.rfCenterPatchData.VerticesNumForRGC(iRGC);
            idx = currentFacesNum + (1:newVerticesNum);
            obj.visualizationCache.rfCenterPatchData.FaceVertexCData(idx,:) = activation(iRGC);
            currentFacesNum = currentFacesNum + newVerticesNum;
        end
    else
        % All mRGC centers in gray
        cMap = [0 0 0; 0.5 0.5 0.5; 0 0 0];
        S.FaceVertexCData = obj.visualizationCache.rfCenterPatchData.FaceVertexCData*0+0.5;
    end

    % Plot the RFs
    S.Vertices = obj.visualizationCache.rfCenterPatchData.Vertices;
    S.Faces = obj.visualizationCache.rfCenterPatchData.Faces;
    S.FaceVertexCData = obj.visualizationCache.rfCenterPatchData.FaceVertexCData;

    S.FaceColor = 'flat';
    if (~isempty(activation))
        S.EdgeColor = 'none';
    else
        S.EdgeColor = [1 1 1];
    end
    S.FaceAlpha = 1.0;
    S.LineWidth = 1.0;
    patch(S, 'Parent', ax)


    if (identifyInputCones)
        hold(ax, 'on')

        % Plot the connections from the RF center to the input L-cones
        idx = find(obj.visualizationCache.rfCenterConeConnectionLineSegments.coneTypes == cMosaic.LCONE_ID);
        plot(ax, ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
            'Color', [0 0 0],...
            'LineWidth', 3.0);  
        plot(ax, ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
            'Color', [1 0 0],...
            'LineWidth', 1.5); 

        idx = find(obj.visualizationCache.rfCenterConeConnectionLineSegments.coneTypes == cMosaic.MCONE_ID);
        plot(ax, ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
            'Color', [0 0 0],...
            'LineWidth', 3.0); 
        plot(ax, ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
            'Color', [0 1 0],...
            'LineWidth', 1.5); 

        % Plot the L-cones
        plot(ax, obj.inputConeMosaic.coneRFpositionsDegs(obj.visualizationCache.lConeIndicesConnectedToRGCsurrounds,1), ...
                 obj.inputConeMosaic.coneRFpositionsDegs(obj.visualizationCache.lConeIndicesConnectedToRGCsurrounds,2), ...
                 'r.', 'MarkerSize', 10);

        % Plot the M-cones
        plot(ax, obj.inputConeMosaic.coneRFpositionsDegs(obj.visualizationCache.mConeIndicesConnectedToRGCsurrounds ,1), ...
                 obj.inputConeMosaic.coneRFpositionsDegs(obj.visualizationCache.mConeIndicesConnectedToRGCsurrounds ,2), ...
                 'g.',  'MarkerSize', 10);

        % Plot the S-cones
        plot(ax, obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.sConeIndices,1), ...
                 obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.sConeIndices,2), ...
                 'b.',  'MarkerSize', 10);
    end

    if (~isempty(labelRGCsWithIndices))
        hold(ax, 'on')
        for iRGC = 1:numel(labelRGCsWithIndices)
            theRGCindex = labelRGCsWithIndices(iRGC);
            S = obj.visualizationCache.theContourData{theRGCindex}{1};
            S.FaceVertexCData = 0.5;
            S.FaceColor = 'flat';
            S.EdgeColor = [1 1 0];
            S.FaceAlpha = 0.0;
            S.LineWidth = 1.0;
            S.LineStyle = '-';
            patch(S, 'Parent', ax);

        end
    end


    % Finalize plot
    axis(ax, 'equal');
    if (~isempty(XLims))
        set(ax, 'XLim', XLims);
    end
    if (~isempty(YLims))
        set(ax, 'YLim', YLims);
    end
    set(ax, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);
   
    set(ax, 'FontSize', fontSize);

    colormap(ax, cMap);
    if (isempty(backgroundColor))
        set(ax, 'CLim', [0 1], 'Color', 'none');
    else
        set(ax, 'CLim', [0 1], 'Color', backgroundColor);
    end

    box(ax, 'on')
    xlabel(ax, 'degrees'); 
    ylabel(ax, 'degrees');

    % Colorbar and colorbar ticks
    if (~isempty(activation))
        if (verticalColorBar) || (horizontalColorBar) || (verticalColorBarInside) || (horizontalColorBarInside)
            colorBarTicks = [0.00 0.25 0.5 0.75 1.0];
            colorBarTickLabels = cell(1, numel(colorBarTicks));
            colorBarTickLevels = activationRange(1) + (activationRange(2)-activationRange(1)) * colorBarTicks;
            
            for k = 1:numel(colorBarTicks)
                if (max(abs(colorBarTickLevels)) >= 10)
                    colorBarTickLabels{k} = sprintf('%2.0f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                elseif (max(abs(colorBarTickLevels)) >= 1)
                    colorBarTickLabels{k} = sprintf('%2.1f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                elseif (max(abs(colorBarTickLevels)) >= 0.1)
                    colorBarTickLabels{k} = sprintf('%2.2f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                else
                    colorBarTickLabels{k} = sprintf('%2.3f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                end
            end
            
            if (isempty(colorbarFontSize))
                colorbarFontSize = fontSize/2;
            end
    
  
            if (verticalColorBar)
                colorbar(ax, 'eastOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor);
            elseif (verticalColorBarInside)
                colorbar(ax, 'east', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor,  'FontWeight', 'Bold', 'FontSize', colorbarFontSize, 'FontName', 'Spot mono');
            elseif (horizontalColorBar)
                colorbar(ax,'northOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor);
            elseif (horizontalColorBarInside)
                colorbar(ax,'north', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor,  'FontWeight', 'Bold', 'FontSize', colorbarFontSize, 'FontName', 'Spot mono');
            end
        else
            colorbar(ax, 'off');
        end
    else
        colorbar(ax, 'off');
    end

    if (plotTitle)    
       title(ax, plotTitle, 'Color', plotTitleColor);
    end


    fprintf('\nDrawning mRGCmosaic patch. Please wait ...');
    tic
    drawnow;
    fprintf('Done in %2.2f seconds\n', toc);

end


function generateVisualizationCache(obj, xSupport, ySupport, spatialSupportSamples)
    if (isfield(obj.visualizationCache, 'rfCenterPatchData')) && ...
       (~isempty(obj.visualizationCache.rfCenterPatchData))
        % Already in visualizationCache, so return
        return;
    end

    fprintf('\nComputing RF center outline contours. Please wait ...');
    tic

    coneCharacteristicRadiiDegs = ...
        obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
        obj.inputConeMosaic.coneApertureDiametersDegs;


    verticesNumForRGC = zeros(1, obj.rgcsNum);
    theContourData = cell(1, obj.rgcsNum);
 

    lineSegmentIndex = 0;
    for iRGC = 1:obj.rgcsNum
        % Retrieve the center cone indices & weights
        centerConnectivityVector = full(squeeze(obj.centerConePoolingMatrix(:, iRGC)));
        centerConeIndices = find(centerConnectivityVector > 0.0001);

        theConePoolingWeights = centerConnectivityVector(centerConeIndices);
        theConePositions = obj.inputConeMosaic.coneRFpositionsDegs(centerConeIndices,:);
        theConeCharacteristicRadii = coneCharacteristicRadiiDegs(centerConeIndices);
        theConeTypes = obj.inputConeMosaic.coneTypes(centerConeIndices);

        inputConesNum = size(theConePositions,1);
        if (inputConesNum > 1)
            for iCone = 1:inputConesNum
                lineSegmentIndex = lineSegmentIndex + 1;
                rfCenterConeConnectionLineSegments.Xpos(:, lineSegmentIndex) = ...
                    [obj.rgcRFpositionsDegs(iRGC,1) theConePositions(iCone,1)]';
                rfCenterConeConnectionLineSegments.Ypos(:, lineSegmentIndex) = ...
                    [obj.rgcRFpositionsDegs(iRGC,2) theConePositions(iCone,2)]';
                rfCenterConeConnectionLineSegments.coneTypes(lineSegmentIndex) = theConeTypes(iCone);
            end
        end


        theContourData{iRGC} = subregionOutlineContourFromPooledCones(...
                 theConePositions, theConeCharacteristicRadii , theConePoolingWeights, ...
                 xSupport, ySupport, spatialSupportSamples);

        s = theContourData{iRGC}{1};
        verticesNumForRGC(iRGC) = size(s.vertices,1);
    end

    maxNoVertices = max(verticesNumForRGC);
    totalVerticesNum = sum(verticesNumForRGC);

    facesList = nan(obj.rgcsNum, maxNoVertices);
    verticesList = zeros(totalVerticesNum,2);
    colorVertexCData = zeros(totalVerticesNum,1);

    currentFacesNum = 0;

    for iRGC = 1:obj.rgcsNum
        s = theContourData{iRGC}{1};
        newVerticesNum = size(s.vertices,1);
        idx = currentFacesNum+(1:newVerticesNum);
        verticesList(idx,:) =  s.vertices;
        colorVertexCData(idx,:) = 0.5;
        facesList(iRGC,1:newVerticesNum) = currentFacesNum + s.faces;
        currentFacesNum = currentFacesNum + newVerticesNum;
    end

    obj.visualizationCache.theContourData = theContourData;
    obj.visualizationCache.rfCenterConeConnectionLineSegments = rfCenterConeConnectionLineSegments;
    obj.visualizationCache.rfCenterPatchData.VerticesNumForRGC = verticesNumForRGC;
    obj.visualizationCache.rfCenterPatchData.Vertices = verticesList;
    obj.visualizationCache.rfCenterPatchData.Faces = facesList;
    obj.visualizationCache.rfCenterPatchData.FaceVertexCData = colorVertexCData;

    % Find all input cone indices that are connected to the RF surrounds
    surroundConnectedConeIndices = find(sum(obj.surroundConePoolingMatrix,2) > 0);
    xx = squeeze(obj.inputConeMosaic.coneRFpositionsDegs(surroundConnectedConeIndices,1));
    yy = squeeze(obj.inputConeMosaic.coneRFpositionsDegs(surroundConnectedConeIndices,2));
    obj.visualizationCache.surroundConesXrange = [min(xx) max(xx)];
    obj.visualizationCache.surroundConesYrange = [min(yy) max(yy)];
    idx = find(obj.inputConeMosaic.coneTypes(surroundConnectedConeIndices) == cMosaic.LCONE_ID);
    obj.visualizationCache.lConeIndicesConnectedToRGCsurrounds = surroundConnectedConeIndices(idx);
    idx = find(obj.inputConeMosaic.coneTypes(surroundConnectedConeIndices) == cMosaic.MCONE_ID);
    obj.visualizationCache.mConeIndicesConnectedToRGCsurrounds = surroundConnectedConeIndices(idx);

    fprintf(' Done in %2.1f seconds\n', toc);

end




function contourData = subregionOutlineContourFromPooledCones(...
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

    contourData = generateContourData(spatialSupportXY, RF, zLevels);
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