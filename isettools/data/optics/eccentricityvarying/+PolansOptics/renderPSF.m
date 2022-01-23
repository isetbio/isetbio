function renderPSF(axesHandle, xSupport, ySupport, thePSF, xyRange, zLevels, cMap, contourLineColor, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('superimposedConeMosaic', [], @(x)(isempty(x)||isa(x, 'coneMosaicHex')));
    p.addParameter('withConeData', [],  @(x)(isempty(x)||isstruct(x)));
    p.addParameter('fontSize', 14, @isnumeric);
    p.addParameter('plotTitle', '', @ischar);
    p.addParameter('xyTicks', -20:2:20, @isnumeric);
    p.addParameter('alpha', 0.3, @isscalar);
    
    p.parse(varargin{:});
    theConeMosaic = p.Results.superimposedConeMosaic;
    theConeData = p.Results.withConeData;
    xyTicks = p.Results.xyTicks;
    fontSize = p.Results.fontSize;
    plotTitle = p.Results.plotTitle;
    alpha = p.Results.alpha;
    
    image(axesHandle, xSupport, ySupport, ones(size(thePSF,1), size(thePSF,2),3));
    axis(axesHandle, 'image');
    hold(axesHandle, 'on');
    
    smallTitle = true;
    if (~isempty(theConeData))
        conePositionsArcMin = theConeData.conePositionsArcMin;
        coneAperturesArcMin = theConeData.coneAperturesArcMin;
        faceColors = [0.9 0.9 0.9];
        faceAlpha = 0.32;
        lineWidth = 1.0;
        edgeColor = [0.6 0.6 0.6];
        renderPatchArray(axesHandle, coneAperturesArcMin, conePositionsArcMin, ...
    faceColors, edgeColor, lineWidth, faceAlpha);
        smallTitle = true;

    elseif (~isempty(theConeMosaic))
        % Retrieve cone positions (microns), cone spacings, and cone types
        cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
        conePositionsArcMin = cmStruct.coneLocs * 60;
        coneAperturesArcMin = cmStruct.coneApertures * 60;
        renderConeApertures(axesHandle, conePositionsArcMin, coneAperturesArcMin, [.2 0.2 0.2]);
        smallTitle = false;
    end
    

    % Render the semi-transparent plot of the RGC RF
    if (isempty(cMap))
        cMap = brewermap(1024, 'reds');
    end
    semiTransparentContourPlot(axesHandle, xSupport, ySupport, thePSF, zLevels, cMap, alpha, contourLineColor);
    imagesc(axesHandle, xSupport, ySupport, thePSF);
    colormap(axesHandle, brewermap(1024, '*spectral'));
    hold(axesHandle, 'on');
    
    xlabel(axesHandle, 'arc min');
    ylabel(axesHandle, 'arc min');
    grid(axesHandle, 'on');
  
    box(axesHandle, 'on');
    set(axesHandle, 'Color', [0.8 0.85 0.85]);
    set(axesHandle, 'XLim', xyRange, 'YLim', xyRange, 'CLim', [0 1], 'Color', [1 1 1]);
    set(axesHandle, 'XTick', xyTicks, 'YTick', xyTicks);
    set(axesHandle, 'FontSize', fontSize);
    if (~isempty(plotTitle))
        if (smallTitle)
            xx = xyRange(1) + 0.1*(xyRange(2)-xyRange(1));
            yy = xyRange(2) - 0.1*(xyRange(2)-xyRange(1));
            title(axesHandle,plotTitle, 'FontSize', 12, 'FontWeight', 'normal', 'Color', [.5 .2 1], 'FontName', 'Menlo');
        else
            title(axesHandle,plotTitle);
        end
    end

end

function renderConeApertures(axesHandle, conePositionsArcMin, coneAperturesArcMin, color)
  
    xOutline = cosd(0:30:360);
    yOutline = sind(0:30:360);
    for iCone = 1:size(conePositionsArcMin,1)
        r = 0.5*coneAperturesArcMin(iCone);
        plot(axesHandle, xOutline*r +  conePositionsArcMin(iCone,1), ...
                         yOutline*r +  conePositionsArcMin(iCone,2), ...
                         'k-', 'Color', color, 'LineWidth', 1.0);
         if (iCone > 1)
             hold(axesHandle, 'on');
         end
    end
                 
end


function renderPatchArray(axesHandle, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth, faceAlpha)

    apertureShape.x = cosd(0:30:360);
    apertureShape.y = sind(0:30:360);
    
    conesNum = numel(apertureRadii);
    if (conesNum == 0)
        return;
    end
    
    verticesPerCone = numel(apertureShape.x);
    
    verticesList = zeros(verticesPerCone * conesNum, 2);
    facesList = [];
    
    if (numel(faceColors) == 3)
        colors = repmat(faceColors, [verticesPerCone*conesNum 1]);
    else
        colors = [];
    end
    
    for coneIndex = 1:conesNum
        idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
        verticesList(idx, 1) = apertureShape.x*0.5*apertureRadii(coneIndex) + rfCoords(coneIndex,1);
        verticesList(idx, 2) = apertureShape.y*0.5*apertureRadii(coneIndex) + rfCoords(coneIndex,2);
        if (numel(faceColors) == conesNum)
            colors = cat(1, colors, repmat(faceColors(coneIndex), [verticesPerCone 1]));
        end
        facesList = cat(1, facesList, idx);
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.FaceAlpha = faceAlpha;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
    hold(axesHandle, 'on');
end



 
% Function to generate a semitransparent controur plot
function semiTransparentContourPlot(axesHandle, xSupport, ySupport, zData, zLevels, cmap, alpha, contourLineColor)
    % Compute contours at desired Z-level
    C = contourc(xSupport, ySupport, zData, zLevels);
    % Go through the contour matrix and plot each contour separately
    dataPoints = size(C,2);
    startPoint = 1;
    hold(axesHandle, 'on');
    while (startPoint < dataPoints)
        theLevel = C(1,startPoint);
        theNormalizedLevel = (theLevel-min(zLevels))/(max(zLevels)-min(zLevels));
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);

        lutIndex = 1+round(theNormalizedLevel*(size(cmap,1)-1));
        patch('Faces', f, 'Vertices', v, 'EdgeColor', [0 0 0], ...
            'FaceColor', cmap(lutIndex,:), ...
            'FaceAlpha', alpha, ...
            'EdgeColor', contourLineColor, ...
            'LineStyle', '-', 'LineWidth', 1.0, ...
        'Parent', axesHandle);
        startPoint = startPoint + theLevelVerticesNum+1;
    end
end
 
