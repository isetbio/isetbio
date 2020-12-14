function visualize(obj, varargin)
    p = inputParser;
    p.addParameter('domain', 'degs', @(x)(ischar(x) && (ismember(x, {'degs', 'microns'}))));
    p.addParameter('withActivation', []);
    p.addParameter('withCrossHairs', false, @islogical);
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('fontSize', 14, @isscalar);
    p.parse(varargin{:});
    
    domain = p.Results.domain;
    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    activation = p.Results.withActivation;
    crossHairs = p.Results.withCrossHairs;
    fontSize = p.Results.fontSize;
    
    switch (domain)
        case 'degs'
            rfPositions = obj.coneRFpositionsDegs;
            rfSpacings = obj.coneRFspacingsDegs;
        case 'microns'
            rfPositions = obj.coneRFpositionsMicrons;
            rfSpacings = obj.coneRFspacingsMicrons;
    end
    
    % Determine X,Y limits
    xRange(1) = min(rfPositions(:,1));
    xRange(2) = max(rfPositions(:,1));
    yRange(1) = min(rfPositions(:,2));
    yRange(2) = max(rfPositions(:,2));
    xx = xRange(2)-xRange(1);
    yy = yRange(2)-yRange(1);
    xRange(1) = xRange(1)-xx*0.02;
    xRange(2) = xRange(2)+xx*0.02;
    yRange(1) = yRange(1)-yy*0.02;
    yRange(2) = yRange(2)+yy*0.02;
    
    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.05 0.05 0.94 0.94]);
    else
        figure(figureHandle);
        if (isempty(axesHandle))
            clf;
            set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.07 0.07 0.92 0.92]);
        end
    end
     
    % Aperture shape (disk)
    iTheta = (0:15:360) / 180 * pi;
    coneApertureShape.x = cos(iTheta);
    coneApertureShape.y = sin(iTheta);
   
    cla(axesHandle);
    hold(axesHandle, 'on');
    if (isempty(activation))
        % Visualize cone types
        % Plot L-cones
        renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.lConeIndices)*0.5, ...
            rfPositions(obj.lConeIndices,:), 1/4*0.9, [0 0 0], 1.0);
        % Plot M-cones
        renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.mConeIndices)*0.5, ...
            rfPositions(obj.mConeIndices,:), 2/4*0.9, [0 0 0], 1.0);
        % Plot S-cones
        renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.sConeIndices)*0.5, ...
            rfPositions(obj.sConeIndices,:), 3/4*0.9, [0 0 0], 1.0);
        % Plot K-cones
        renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.kConeIndices)*0.5, ...
            rfPositions(obj.kConeIndices,:), 4/4*0.9, [0 0 0], 1.0);
    else
        % Visualize activations
        % Plot L-cone activations
        renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.lConeIndices)*0.5, ...
            rfPositions(obj.lConeIndices,:), activation(obj.lConeIndices), [0 0 0], 1.0);
        % Plot M-cone activations
        renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.mConeIndices)*0.5, ...
            rfPositions(obj.mConeIndices,:), activation(obj.mConeIndices), [0 0 0], 1.0);
        % Plot S-cone activations
        renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.sConeIndices)*0.5, ...
            rfPositions(obj.sConeIndices,:), activation(obj.sConeIndices), [0 0 0], 1.0);
        % Plot K-cone activations
        renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.kConeIndices)*0.5, ...
            rfPositions(obj.kConeIndices,:), activation(obj.kConeIndices), [0 0 0], 1.0);
    end
    
    if (crossHairs)
        plot(axesHandle, [xRange(1) xRange(2)], mean(yRange)*[1 1], 'k-');
        plot(axesHandle, mean(xRange)*[1 1], [yRange(1) yRange(2)],  'k-');
    end
    
    hold(axesHandle, 'off');
    
    % Set appropriate colormap
    if (isempty(activation))
        % Colormap for visualization of cone types
        cMap = [obj.lConeColor; obj.mConeColor; obj.sConeColor; obj.kConeColor];
    else
        % Colormap for visualization of activity
        cMap = gray(numel(obj.coneRFspacingsDegs));
    end
    colormap(axesHandle, cMap);
    
    % Finalize plot
    axis(axesHandle, 'equal');
    set(axesHandle, 'XLim', xRange, 'YLim', yRange, 'CLim', [0 1], 'FontSize', fontSize);
    box(axesHandle, 'on');
    
    switch (domain)
        case 'degs'
            xlabel(axesHandle, 'space (degs)');
            ylabel(axesHandle, 'space (degs)');
        case 'microns'
            xlabel(axesHandle, 'space (microns)');
            ylabel(axesHandle, 'space (degs)');
    end
    title(sprintf('L (%2.1f%%), M (%2.1f%%), S (%2.1f%%), K (%2.1f%%)\n', ...
        100*obj.coneDensities(1), ...
        100*obj.coneDensities(2), ...
        100*obj.coneDensities(3), ...
        100*obj.coneDensities(4)));
end

function renderPatchArray(axesHandle, apertureShape, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth)

    conesNum = numel(apertureRadii);
    verticesPerCone = numel(apertureShape.x);
    
    verticesList = zeros(verticesPerCone * conesNum, 2);
    facesList = [];
    
    if (numel(faceColors) == 1)
        colors = repmat(faceColors, [verticesPerCone*conesNum 1]);
    else
        colors = [];
    end
    
    for coneIndex = 1:conesNum
        idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
        verticesList(idx, 1) = apertureShape.x*apertureRadii(coneIndex) + rfCoords(coneIndex,1);
        verticesList(idx, 2) = apertureShape.y*apertureRadii(coneIndex) + rfCoords(coneIndex,2);
        if (numel(faceColors) == conesNum)
            colors = cat(1, colors, repmat(faceColors(coneIndex), [verticesPerCone 1]));
            size(colors)
        end
        facesList = cat(1, facesList, idx);
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.FaceAlpha = 0.7;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end
