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
    
    % Determin X,Y limits
    xRange(1) = min(rfPositions(:,1));
    xRange(2) = max(rfPositions(:,1));
    yRange(1) = min(rfPositions(:,2));
    yRange(2) = max(rfPositions(:,2));
    xx = xRange(2)-xRange(1);
    yy = yRange(2)-yRange(1);
    xRange(1) = xRange(1)-xx*0.05;
    xRange(2) = xRange(2)+xx*0.05;
    yRange(1) = yRange(1)-yy*0.05;
    yRange(2) = yRange(2)+yy*0.05;
    
    % Set Figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.07 0.07 0.92 0.92]);
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
        
    
    hold(axesHandle, 'on');
    
    if (isempty(activation))
        % Visualize cone types
        % Plot L-cones
        obj.renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.lConeIndices)*0.5, ...
            rfPositions(obj.lConeIndices,:), 1/4*0.9, [0 0 0], 1.0);
        % Plot M-cones
        obj.renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.mConeIndices)*0.5, ...
            rfPositions(obj.mConeIndices,:), 2/4*0.9, [0 0 0], 1.0);
        % Plot S-cones
        obj.renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.sConeIndices)*0.5, ...
            rfPositions(obj.sConeIndices,:), 3/4*0.9, [0 0 0], 1.0);
        % Plot K-cones
        obj.renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.kConeIndices)*0.5, ...
            rfPositions(obj.kConeIndices,:), 4/4*0.9, [0 0 0], 1.0);
    else
        % Visualize activations
        % Plot L-cone activations
        obj.renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.lConeIndices)*0.5, ...
            rfPositions(obj.lConeIndices,:), activation(obj.lConeIndices), [0 0 0], 1.0);
        % Plot M-cone activations
        obj.renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.mConeIndices)*0.5, ...
            rfPositions(obj.mConeIndices,:), activation(obj.mConeIndices), [0 0 0], 1.0);
        % Plot S-cone activations
        obj.renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.sConeIndices)*0.5, ...
            rfPositions(obj.sConeIndices,:), activation(obj.sConeIndices), [0 0 0], 1.0);
        % Plot K-cone activations
        obj.renderPatchArray(axesHandle, coneApertureShape, rfSpacings(obj.kConeIndices)*0.5, ...
            rfPositions(obj.kConeIndices,:), activation(obj.kConeIndices), [0 0 0], 1.0);
    end
    
    if (crossHairs)
        plot(axesHandle, [xRange(1) xRange(2)], mean(yRange)*[1 1], 'k-');
        plot(axesHandle, mean(xRange)*[1 1], [yRange(1) yRange(2)],  'k-');
    end
    
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

