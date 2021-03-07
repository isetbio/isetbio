function visualize(obj, varargin)
    p = inputParser;
    p.addParameter('domain', 'degrees', @(x)(ischar(x) && (ismember(x, {'degrees', 'microns'}))));
    p.addParameter('domainVisualizationLimits', [], @(x)((isempty(x))||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('visualizedConeAperture', 'lightCollectingArea', @(x)ismember(x, {'lightCollectingArea', 'geometricArea'}));
    p.addParameter('activation', []);
    p.addParameter('horizontalActivationSliceEccentricity', [], @(x)((isempty(x))||(isscalar(x))));
    p.addParameter('verticalActivationSliceEccentricity', [], @(x)((isempty(x))||(isscalar(x))));
    p.addParameter('activationRange', [],@(x)((isempty(x))||(numel(x)==2)));
    p.addParameter('activationColorMap', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.addParameter('horizontalActivationColorBar', false, @islogical);
    p.addParameter('verticalActivationColorBar', false, @islogical);
    p.addParameter('colorBarTickLabelPostFix', '', @ischar);
    p.addParameter('displayedEyeMovementData', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('currentEMposition', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('crossHairsOnMosaicCenter', false, @islogical);
    p.addParameter('crossHairsOnFovea', false, @islogical);
    p.addParameter('crossHairsOnOpticalImageCenter', false, @islogical);
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('backgroundColor', [0.7 0.7 0.7]);
    p.addParameter('plotTitle', '', @ischar);
    p.parse(varargin{:});
    
    domain = p.Results.domain;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    visualizedConeAperture = p.Results.visualizedConeAperture;
    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    activation = p.Results.activation;
    activationRange = p.Results.activationRange;
    currentEMposition = p.Results.currentEMposition;
    crossHairsOnMosaicCenter = p.Results.crossHairsOnMosaicCenter;
    crossHairsOnOpticalImageCenter = p.Results.crossHairsOnOpticalImageCenter;
    crossHairsOnFovea = p.Results.crossHairsOnFovea;
    displayedEyeMovementData = p.Results.displayedEyeMovementData;
    fontSize = p.Results.fontSize;
    cMap = p.Results.activationColorMap;
    verticalColorBar = p.Results.verticalActivationColorBar;
    horizontalColorBar = p.Results.horizontalActivationColorBar;
    colorBarTickLabelPostFix = p.Results.colorBarTickLabelPostFix;
    horizontalActivationSliceEccentricity = p.Results.horizontalActivationSliceEccentricity;
    verticalActivationSliceEccentricity = p.Results.verticalActivationSliceEccentricity;
    backgroundColor = p.Results.backgroundColor;
    plotTitle = p.Results.plotTitle;

    % Determine what eye movement data have to be displayed
    if (isstruct(displayedEyeMovementData))
       if (ischar(displayedEyeMovementData.trial))
          switch(displayedEyeMovementData.trial)
              case 'all'
                    displayedTrials = 1:size(obj.fixEMobj.emPosArcMin,1);
              otherwise
                error('unknown ''displayedEyeMovementData.trial'': ''%s''.', displayedEyeMovementData.trial);
          end
       else
          displayedTrials = displayedEyeMovementData.trial;
       end
       if (ischar(displayedEyeMovementData.timePoints))
           switch(displayedEyeMovementData.timePoints)
              case 'all'
                    displayedTimePoints = 1:size(obj.fixEMobj.emPosArcMin,2);
              otherwise
                error('unknown ''displayedEyeMovementData.samples'': ''%s''.', displayedEyeMovementData.samples);
          end
       else
           displayedTimePoints = displayedEyeMovementData.timePoints;
       end
    end
    
    % Determine displayed domain (degs or microns)
    switch (domain)
        case 'degrees'
            rfPositions = obj.coneRFpositionsDegs;
            rfSpacings = obj.coneRFspacingsDegs;
            rfProximityThreshold = 1/270;
            if (isstruct(displayedEyeMovementData))
                emPath = 1/60*obj.fixEMobj.emPosArcMin(displayedTrials,displayedTimePoints,:);
            else
                emPath = [];
            end
        case 'microns'
            rfPositions = obj.coneRFpositionsMicrons;
            rfSpacings = obj.coneRFspacingsMicrons;
            rfProximityThreshold = 1;
            if (isstruct(displayedEyeMovementData))
                emPath = obj.fixEMobj.emPosMicrons(displayedTrials,displayedTimePoints,:);
            else
                emPath = [];
            end
    end
    
    
    % Show the emPath on the center of the mosaic
    emPath = bsxfun(@plus, emPath, reshape(mean(rfPositions,1), [1 1 2]));
    
    % Determine X,Y limits
    if (isempty(domainVisualizationLimits))
        xRange(1) = min(rfPositions(:,1));
        xRange(2) = max(rfPositions(:,1));
        if (xRange(2) == xRange(1))
            switch (domain)
                case 'degrees'
                    xRange = xRange(1) + 0.02*[-1 1];
                case 'microns'
                    xRange = xRange(1) + 2*[-1 1];
            end
        end
        yRange(1) = min(rfPositions(:,2));
        yRange(2) = max(rfPositions(:,2));
        if (yRange(2) == yRange(1))
            switch (domain)
                case 'degrees'
                    yRange = yRange(1) + 0.02*[-1 1];
                case 'microns'
                    yRange = yRange(1) + 2*[-1 1];
            end
        end
        xx = xRange(2)-xRange(1);
        yy = yRange(2)-yRange(1);
        xRange(1) = xRange(1)-xx*0.02;
        xRange(2) = xRange(2)+xx*0.02;
        yRange(1) = yRange(1)-yy*0.02;
        yRange(2) = yRange(2)+yy*0.02;
    else
        xRange(1) = domainVisualizationLimits(1);
        xRange(2) = domainVisualizationLimits(2);
        yRange(1) = domainVisualizationLimits(3);
        yRange(2) = domainVisualizationLimits(4);
    end

    
    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.09 0.07 0.90 0.92]);
    else
        %figure(figureHandle);
        if (isempty(axesHandle))
            clf;
            set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.07 0.07 0.92 0.92]);
        end
    end
    
    % Number of cones
    conesNum = numel(rfSpacings);
    
    % Aperture shape (disk)
    if (conesNum > 10000)
        deltaAngle = 60;
    elseif (conesNum > 5000)
        deltaAngle = 45;
    elseif (conesNum > 1000)
       deltaAngle = 30;
    elseif (conesNum > 500)
        deltaAngle = 20;
    elseif (conesNum > 250)
        deltaAngle = 15;
    elseif (conesNum > 100)
        deltaAngle = 10;
    else
        deltaAngle = 5;
    end
    iTheta = (0:deltaAngle:360) / 180 * pi;
    coneApertureShape.x = cos(iTheta);
    coneApertureShape.y = sin(iTheta);
   
    
    cla(axesHandle);
    hold(axesHandle, 'on');
    
    % Visualize cone aperture multiplier
    if (strcmp(visualizedConeAperture', 'geometricArea'))
       visualizeApertureMultiplier = 1.0;
    else
       visualizeApertureMultiplier = obj.coneApertureToDiameterRatio;
    end
        
    if (isempty(activation))
        % Visualize cone types
        % Plot L-cones

        renderPatchArray(axesHandle, coneApertureShape, visualizeApertureMultiplier*rfSpacings(obj.lConeIndices)*0.5, ...
            rfPositions(obj.lConeIndices,:), 1/4*0.9, 'none', 1.0);
   
        % Plot M-cones
        renderPatchArray(axesHandle, coneApertureShape, visualizeApertureMultiplier*rfSpacings(obj.mConeIndices)*0.5, ...
            rfPositions(obj.mConeIndices,:), 2/4*0.9, 'none', 1.0);
        % Plot S-cones
        renderPatchArray(axesHandle, coneApertureShape, visualizeApertureMultiplier*rfSpacings(obj.sConeIndices)*0.5, ...
            rfPositions(obj.sConeIndices,:), 3/4*0.9, 'none', 1.0);
        % Plot K-cones
        renderPatchArray(axesHandle, coneApertureShape, visualizeApertureMultiplier*rfSpacings(obj.kConeIndices)*0.5, ...
            rfPositions(obj.kConeIndices,:), 4/4*0.9, 'none', 1.0);
    else
        if (isempty(activationRange))
            activationRange(1) = min(activation(:));
            activationRange(2) = max(activation(:));
        end
        if (activationRange(1) == activationRange(2))
            activationRange(1) = activationRange(1)-0.1;
            activationRange(2) = activationRange(2)+0.1;
        end
            
        activation = (activation - activationRange(1))/(activationRange(2)-activationRange(1));
        activation(activation<0) = 0;
        activation(activation>1) = 1;
        
        % Visualize activations
        % Plot L-cone activations
        renderPatchArray(axesHandle, coneApertureShape, visualizeApertureMultiplier*rfSpacings(obj.lConeIndices)*0.5, ...
            rfPositions(obj.lConeIndices,:), activation(obj.lConeIndices), [0 0 0], 0.1);
        % Plot M-cone activations
        renderPatchArray(axesHandle, coneApertureShape, visualizeApertureMultiplier*rfSpacings(obj.mConeIndices)*0.5, ...
            rfPositions(obj.mConeIndices,:), activation(obj.mConeIndices), [0 0 0], 0.1);
        % Plot S-cone activations
        renderPatchArray(axesHandle, coneApertureShape, visualizeApertureMultiplier*rfSpacings(obj.sConeIndices)*0.5, ...
            rfPositions(obj.sConeIndices,:), activation(obj.sConeIndices), [0 0 0], 0.1);
        % Plot K-cone activations
        renderPatchArray(axesHandle, coneApertureShape, visualizeApertureMultiplier*rfSpacings(obj.kConeIndices)*0.5, ...
            rfPositions(obj.kConeIndices,:), activation(obj.kConeIndices), [0 0 0], 0.1);
        
        if (~isempty(verticalActivationSliceEccentricity))
            d = abs(rfPositions(:,2)-verticalActivationSliceEccentricity);
            idx = find(d < rfProximityThreshold);
            activationSlice = squeeze(activation(idx));
            h = stem(rfPositions(idx,1), rfPositions(idx,2) + activationSlice*0.2*(yRange(2)-yRange(1)), ...
                'filled',  'g-', 'LineWidth', 1.5);
            h.BaseValue = verticalActivationSliceEccentricity;
        end
        
    end
    
    % Add crosshairs
    if (crossHairsOnMosaicCenter) || (crossHairsOnOpticalImageCenter) || (crossHairsOnFovea)
        if (isempty('withActivation'))
            crossHairsColor = [0 0 0];
        else
            crossHairsColor = [0 0 0];
        end
        
        if (crossHairsOnMosaicCenter)
            % Crosshairs centered on the middle of the mosaic
            xx1 = [xRange(1) xRange(2)];
            yy1 = mean(yRange)*[1 1];
            xx2 = mean(xRange)*[1 1];
            yy2 = [yRange(1) yRange(2)];
        elseif (crossHairsOnFovea)
            % Crosshairs centered on [0 0]
            xx1 = [xRange(1) xRange(2)];
            yy1 = [0 0];
            xx2 = [0 0];
            yy2 = [yRange(1) yRange(2)];
        else
            % Crosshairs centered on the middle of the optical image, i.e.,
            % (0,0) or current em position, if that is passed in
            xx1 = [xRange(1) xRange(2)];
            if (~isempty(currentEMposition))
                yy1 = -[1 1] * currentEMposition(2);
                xx2 = -[1 1] * currentEMposition(1);
            else
                yy1 = [0 0];
                xx2 = [0 0];
            end
            yy2 = [yRange(1) yRange(2)];
        end
        plot(axesHandle, xx1, yy1, '-', 'Color', crossHairsColor, 'LineWidth', 1.5);
        plot(axesHandle, xx2, yy2,  '-', 'Color', crossHairsColor,'LineWidth', 1.5);
    end
    
    % Superimpose eye movement path(s)
    if (~isempty(emPath))
        assert(size(emPath,3) == 2, sprintf('The third dimension of an emPath must be 2.'));
        nTrials = size(emPath,1);
        lineColors = brewermap(nTrials, 'Spectral');
        for emTrialIndex = 1:nTrials
            plot(axesHandle, emPath(emTrialIndex,:,1), emPath(emTrialIndex,:,2), ...
                'k-', 'LineWidth', 4.0);
        end
        for emTrialIndex = 1:nTrials
            plot(axesHandle, emPath(emTrialIndex,:,1), emPath(emTrialIndex,:,2), ...
                '-', 'LineWidth', 2, 'Color', squeeze(lineColors(emTrialIndex,:)));
        end
    end
    
    hold(axesHandle, 'off');
    
    % Set appropriate colormap
    if (isempty(activation))
        % Colormap for visualization of cone types
        cMap = [obj.lConeColor; obj.mConeColor; obj.sConeColor; obj.kConeColor];
    else
        % Colormap for visualization of activity
        if (isempty(cMap))
            cMap = gray(numel(obj.coneRFspacingsDegs)); %brewermap(numel(obj.coneRFspacingsDegs), '*greys');
        end

        if (ischar(backgroundColor) && strcmp(backgroundColor, 'mean of color map'))
            midRow = round(size(cMap,1)/2);
            backgroundColor = squeeze(cMap(midRow,:));
        else
            backgroundColor = squeeze(cMap(1,:));
        end
    end
    colormap(axesHandle, cMap);
    
    if (verticalColorBar) || (horizontalColorBar)
        colorBarTicks = [0 0.25 0.5 0.75 1];
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

        if (verticalColorBar)
            colorbar('eastOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels);
        elseif (horizontalColorBar)
            colorbar('northOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels);
        end
    end
    
    % Finalize plot
    set(axesHandle, 'Color', backgroundColor);
    axis(axesHandle, 'xy');
    axis(axesHandle, 'equal');
    set(axesHandle, 'XLim', xRange, 'YLim', yRange, 'CLim', [0 1], 'FontSize', fontSize);        
    
    if (isempty(domainVisualizationTicks))
        xo = (xRange(1)+xRange(2))/2;
        xx = xRange(2)-xRange(1);
        yo = (yRange(1)+yRange(2))/2;
        yy = yRange(2)-yRange(1);
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
    
    set(axesHandle, 'XTick', domainVisualizationTicks.x, ...
                    'YTick', domainVisualizationTicks.y);
     
    box(axesHandle, 'on');
    set(figureHandle, 'Color', [1 1 1]);
    
    switch (domain)
        case 'degrees'
            xlabel(axesHandle, 'space (degrees)');
            ylabel(axesHandle, 'space (degrees)');
        case 'microns'
            xlabel(axesHandle, 'space (microns)');
            ylabel(axesHandle, 'space (microns)');
    end
    if (isempty(plotTitle))
        title(axesHandle,sprintf('L (%2.1f%%), M (%2.1f%%), S (%2.1f%%), K (%2.1f%%), N = %d', ...
            100*obj.coneDensities(1), ...
            100*obj.coneDensities(2), ...
            100*obj.coneDensities(3), ...
            100*obj.coneDensities(4), ...
            conesNum));
    else
        title(axesHandle,plotTitle);
    end
    drawnow;
end

function renderPatchArray(axesHandle, apertureShape, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth)

    conesNum = numel(apertureRadii);
    if (conesNum == 0)
        return;
    end
    
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
        end
        facesList = cat(1, facesList, idx);
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.FaceAlpha = 1.0;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end
