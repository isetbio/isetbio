function hFig = visualizeActivationMaps(obj, activation, varargin)
% Separately visualize mosaic activations for each submosaic and the whole
%
% Syntax:
%   hFig = visualizeActivationMaps(obj, activation, varargin)
%
% Description:
%    Visualize mosaic activations separately for each submosaic and for the
%    entire mosaic.
%
% Inputs:
%    obj        - The cone mosaic hex object
%    activation - The activation map
%    varargin   - (Optional) Additional information for creating the
%                 visualized activations
%
% Outputs:
%    hFig       - The figure handle
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/15  NPC  ISETBIO TEAM, 2015
%    02/16/18  jnm  Formatting

    p = inputParser;
    p.addParameter('figureSize', [920 875], @isnumeric);
    p.addParameter('signalName', ' ', @ischar);
    p.addParameter('mapType', 'modulated disks', @ischar);
    p.addParameter('visualizedConeAperture', 'lightCollectingArea', ...
        @(x)ismember(x, {'lightCollectingArea', 'geometricArea'}));
    p.addParameter('colorMap', jet(1024), @isnumeric);
    p.addParameter('signalRange', [], @isnumeric);
    p.addParameter('xRange', [], @isnumeric);
    p.addParameter('yRange', [], @isnumeric);
    p.addParameter('separateLMSmosaics', true, @islogical);
    p.addParameter('activationTime', [], @isnumeric);
    p.addParameter('zoomInFactor', 1.0, @isnumeric);
    p.addParameter('visualizedInstanceIndex', nan, @isnumeric);

    p.parse(varargin{:});

    if (any(size(activation) ~= size(obj.pattern)))
       activation = obj.reshapeHex2DmapToHex3Dmap(activation);
    end

    if strcmp(p.Results.mapType, 'modulated disks') || ...
            strcmp(p.Results.mapType, 'modulated hexagons')
        hFig = visualizeMosaicActivationsMapsAsModulatedPixels(obj, ...
            activation, p.Results.mapType, ...
            p.Results.visualizedConeAperture, p.Results.colorMap, ...
            p.Results.signalName, p.Results.signalRange, ...
            p.Results.xRange, p.Results.yRange, ...
            p.Results.separateLMSmosaics, p.Results.activationTime, ...
            p.Results.zoomInFactor, p.Results.visualizedInstanceIndex, ...
            p.Results.figureSize);
    elseif strcmp(p.Results.mapType, 'density plot')
        hFig = visualizeMosaicActivationsAsDensityMaps(obj, activation, ...
            p.Results.colorMap, p.Results.signalName, ...
            p.Results.signalRange,p.Results.figureSize);
    else
        error('visualizeActivationMaps:: Unknown map type');
    end
end

function hFig = visualizeMosaicActivationsMapsAsModulatedPixels(obj, ...
    activation, mapType, visualizedConeAperture, cMap, signalName, signalRange, xRange, yRange, ...
    separateLMSmosaics, activationTime, zoomInFactor, instanceIndex, ...
    figureSize)
% Visualize mosaic activations as disk mosaics
%
% Syntax:
%   hFig = visualizeMosaicActivationsMapsAsModulatedPixels(obj, ...
%       activation, mapType, cMap, signalName, signalRange, xRange, ...
%       yRange, separateLMSmosaics, activationTime, zoomInFactor, ...
%       instanceIndex, figureSize)
%
% Description:
%    Visualize mosaic activations as disk mosaics
%
% Inputs:
%    obj                - The cone mosaic hex object
%    activation         - The activation to map
%    mapType            - The type of map
%    cMap               - The color map
%    signalName         - The name of the signal
%    signalRange        - The range for the signal
%    xRange             - The x-axis range for the graph
%    yRange             - The y-axis range for the graph
%    separateLMSmosaics - Boolean indicating whether to plot the LMS
%                         mosaics separately.
%    activationTime     - The latency of the activation map, in ms
%    zoomInFactor       - The factor by which to zoom in (Multiply ranges)
%    instanceIndex      - The index for the instance (of the response)
%    figureSize         - The size of the figure
%
% Outputs:
%    hFig               - The figure handle
%
% Optional key/value pairs:
%    None.
%
% History:
%
%    5/26/20   NPC  Fixed iRow issue, which was causing the mosaic plotting Y-coord flip 
%                   (rows grow top -> bottom, whereas Y-coords grow bottom -> top)
%

    sampledHexMosaicXaxis = squeeze(obj.patternSupport(1, :, 1)) + ...
        obj.center(1);
    sampledHexMosaicYaxis = squeeze(obj.patternSupport(:, 1, 2)) + ...
        obj.center(2);

    if (strcmp(visualizedConeAperture, 'lightCollectingArea'))
        % Note that pigment.pdWidth defines the size of a square collective
        % aperture. Here we compute the equivalent circular aperture
    	aperture = diameterForCircularApertureFromWidthForSquareAperture(...
        obj.pigment.pdWidth);
        elseif (strcmp(visualizedConeAperture, 'geometricArea'))
        aperture = diameterForCircularApertureFromWidthForSquareAperture(...
        obj.pigment.width);
    end
    
    if strcmp(mapType, 'modulated disks') 
        iTheta = (0:15:360) / 180 * pi;
        apertureOutline.x = aperture / 2.0 * cos(iTheta);
        apertureOutline.y = aperture / 2.0 * sin(iTheta);
    elseif strcmp(mapType, 'modulated hexagons')
        iTheta = (0:60:360) / 180 * pi;
        apertureOutline.x = aperture / 2.0 * cos(iTheta);
        apertureOutline.y = aperture / 2.0 * sin(iTheta);
    end

    activeConesActivations = activation(obj.pattern > 1);
    if (isempty(signalRange))
        activationRange = [min(activeConesActivations(:)) ...
            max(activeConesActivations(:))];
    else
        activationRange = signalRange;
    end

    hFig = figure();
    clf;
    set(hFig, 'Position', [10 10 figureSize(1) figureSize(2)], ...
        'Color', [0 0 0]); % , 'MenuBar', 'None');

    if (separateLMSmosaics)
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 2, ...
               'colsNum', 3, ...
               'heightMargin', 0.01, ...
               'widthMargin', 0.01, ...
               'leftMargin', 0.01, ...
               'rightMargin', 0.006, ...
               'bottomMargin', 0.005, ...
               'topMargin', 0.005);

        for subplotIndex = 1:4
            if (subplotIndex < 4)
                subplot('Position', subplotPosVectors(1, subplotIndex).v);
            else
                subplot('Position', subplotPosVectors(2, 2).v);
            end
            set(gca, 'Color', [0 0 0]);
            switch subplotIndex
                case 1
                    idx = find(obj.pattern == 2);
                    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
                    subplotTitle = 'L-cone submosaic activation';
                    showXticks = false;
                    showYticks = false;
                    edgeColor = 'none';
                    lineWidth = 0.1;
                case 2
                    idx = find(obj.pattern == 3);
                    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
                    subplotTitle = 'M-cone submosaic activation';
                    showXticks = false;
                    showYticks = false;
                    edgeColor = 'none';
                    lineWidth = 0.1;
                case 3
                    idx = find(obj.pattern == 4);
                    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
                    subplotTitle = 'S-cone submosaic activation';
                    showXticks = false;
                    showYticks = false;
                    edgeColor = 'none';
                    lineWidth = 0.1;
                case 4
                    idx = find(obj.pattern > 1);
                    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
                    subplotTitle = 'total mosaic activation';
                    showXticks = true;
                    showYticks = true;
                    edgeColor = 'none';
                    lineWidth = 0.1;
            end % switch

            cMapLevels = size(cMap, 1);
            activationsNlevels = round((activation(idx) - ...
                activationRange(1)) / (activationRange(2) - ...
                activationRange(1)) * cMapLevels);
            faceColorsNormalizedValues = activationsNlevels / cMapLevels;
            renderPatchArray(apertureOutline, ...
                sampledHexMosaicXaxis(iCols), ...
                sampledHexMosaicYaxis(end-iRows+1), ...
                faceColorsNormalizedValues, edgeColor, lineWidth);
            set(gca, 'CLim', [0 1]);
            axis 'image';
            axis 'xy';
            xTicks = obj.center(1) + (-75:75:75);
            yTicks = obj.center(2) + (-75:75:75);
            xTickLabels = sprintf('%2.0f um\n', xTicks * 1e6);
            yTickLabels = sprintf('%2.0f um\n', yTicks * 1e6);
            set(gca, 'XTick', xTicks, 'YTick', yTicks, ...
                'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels, ...
                'XColor', [0.0 0.0 0.0], 'YColor', [0.0 0.0 0.0]);
            if (~showXticks), set(gca, 'XTick', []); end
            if (~showYticks), set(gca, 'YTick', []); end
            box on;
            grid off;
            if (isempty(xRange))
                xRange = [sampledHexMosaicXaxis(1) - obj.pigment.width ...
                    sampledHexMosaicXaxis(end) + obj.pigment.width];
            else
                xRange(1) = xRange(1) - obj.pigment.width;
                xRange(2) = xRange(2) + obj.pigment.width;
            end
            if (isempty(yRange))
                yRange = [sampledHexMosaicYaxis(1) - obj.pigment.width ...
                    sampledHexMosaicYaxis(end) + obj.pigment.width];
            else
                yRange(1) = yRange(1) - obj.pigment.width;
                yRange(2) = yRange(2) + obj.pigment.width;
            end

            set(gca, 'XLim', xRange * zoomInFactor);
            set(gca, 'YLim', yRange * zoomInFactor);
            set(gca, 'FontSize', 18, 'FontName', 'Menlo');
            title(subplotTitle, 'FontSize', 18, 'Color', [1 1 1], ...
                'FontName', 'Menlo');

            if (subplotIndex == 4)
                ticks = 0:0.1:1.0;
                tickLabels = sprintf('%2.0f\n', ...
                    round(activationRange(1) + ticks * ...
                    (activationRange(2) - activationRange(1))));
                % Add colorbar
                originalPosition = get(gca, 'position');
                hCbar = colorbar('eastoutside', 'peer', gca, ...
                    'Ticks', ticks, 'TickLabels', tickLabels);
                hCbar.Orientation = 'vertical';
                hCbar.Label.String = signalName;
                hCbar.FontSize = 16;
                hCbar.FontName = 'Menlo';
                hCbar.Color = [0.9 0.9 0.5];
                % The addition changes the figure size, so undo this change
                newPosition = get(gca, 'position');
                set(gca, 'position', [newPosition(1) newPosition(2) ...
                    originalPosition(3) originalPosition(4)]);
            end
        end
    else

        % Animation of color 
        zoomInFactorToBeginAlpha = 0.65;
        zoomInFactorMin = 0.3;
        if (zoomInFactor <= zoomInFactorToBeginAlpha) 
            alpha = (zoomInFactorToBeginAlpha - zoomInFactor) / ...
                (zoomInFactorToBeginAlpha - zoomInFactorMin);
        else
            alpha = 0.0;
        end
        lineWidth = 0.1 + alpha;
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 1, ...
               'colsNum', 1, ...
               'heightMargin', 0.0, ...
               'widthMargin', 0.0, ...
               'leftMargin', 0.009, ...
               'rightMargin', 0.06, ...
               'bottomMargin', 0.04, ...
               'topMargin', 0.03);
        
        subplot('Position', subplotPosVectors(1, 1).v);
        set(gca, 'Color', [0 0 0]);
        showXticks = true;
        showYticks = false;

        if (isempty(xRange))
            xRange = [sampledHexMosaicXaxis(1) - obj.pigment.width ...
                sampledHexMosaicXaxis(end) + obj.pigment.width];
        else
            xRange(1) = xRange(1) - obj.pigment.width;
            xRange(2) = xRange(2) + obj.pigment.width;
        end
        if (isempty(yRange))
            yRange = [sampledHexMosaicYaxis(1) - obj.pigment.width, ...
                sampledHexMosaicYaxis(end) + obj.pigment.width];
        else
            yRange(1) = yRange(1) - obj.pigment.width;
            yRange(2) = yRange(2) + obj.pigment.width;
        end

        xRange = xRange * zoomInFactor;
        yRange = yRange * zoomInFactor;
        set(gca, 'XLim', xRange, 'YLim', yRange, 'CLim', [0 1]);

        cMapLevels = size(cMap, 1);

        for coneType = 2:4
            idx = find(obj.pattern == coneType);
            [iRows, iCols] = ind2sub(size(obj.pattern), idx);
            coneXcoords = sampledHexMosaicXaxis(iCols);
            coneYcoords = sampledHexMosaicYaxis(iRows);
            activationsNlevels = round(...
                (activation(idx) - activationRange(1)) / ...
                (activationRange(2) - activationRange(1)) * cMapLevels);
            faceColorsNormalizedValues = activationsNlevels / cMapLevels;
            faceColorsNormalizedValues(faceColorsNormalizedValues < 0) = 0;
            faceColorsNormalizedValues(faceColorsNormalizedValues > 1) = 1;
           
            coneXcoords = coneXcoords(:);
            coneYcoords = coneYcoords(:);
            idx = find( (coneXcoords > xRange(1)) & ...
                (coneXcoords < xRange(2)) & (coneYcoords > yRange(1)) & ...
                (coneYcoords < yRange(2)));
            coneXcoords = coneXcoords(idx);
            coneYcoords = coneYcoords(idx);
            faceColorsNormalizedValues = faceColorsNormalizedValues(idx);

            if (coneType == 2)
                edgeColor = [1.0 0.1 0.1] * alpha + ...
                    (1 - alpha) * [0.2 0.2 0.2];
            elseif (coneType == 3)
                edgeColor = [0.1 1.0 0.1] * alpha + ...
                    (1 - alpha) * [0.2 0.2 0.2];
            else
                edgeColor = [0.1 0.6 1.0] * alpha + ...
                    (1 - alpha) * [0.2 0.2 0.2];
            end

            if (lineWidth < 0.2), edgeColor = 'none'; end

            renderPatchArray(apertureOutline, coneXcoords, coneYcoords, ...
                faceColorsNormalizedValues, edgeColor, lineWidth);
            if (coneType == 2), hold on; end
        end
        hold off;

        xTicks = obj.center(1) + 1e-6 * (-75:75:75);
        yTicks = obj.center(2) + 1e-6 * (-75:75:75);
        xTickLabels = sprintf('%2.0f um\n', xTicks * 1e6);
        yTickLabels = sprintf('%2.0f um\n', yTicks * 1e6);

        axis 'equal';
        set(gca, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8]);
        set(gca, 'XTick', xTicks, 'YTick', yTicks, ...
            'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels, ...
            'FontSize', 18, 'FontName', 'Menlo');
        if (~showXticks), set(gca, 'XTick', []); end
        if (~showYticks), set(gca, 'YTick', []); end
        box on;
        grid off;

        % Add colorbar
        ticks = 0:0.1:1.0;
        tickLabels = sprintf('%2.0f\n', round(activationRange(1) + ...
            ticks * (activationRange(2) - activationRange(1))));
        originalPosition = get(gca, 'position');
        hCbar = colorbar('eastoutside', 'peer', gca, 'Ticks', ticks, ...
            'TickLabels', tickLabels);
        hCbar.Orientation = 'vertical';
        hCbar.Label.String = signalName;
        hCbar.FontSize = 16;
        hCbar.FontName = 'Menlo';
        hCbar.Color = [0.9 0.9 0.5];
        % The addition changes the figure size, so undo this change
        newPosition = get(gca, 'position');
        set(gca, 'position', [newPosition(1), newPosition(2), ...
            originalPosition(3), originalPosition(4)]);
        
        if (~isnan(instanceIndex))
            if (~isempty(activationTime))
                title(sprintf(['t = %7.1f msec              (response ' ...
                    'instance: #%3d)'], activationTime * 1000, ...
                    instanceIndex), 'FontSize', 18, 'Color', [1 1 1], ...
                    'FontName', 'Menlo');
            else
                title(sprintf('instance: %d', instanceIndex), ...
                    'FontSize', 18, 'Color', [1 1 1], 'FontName', 'Menlo');
            end
        else
            if (~isempty(activationTime))
                title(sprintf('t = %02.1f msec', activationTime * 1000),...
                    'FontSize', 18, 'Color', [1 1 1], 'FontName', 'Menlo');
            end
        end
    end
    colormap(cMap);
    drawnow
end

function renderPatchArray(pixelOutline, xCoords, yCoords, ...
    faceColorsNormalizedValues, edgeColor, lineWidth)
% Render the patch array
%
% Syntax:
%   renderPatchArray(pixelOutline, xCoords, yCoords, ...
%       faceColorsNormalizedValues, edgeColor, lineWidth)
%
% Description:
%    Render the patch array using the provided variables
%
% Inputs:
%    pixelOutline               - The outline
%    xCoords                    - X-axis coordinates
%    yCoords                    - Y-axis coordinates
%    faceColorsNormalizedValues - Normalized values for the face colors
%    edgeColor                  - Color of the edges
%    lineWidth                  - Width of the boundary line
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    verticesPerCone = numel(pixelOutline.x);
    verticesList = zeros(verticesPerCone * numel(xCoords), 2);
    facesList = [];
    colors = [];
    for coneIndex = 1:numel(xCoords)
        idx = (coneIndex-1)*verticesPerCone + (1:verticesPerCone);
        verticesList(idx, 1) = pixelOutline.x(:) + xCoords(coneIndex);
        verticesList(idx, 2) = pixelOutline.y(:) + yCoords(coneIndex);
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
    patch(S);
end

function hFig = visualizeMosaicActivationsAsDensityMaps(...
    obj, activation, cMap, signalName, signalRange, figureSize)
% Visualize mosaic activations as density maps
%
% Syntax:
%   hFig = visualizeMosaicActivationsAsDensityMaps(obj, activation, ...
%       cMap, signalName, signalRange, figureSize)
%
% Description:
%    Visualize the mosaic activations as density maps
%
% Inputs:
%    obj                - The cone mosaic hex object
%    activation         - The activation you wish to map
%    cMap               - The color map
%    signalName         - The name of the signal
%    signalRange        - The range associated with the signal
%    figureSize         - The size of the figure
%
% Outputs:
%    hFig               - The figure handle
%
% Optional key/value pairs:
%    None.
%

    % Compute activation image maps
    [activationImage, activationImageLMScone, sampledHexMosaicXaxis, ...
        sampledHexMosaicYaxis] = obj.computeActivationDensityMap(...
        activation);

    activeConesActivations = activation(obj.pattern>1);
    % activationRange = prctile(activeConesActivations, [10 90]);
    if (isempty(signalRange))
        activationRange = [min(activeConesActivations(:)), ...
            max(activeConesActivations(:))];
    else
         activationRange = signalRange;
    end

    hFig = figure();
    set(hFig, 'Position', cat(2, [10 10], figureSize), 'Color', [0 0 0]);
    % , 'MenuBar', 'None');

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 3, ...
           'heightMargin', 0.01, ...
           'widthMargin', 0.01, ...
           'leftMargin', 0.01, ...
           'rightMargin', 0.005, ...
           'bottomMargin', 0.005, ...
           'topMargin', 0.005);

    for subplotIndex = 1:4
        if (subplotIndex < 4)
            subplot('Position', subplotPosVectors(1, subplotIndex).v);
        else
            subplot('Position', subplotPosVectors(2, 2).v);
        end
        set(gca, 'Color', [0 0 0]);

        switch subplotIndex
            case 1
                activationMapImage = ...
                    squeeze(activationImageLMScone(:, :, 1));
                subplotTitle = 'L-cone submosaic activation';
                showXticks = false;
                showYticks = false;
            case 2
                activationMapImage = ...
                    squeeze(activationImageLMScone(:, :, 2));
                subplotTitle = 'M-cone submosaic activation';
                showXticks = false;
                showYticks = false;
            case 3
                activationMapImage = ...
                    squeeze(activationImageLMScone(:, :, 3));
                subplotTitle = 'S-cone submosaic activation';
                showXticks = false;
                showYticks = false;
            case 4
                activationMapImage = activationImage;
                subplotTitle = 'total mosaic activation';
                showXticks = true;
                showYticks = true;
        end

        imagesc(sampledHexMosaicXaxis, sampledHexMosaicYaxis, ...
            activationMapImage);
        axis 'image';
        axis 'xy';
        xTicks = [sampledHexMosaicXaxis(1), obj.center(1), ...
            sampledHexMosaicXaxis(end)];
        yTicks = [sampledHexMosaicYaxis(1), obj.center(2), ...
            sampledHexMosaicYaxis(end)];
        xTickLabels = sprintf('%2.0f um\n', xTicks * 1e6);
        yTickLabels = sprintf('%2.0f um\n', yTicks * 1e6);
        set(gca, 'XTick', xTicks, 'YTick', yTicks, ...
            'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels, ...
            'XColor', [0.0 0.0 0.0], 'YColor', [0.0 0.0 0.0]);
        set(gca, 'CLim', activationRange);
        if (~showXticks), set(gca, 'XTick', []); end
        if (~showYticks), set(gca, 'YTick', []); end
        box on;
        grid off;
        set(gca, 'XLim', [sampledHexMosaicXaxis(1) - obj.pigment.width ...
            sampledHexMosaicXaxis(end) + obj.pigment.width]);
        set(gca, 'YLim', [sampledHexMosaicYaxis(1) - obj.pigment.width ...
        sampledHexMosaicYaxis(end) + obj.pigment.width]);
        set(gca, 'FontSize', 18, 'FontName', 'Menlo');
        title(subplotTitle, 'FontSize', 18, 'Color', [1 1 1], ...
            'FontName', 'Menlo');

        if (subplotIndex == 4)
            % Add colorbar
            originalPosition = get(gca, 'position');
            hCbar = colorbar('eastoutside', 'peer', gca);
            % , ... 'Ticks', cbarStruct.ticks, 'TickLabels', ...
            %  cbarStruct.tickLabels);
            hCbar.Orientation = 'vertical';
            hCbar.Label.String = signalName;
            hCbar.FontSize = 16;
            hCbar.FontName = 'Menlo';
            hCbar.Color = [0.9 0.9 0.5];
            % The addition changes the figure size, so undo this change
            newPosition = get(gca, 'position');
            set(gca, 'position', [newPosition(1), newPosition(2), ...
                originalPosition(3), originalPosition(4)]);
        end

    end

    cMap(1, :) = 0;
    colormap(cMap);
    drawnow
end
