function hFig = visualizeCurrentConnectivity(obj, figNo, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('titleString', '', @ischar);
    p.addParameter('visualizedDestinationRFindices', [], @isnumeric);
    p.addParameter('onlyVisualizeConnectivity', false, @islogical);
    p.parse(varargin{:});
    titleString = p.Results.titleString;
    visualizedDestinationRFindices = p.Results.visualizedDestinationRFindices;
    onlyVisualizeConnectivity = p.Results.onlyVisualizeConnectivity;


    hFig = figure(figNo); clf;
    if (onlyVisualizeConnectivity)
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 1, ...
           'heightMargin',  0.07, ...
           'widthMargin',    0.04, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.03);
    else
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 2, ...
           'heightMargin',  0.07, ...
           'widthMargin',    0.04, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.03);
    end
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1650 990]);

    if (onlyVisualizeConnectivity)
        minXYdest = min(obj.destinationLattice.RFpositionsMicrons,[],1);
        maxXYdest = max(obj.destinationLattice.RFpositionsMicrons,[],1);
        maxDestSpacing = max(obj.destinationLattice.RFspacingsMicrons);
        XLims = [minXYdest(1) maxXYdest(1)] + maxDestSpacing*[-1 1];
        YLims = [minXYdest(2) maxXYdest(2)] + maxDestSpacing*[-1 1];
    else
        % The input lattices on the rop-right plot
        ax = subplot('Position', subplotPosVectors(1,1).v);
        [~,~,XLims, YLims] = obj.visualizeInputLattices(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'thetaSamples', 30);
        set(ax, 'FontSize', 16);
    end

    if (isempty(titleString))
        titleString = 'lattice wiring (stars: destination RFs with 0 inputs)';
    end

    % The connectivity (pooling of destination lattice RFs) in the bottom-right plot
    if (~onlyVisualizeConnectivity)
        ax = subplot('Position', subplotPosVectors(1,2).v);
    else
        ax = subplot('Position', subplotPosVectors(1,1).v);
    end
    visualizeDestinationLatticePooling(obj, ...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'titleString', titleString, ...
        'titleWithPoolingStats', true, ...
        'visualizedDestinationRFindices', visualizedDestinationRFindices, ...
        'thetaSamples', 10, ...
        'XLims', XLims, ...
        'YLims', YLims);

    % Report statistics of spatial/chromatic cost
    if (~onlyVisualizeConnectivity)
        if (~isempty(obj.connectivityMatrix))
            if (~isempty(obj.destinationRFspacingsFromCentroids))

                % Compute costs across entire destinationRF mosaic
                theTotalPoolingCosts = obj.totalPoolingCosts();

                ax1 = subplot('Position', subplotPosVectors(2,1).v);
                ax2 = subplot('Position', subplotPosVectors(2,2).v);
                obj.visualizeCostComponentStatistics(ax1, ax2, theTotalPoolingCosts);
            end
        end
    end
    drawnow;
end

function visualizeDestinationLatticePooling(obj, varargin)
    
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('titleString', '', @(x)(isempty(x) || (ischar(x))));
    p.addParameter('titleWithPoolingStats', false, @islogical);
    p.addParameter('thetaSamples', 20, @isnumeric);
    p.addParameter('connectPooledSourceRFs', true, @islogical);
    p.addParameter('visualizedDestinationRFindices', [], @isnumeric);
    p.addParameter('displayDestinationRFID', false, @islogical);
    p.addParameter('XLims', [], @isnumeric);
    p.addParameter('YLims', [], @isnumeric);
    p.parse(varargin{:});
    
    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    titleString = p.Results.titleString;
    titleWithPoolingStats = p.Results.titleWithPoolingStats;
    thetaSamples = p.Results.thetaSamples;
    connectPooledSourceRFs = p.Results.connectPooledSourceRFs;
    displayDestinationRFID = p.Results.displayDestinationRFID;
    visualizedDestinationRFindices = p.Results.visualizedDestinationRFindices;

    XLims = p.Results.XLims;
    YLims = p.Results.YLims;

    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1], 'Position', [10 10 850 800]);
        end
        ax = subplot('Position', [0.05 0.07 0.93 0.9]);
    end
    
    hold(ax, 'off')

    % Small source RF outline so we can see the pooling better
    thetas = linspace(0,360,4);
    sourceRFoutline = 0.3*[cosd(thetas); sind(thetas)]';

    
    % Visualize the pooling of source lattice RFs by the destination lattice RFs
    cMap = [0.7 0.85 0.99; 0 0 0];

    if (isempty(visualizedDestinationRFindices))
        visualizedDestinationRFindices = 1:size(obj.connectivityMatrix,2);
    end
    inputsPerDestinationRF = nan(1, numel(visualizedDestinationRFindices));


    for visualizedDestinationRFindex = 1:numel(visualizedDestinationRFindices)

        % Retrieve the visualized destinationRF index
        destinationRFindex = visualizedDestinationRFindices(visualizedDestinationRFindex);

        % Retrieve connection data
        % Indices of source lattice RFs pooled by this destination RF
        indicesOfConnectedSourceRFs = find(squeeze(obj.connectivityMatrix(:, destinationRFindex))>0);
        
        if (isempty(indicesOfConnectedSourceRFs))
            % No sourceRFs pooled by this destination RF. Identify with a black star
            plot(ax, obj.destinationLattice.RFpositionsMicrons(destinationRFindex,1), ...
                     obj.destinationLattice.RFpositionsMicrons(destinationRFindex,2), ...
                     'h', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 16, 'LineWidth', 1);
            continue;
        end

        % Update inputsPerDestinationRF
        inputsPerDestinationRF(visualizedDestinationRFindex) = numel(indicesOfConnectedSourceRFs);

        % Positions and spacings of these pooled source lattice RFs
        positionsOfConnectedSourceRFs = obj.sourceLattice.RFpositionsMicrons(indicesOfConnectedSourceRFs,:);
        spacingsOfConnectedSourceRFs = obj.sourceLattice.RFspacingsMicrons(indicesOfConnectedSourceRFs);

        % Weights of these pooled source lattice RFs
        weightsOfConnectedSourceRFs  = full(obj.connectivityMatrix(indicesOfConnectedSourceRFs, destinationRFindex));

        % Generate RF visualization struct
        rfVisualizationStruct = destinationRFoutlineFromConnectedSourceRFs(...
                positionsOfConnectedSourceRFs, spacingsOfConnectedSourceRFs, weightsOfConnectedSourceRFs);

        % Plot contour of input source RFs
        zLevels(1) = 0.05*min(weightsOfConnectedSourceRFs);
        zLevels(2) = max(weightsOfConnectedSourceRFs);
        faceAlpha = 0.7;
        lineStyle = '-';
        lineWidth = 1.0;
        renderTransparentContourPlot(ax, ...
            rfVisualizationStruct.spatialSupportXY, ...
            rfVisualizationStruct.rfMap, ...
            zLevels, faceAlpha, cMap, lineStyle, lineWidth);
    end % visualizedDestinationRFindex

    if (connectPooledSourceRFs)
        for visualizedDestinationRFindex = 1:numel(visualizedDestinationRFindices)

            % Retrieve the visualized destinationRF index
            destinationRFindex = visualizedDestinationRFindices(visualizedDestinationRFindex);

            % Retrieve connection data
            % Indices of source lattice RFs pooled by this destination RF
            indicesOfConnectedSourceRFs = find(squeeze(obj.connectivityMatrix(:, destinationRFindex))>0);

            % Positions and weights of pooled source lattice RFs
            positionsOfConnectedSourceRFs = obj.sourceLattice.RFpositionsMicrons(indicesOfConnectedSourceRFs,:);
            weightsOfConnectedSourceRFs = full(obj.connectivityMatrix(indicesOfConnectedSourceRFs, destinationRFindex));

            centroid = obj.destinationRFcentroidsFromInputs(destinationRFindex,:);
            if (displayDestinationRFID)
               text(ax, centroid(1), centroid(2), sprintf('%d', destinationRFindex), 'Color', [0 1 0], 'FontSize', 12, 'BackgroundColor', [0 0 0]);
            end

            if (numel(indicesOfConnectedSourceRFs)>1)
                for inputIndex = 1:numel(indicesOfConnectedSourceRFs)
                    xx = [centroid(1) positionsOfConnectedSourceRFs(inputIndex ,1)];
                    yy = [centroid(2) positionsOfConnectedSourceRFs(inputIndex ,2)];
                    plot(ax, xx, yy, 'k-', 'LineWidth', weightsOfConnectedSourceRFs(inputIndex)/max(weightsOfConnectedSourceRFs)*2.0);
                end
            elseif (numel(indicesOfConnectedSourceRFs) == 1)
                plot(ax, centroid(1), centroid(2), 'k.');
            end
        end % visualizedDestinationRFindex
    end

    % Plot the source lattice RFs
    obj.visualizeSourceLatticeRFs(ax, sourceRFoutline);

    % Finalize
    axis(ax, 'equal');

    if (~isempty(XLims)) && (~isempty(YLims))
        set(ax, 'XLim', XLims, 'YLim', YLims);
    end
    set(ax, 'FontSize', 16);

    if (titleWithPoolingStats)
        poolingStatsString = sprintf('inputs/destination RF: min:%d, max:%d, mean:%2.3f', ...
            min(inputsPerDestinationRF, [], 'omitnan'), ...
            max(inputsPerDestinationRF, [], 'omitnan'), ...
            mean(inputsPerDestinationRF, 'omitnan'));
        titleString = sprintf('%s\n%s', titleString, poolingStatsString);
    end

    if (~isempty(titleString))
        title(ax, titleString);
    end

    box(ax, 'on')
    xlabel(ax, 'microns'); 
    ylabel(ax, 'microns');
    
    drawnow;
end

function rfVisualizationStruct = destinationRFoutlineFromConnectedSourceRFs(...
            inputPositions, inputSpacings, inputWeights, varargin)
    p = inputParser;
    p.addOptional('xSupport', [], @(x)(isempty(x) || (isnumeric(x))));
    p.addOptional('ySupport', [], @(x)(isempty(x) || (isnumeric(x))));
    p.addOptional('spatialSupportSamples', 60, @(x)(isempty(x) || (isnumeric(x))));
    p.parse(varargin{:});
    xSupport = p.Results.xSupport;
    ySupport = p.Results.ySupport;
    spatialSupportSamples = p.Results.spatialSupportSamples;

    % Compute spatial support
    xSep = max(inputSpacings);

    if (isempty(xSupport))
        xx = inputPositions(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = inputPositions(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end
    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);
    RF = zeros(size(X));

    % Output
    theInputLineSpreadFunctionsXY = cell(2, numel(inputWeights));

    for inputIndex = 1:numel(inputWeights)
        % Characteristic radius of the input RF
        rC = 0.204*sqrt(2.0)*inputSpacings(inputIndex);
        % Compute aperture2D x weight
        XX = X-inputPositions(inputIndex,1);
        YY = Y-inputPositions(inputIndex,2);
        theAperture2D = inputWeights(inputIndex) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        % Accumulate 2D apertures
        RF = RF + theAperture2D;
        % 1D Line spread functions (XY) for each cone aperture
        theInputLineSpreadFunctionsXY{1,inputIndex} = sum(theAperture2D,1);
        theInputLineSpreadFunctionsXY{2,inputIndex} = sum(theAperture2D,2);
    end
    
    % Assemble rfVisualizationStruct
    rfVisualizationStruct.spatialSupportXY = spatialSupportXY;
    rfVisualizationStruct.rfMap = RF;
    rfVisualizationStruct.theInputLineSpreadFunctionsXY = theInputLineSpreadFunctionsXY;
end

function renderTransparentContourPlot(axesHandle, spatialSupportXY, zData, ...
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