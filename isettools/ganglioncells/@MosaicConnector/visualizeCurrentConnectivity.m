function hFig = visualizeCurrentConnectivity(obj, figNo, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('titleString', '', @ischar);
    p.parse(varargin{:});
    titleString = p.Results.titleString;


    hFig = figure(figNo); clf;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 2, ...
           'heightMargin',  0.07, ...
           'widthMargin',    0.04, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.03);
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1650 990]);

    % The input lattices on the rop-right plot
    ax = subplot('Position', subplotPosVectors(1,1).v);
    [~,~,XLims, YLims] = obj.visualizeInputLattices(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'thetaSamples', 30);
    set(ax, 'FontSize', 16);

    if (isempty(titleString))
        titleString = 'lattice wiring (stars: destination RFs with 0 inputs)';
    end

    % The connectivity (pooling of destination lattice RFs) in the bottom-right plot
    ax = subplot('Position', subplotPosVectors(1,2).v);
    obj.visualizeDestinationLatticePooling(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'titleString', titleString, ...
        'titleWithPoolingStats', true, ...
        'XLims', XLims, ...
        'YLims', YLims);

    % Report statistics of spatial/chromatic cost
    if (~isempty(obj.connectivityMatrix))
        if (~isempty(obj.destinationRFspacingsFromCentroids))

            % Compute costs across entire destinationRF mosaic
            theCostComponentsMatrix = obj.totalInputMaintenanceCost();

            ax1 = subplot('Position', subplotPosVectors(2,1).v);
            ax2 = subplot('Position', subplotPosVectors(2,2).v);
            obj.visualizeCostComponentStatistics(ax1, ax2, theCostComponentsMatrix);
        end
    end

    drawnow;


end

