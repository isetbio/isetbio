function visualizeEffectiveConeToRGCDensityMap(obj, figNo)

    % Compute the spatial support for the map
    X = obj.RGCRFpositionsMicrons(:,1);
    Y = obj.RGCRFpositionsMicrons(:,2);

    % Compute the density map for the cones
    coneDensityMap = obj.coneToRGCDensityRatioComputeStruct.coneDensityFunctionHandle(X,Y);

    % Compute the density map for the RGCs
    rgcDensityMap = obj.coneToRGCDensityRatioComputeStruct.RGCDensityFunctionHandle(X,Y);

    % Compute the ratio of densities map
    ratioMap = coneDensityMap ./ rgcDensityMap;

    % Limits
    coneDensityLimits = prctile([coneDensityMap(:);rgcDensityMap(:)], [1 99]);
    rgcDensityLimits = coneDensityLimits;
    ratioLimits = prctile(ratioMap(:), [5 95]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 1, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.0);

    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1650 1000]);

    ax = subplot('Position', subplotPosVectors(1,1).v);
    [~,~,XLims, YLims] = obj.visualizeInputMosaics(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'thetaSamples', 30, ...
        'titleString', 'starting positions');
   

    ax = subplot('Position', subplotPosVectors(2,1).v);
    for i = 1:numel(X)
       text(ax, X(i), Y(i), sprintf('%2.1f', ratioMap(i)))
    end
    set(ax, 'FontSize', 16);
    axis(ax, 'equal'); axis(ax, 'xy');
    title(ax, 'theoretical cone-to-RGC ratio');
    set(ax, 'XLim', XLims, 'YLim', YLims);
    box(ax, 'on')
    drawnow
end