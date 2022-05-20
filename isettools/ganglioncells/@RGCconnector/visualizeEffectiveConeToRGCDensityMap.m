function visualizeEffectiveConeToRGCDensityMap(obj)

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

    hFig = figure(998); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1650 700]);
    ax = subplot(1,2,1);
    [~,~,XLims, YLims] = obj.visualizeInputMosaics(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'thetaSamples', 30, ...
        'titleString', 'starting positions');
    set(ax, 'FontSize', 16)

    ax = subplot(1,2,2);
    for i = 1:numel(X)
       text(ax, X(i), Y(i), sprintf('%2.1f', ratioMap(i)))
    end
    set(ax, 'FontSize', 16);
    axis(ax, 'equal'); axis(ax, 'xy');
    set(ax, 'XLim', XLims, 'YLim', YLims);

end