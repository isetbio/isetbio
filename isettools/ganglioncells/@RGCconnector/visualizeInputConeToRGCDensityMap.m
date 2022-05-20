function visualizeInputConeToRGCDensityMap(obj)

    % Compute the spatial support for the map
    [X,Y] = meshgrid(obj.coneToRGCDensityRatioComputeStruct.xSupport, obj.coneToRGCDensityRatioComputeStruct.ySupport);
    
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

    hFig = figure(999); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1024 660]);
    ax = subplot(2,2,1);
    [~,~,XLims, YLims] = obj.visualizeInputMosaics(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'thetaSamples', 30, ...
        'titleString', 'starting positions');
    set(ax, 'FontSize', 16)

    

    ax = subplot(2,2,2);
    imagesc(ax,obj.coneToRGCDensityRatioComputeStruct.xSupport, obj.coneToRGCDensityRatioComputeStruct.ySupport, ratioMap);
    axis(ax, 'equal'); axis(ax, 'xy');
    set(ax, 'CLim', ratioLimits, 'XLim', XLims, 'YLim', YLims);
    set(ax, 'FontSize', 16);
    title(ax,'cone/RGC density ratio map');
    colormap(ax,brewermap(1024, '*YlGnBu'));
    colorbar(ax);
    

    ax = subplot(2,2,3);
    imagesc(ax,obj.coneToRGCDensityRatioComputeStruct.xSupport, obj.coneToRGCDensityRatioComputeStruct.ySupport, coneDensityMap);
    axis(ax, 'equal'); axis(ax, 'xy');
    set(ax, 'CLim', coneDensityLimits, 'XLim', XLims, 'YLim', YLims);
    title(ax,'cone density map');
    set(ax, 'FontSize', 16);
    colormap(ax,brewermap(1024, '*YlGnBu'));
    colorbar(ax);

    ax = subplot(2,2,4);
    imagesc(ax,obj.coneToRGCDensityRatioComputeStruct.xSupport, obj.coneToRGCDensityRatioComputeStruct.ySupport, rgcDensityMap);
    axis(ax, 'equal'); axis(ax, 'xy');
    set(ax, 'CLim', rgcDensityLimits, 'XLim', XLims, 'YLim', YLims);
    title(ax,'rgc density map');
    set(ax, 'FontSize', 16);
    colormap(ax,brewermap(1024, '*YlGnBu'));
    colorbar(ax);
end