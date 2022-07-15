function visualizeLocalConeToRGCDensityMaps(localConeToRGCDensityRatioStruct, RGCRFposMicrons, theInputConeMosaic)

    % Visualize 
    [X,Y] = meshgrid(localConeToRGCDensityRatioStruct.xSupport, localConeToRGCDensityRatioStruct.ySupport);
    coneDensityMap = localConeToRGCDensityRatioStruct.coneDensityFunctionHandle(X,Y);
    rgcDensityMap = localConeToRGCDensityRatioStruct.RGCDensityFunctionHandle(X,Y);
    ratioMap = coneDensityMap ./ rgcDensityMap;

    coneDensityLimits = prctile(coneDensityMap(:), [1 99]);
    rgcDensityLimits = prctile(rgcDensityMap(:), [1 99]);
    ratioLimits = prctile(ratioMap(:), [5 95]);


    hFig = figure(999); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1024 660]);
    ax = subplot(2,2,1);
    RGCRFconnector.plotRGCRFpos(RGCRFposMicrons ,...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'inputConeMosaic', theInputConeMosaic, ...
        'thetaSamples', 30, ...
        'titleString', 'starting positions');
    set(ax, 'FontSize', 16)

    

    ax = subplot(2,2,2);
    imagesc(ax,localConeToRGCDensityRatioStruct.xSupport, localConeToRGCDensityRatioStruct.ySupport, ratioMap);
    set(ax, 'CLim', ratioLimits)
    axis 'equal'
    axis 'xy'
    set(ax, 'XLim', [localConeToRGCDensityRatioStruct.xSupport(1) localConeToRGCDensityRatioStruct.xSupport(end)], ...
            'YLim', [localConeToRGCDensityRatioStruct.ySupport(1) localConeToRGCDensityRatioStruct.ySupport(end)]);
    title(ax,'cone/RGC density ratio map');
    set(ax, 'FontSize', 16)
    colormap(ax,brewermap(1024, '*YlGnBu'));
    colorbar(ax);
    

    ax = subplot(2,2,3);
    imagesc(ax,localConeToRGCDensityRatioStruct.xSupport, localConeToRGCDensityRatioStruct.ySupport, coneDensityMap);
    set(ax, 'CLim', coneDensityLimits);
    axis 'equal'
    axis 'xy'
    set(ax, 'XLim', [localConeToRGCDensityRatioStruct.xSupport(1) localConeToRGCDensityRatioStruct.xSupport(end)], ...
            'YLim', [localConeToRGCDensityRatioStruct.ySupport(1) localConeToRGCDensityRatioStruct.ySupport(end)]);
    title(ax,'cone density map');
    set(ax, 'FontSize', 16)
    colormap(ax,brewermap(1024, '*YlGnBu'));
    colorbar(ax);

    ax = subplot(2,2,4);
    imagesc(ax,localConeToRGCDensityRatioStruct.xSupport, localConeToRGCDensityRatioStruct.ySupport, rgcDensityMap);
    set(ax, 'CLim', rgcDensityLimits);
    axis 'equal'
    axis 'xy'
    set(ax, 'XLim', [localConeToRGCDensityRatioStruct.xSupport(1) localConeToRGCDensityRatioStruct.xSupport(end)], ...
            'YLim', [localConeToRGCDensityRatioStruct.ySupport(1) localConeToRGCDensityRatioStruct.ySupport(end)]);
    title(ax,'rgc density map');
    set(ax, 'FontSize', 16)
    colormap(ax,brewermap(1024, '*YlGnBu'));
    colorbar(ax);
end
