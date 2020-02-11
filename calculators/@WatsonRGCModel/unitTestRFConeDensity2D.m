function unitTestRFConeDensity2D()

    eccMinDegs = 0.15;
    eccMaxDegs = 10;
    eccSamplesNum = 50;
    eccDegsInREVisualSpace = logspace(log10(eccMinDegs), log10(eccMaxDegs), eccSamplesNum);
    eccDegsInREVisualSpace = [-fliplr(eccDegsInREVisualSpace) 0 eccDegsInREVisualSpace];
    eccUnits = 'visual deg';
    densityUnits = 'visual deg^2';
    
    obj = WatsonRGCModel();
    [Z, spatialSupport] = obj.compute2DConeRFDensity(eccDegsInREVisualSpace, eccUnits, densityUnits);
    
    figure(2); clf;
    X = squeeze(spatialSupport(1,:));
    Y = squeeze(spatialSupport(2,:));
    zLevels = logspace(log10(300), log10(10000), 10); %  [300 500 750 1000 1500 2000 3000 5000 10000];
    contourf(X,Y, log10(Z), log10(zLevels));
    cMap = brewermap(1024, 'greys');
    axis 'square';
    colormap(cMap);
    xLims = eccMaxDegs*[-1 1];
    yLims = xLims;
    xTicks = -100:2:100;
    yTicks = xTicks;
    set(gca, 'XLim', xLims, 'YLim', yLims, ...
        'XTick', xTicks, ...
        'YTick', yTicks, ...
        'FontSize', obj.figurePrefs.fontSize);
    colorbar
    % Finalize figure
    grid(gca, obj.figurePrefs.grid);
end
