function rfDensity2DMaps(figNo, neuronType, varargin)

    validDensityUnits = {'deg^2', 'mm^2'};
    validSpatialSupportUnits = {'degs', 'mm', 'microns'};
    
    % Configure inputs parser
    p = inputParser;
    p.addRequired('figNo', @isnumeric);
    p.addRequired('neuronType', @ischar);
    p.addParameter('spatialSupport', @isstruct);
    p.addParameter('spatialSupportUnits', @(x)(ischar(x) && (ismember(x, validSpatialSupportUnits))));
    p.addParameter('densityRange', @isnumeric);
    p.addParameter('densityUnits', @(x)(ischar(x) && (ismember(x, validDensityUnits))));
    
 
    % Parse input
    p.parse(figNo, neuronType, varargin{:});
    spatialSupport = p.Results.spatialSupport;
    spatialSupportUnits = p.Results.spatialSupportUnits;
    densityRange = p.Results.densityRange;
    densityUnits = p.Results.densityUnits;
    
    % Spatial sampling grid
    switch (spatialSupportUnits)
        case 'microns'
            ; % do nothing
        case 'degs'
            % Transform from degs to microns
            spatialSupport.maxEcc = RGCmodels.Watson.convert.rhoDegsToMMs(spatialSupport.maxEcc)*1e3;
        case 'mm'
            % Transform from mm to microns
            spatialSupport.maxEcc = spatialSupport.maxEcc * 1e3;
    end
    
    
    
    rfPosXmicrons = linspace(0, spatialSupport.maxEcc, spatialSupport.eccSamples);
    rfPosXmicrons = cat(2,-fliplr(rfPosXmicrons(2:end)), rfPosXmicrons);
    rfPosYmicrons = rfPosXmicrons;
    [X,Y] = meshgrid(rfPosXmicrons, rfPosYmicrons);
    rfPosMicrons = zeros(numel(X),2);
    
    rfPosMicrons(:,1) = X(:);
    rfPosMicrons(:,2) = Y(:);
    
    % Compute density at spatial sampling grid
    % Left eye
    whichEye = RGCmodels.Watson.constants.leftEye;
    [~, ~, leftRetinaDensityMapDegs2, leftRetinaDensityMapMMs2] = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(whichEye, rfPosMicrons, neuronType);
    switch (densityUnits)
        case 'deg^2'
             leftRetinaDensityMap = reshape(leftRetinaDensityMapDegs2, size(X));
        case 'mm^2'
             leftRetinaDensityMap = reshape(leftRetinaDensityMapMMs2, size(X));
    end

    
    % Right eye
    whichEye = RGCmodels.Watson.constants.rightEye;
    [~, ~, rightRetinaDensityMapDegs2, rightRetinaDensityMapMMs2] = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(whichEye, rfPosMicrons, neuronType);
    switch (densityUnits)
        case 'deg^2'
             rightRetinaDensityMap = reshape(rightRetinaDensityMapDegs2, size(X));
        case 'mm^2'
             rightRetinaDensityMap = reshape(rightRetinaDensityMapMMs2, size(X));
    end

    % Transform rfPosXmicrons to desired units
    switch (spatialSupportUnits)
        case 'microns'
            ; % do nothing
        case 'degs'
            % Transform from microns to degs
            rfPosXmicrons = RGCmodels.Watson.convert.rhoMMsToDegs(rfPosXmicrons/1e3);
            rfPosYmicrons = RGCmodels.Watson.convert.rhoMMsToDegs(rfPosYmicrons/1e3);
        case 'mm'
            % Transform from microns to mm
            rfPosXmicrons = rfPosXmicrons / 1e3;
            rfPosYmicrons = rfPosYmicrons / 1e3;
    end
    
    % Plot
    hFig = figure(figNo);
    set(hFig, 'Position', [10 10 1500 800], 'Color', [1 1 1], 'Name', neuronType);
    
    % Left eye
    ax = subplot(1,2,1);
    plotData(ax, rfPosXmicrons, rfPosYmicrons, leftRetinaDensityMap, densityRange, 'left retina', densityUnits, spatialSupportUnits);
    
    % Right eye
    ax = subplot(1,2,2);
    plotData(ax, rfPosXmicrons, rfPosYmicrons, rightRetinaDensityMap, densityRange, 'right retina', densityUnits, spatialSupportUnits);
end

function plotData(ax, rfPosXmicrons, rfPosYmicrons, densityMap, densityRange, titleLabel, densityUnits, spatialSupportUnits)
    contourf(ax, rfPosXmicrons, rfPosYmicrons, log10(densityMap), linspace(log10(densityRange(1)), log10(densityRange(2)),20));
    hold(ax, 'on');
    m1 = min(densityMap(:));
    [m2, idx] = max(densityMap(:));
    [peakRow, peakCol] = ind2sub(size(densityMap), idx);
    xMin = min(rfPosXmicrons);
    xMax = max(rfPosXmicrons);
    yMin = min(rfPosYmicrons);
    yMax = max(rfPosYmicrons);
    plot(ax,rfPosXmicrons, yMin+10*(yMax-yMin)*(densityMap(peakRow,:)-m1)/(m2-m1), 'k-', 'LineWidth', 1.5);
    axis(ax,'equal')
    set(ax, 'FontSize', 16, 'XLim', [xMin, xMax], 'YLim', [yMin, yMax], 'CLim', [log10(densityRange(1)) log10(densityRange(2))]);
    colormap(brewermap(1024, '*YlGnBu'));
    c = colorbar('Ticks',log10([1 3 10 30 100 300 1000 3000 10000 30000]),...
             'TickLabels',{'1', '3', '10','30','100','300','1000', '3000', '10000', '30000'});
    c.Label.String = sprintf('rfs / %s', densityUnits);
    xlabel(ax, sprintf('space (%s)', spatialSupportUnits));
    title(ax,titleLabel);
end

