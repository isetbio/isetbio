function optimizationPositionsAndSizesGrids = surroundOptimizationGrid(...
	peripheralSampling, minGridSize, maxGridSize, ...
	whichZernikeDataBase, whichEye, sourceLatticeSizeDegs, ...
	mosaicEccDegs, mosaicSizeDegs, varargin)

	% Parse optional input
    p = inputParser;
    p.addParameter('withExtremePositions',false, @isscalar)
    p.parse(varargin{:});
	addTheEightExtrementPositions = p.Results.withExtremePositions;

	assert(ismember(peripheralSampling, {'high res', 'medium res', 'low res', 'very low res'}), ...
		'invalid peripheralSampling. Valid values are: {''high res'', ''medium res'', ''low res'', ''very low res''}');

	neuronType = 'midget ganglion cells';
	
    mosaicOutline.x = mosaicEccDegs(1)+0.5*mosaicSizeDegs(1)*[-1 -1 1 1 -1];
    mosaicOutline.y = mosaicEccDegs(2)+0.5*mosaicSizeDegs(2)*[-1 1 1 -1 -1];
    [in, on] = inpolygon(0,0,mosaicOutline.x,mosaicOutline.y);
    mosaicContainsFovea = in | on;
    

    if (mosaicContainsFovea)

        switch (whichZernikeDataBase)
		    case 'Polans2015'
			    opticsXpositions = PolansOptics.constants.measurementHorizontalEccentricities;
			    opticsYpositions = PolansOptics.constants.measurementVerticalEccentricities;
		    case 'Artal2012'
			    opticsXpositions = ArtalOptics.constants.measurementHorizontalEccentricities;
			    opticsYpositions = ArtalOptics.constants.measurementVerticalEccentricities;
		    case 'Thibos'
			    opticsXpositions = 0;
			    opticsYpositions = 0;
		    otherwise
			    error('Unknown ZernikeDataBase: ''%s''.', whichZernikeDataBase);
    	end
        
        customMMsToDegsConversionFunction = @RGCmodels.Watson.convert.rhoMMsToDegs;

	    % Import full lattice positions up to +/- 30 degrees
	    maxEccDegsWidth = 68;
	    maxEccDegsHeight = 60;
	    eccMicrons = [0 0];
	    sizeMicrons = 1e3*RGCmodels.Watson.convert.rhoDegsToMMs([maxEccDegsWidth maxEccDegsHeight]);

	    fprintf('Importing full MRGC positions lattice (%d degs, for %s)\n', sourceLatticeSizeDegs, whichEye);
	    [~, mRGCRFposDegs] = retinalattice.import.finalMRGCPositions(...
         sourceLatticeSizeDegs, eccMicrons, sizeMicrons, ...
         whichEye, customMMsToDegsConversionFunction);

	    fprintf('Computing spacings for all MRGCs\n');
	    theRGCRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(mRGCRFposDegs);


        XLims = [-1 1] * (max(abs(squeeze(mRGCRFposDegs(:,1))))+0.5);
        YLims = [-1 1] * (max(abs(squeeze(mRGCRFposDegs(:,2))))+0.5);
   
        % Generate sampling grid based on cell density
	    [coneDensityXYGrid, neighboringSampleDistancesDegs] = generateSamplingGrid(mRGCRFposDegs, theRGCRFspacingsDegs);
	    
    
	    demoMode = ~true;
	    if (demoMode)
		    demoAllResolutions(opticsXpositions, opticsYpositions,  coneDensityXYGrid, neighboringSampleDistancesDegs, XLims, YLims);
        end

	    % Generate
	    optimizationPositionGrids = generateForPeripheralSampling(peripheralSampling, ...
		    opticsXpositions, opticsYpositions,  coneDensityXYGrid, neighboringSampleDistancesDegs, XLims, YLims);

    else
        % Generate
        optimizationPositionGrids = generateForUniformSampling(mosaicEccDegs, peripheralSampling);
    end

   % Only keep sampling positions within the extent of the mosaic
    xMin = mosaicEccDegs(1)-0.5*mosaicSizeDegs(1);
    xMax = mosaicEccDegs(1)+0.5*mosaicSizeDegs(1);
    yMin = mosaicEccDegs(2)-0.5*mosaicSizeDegs(2);
    yMax = mosaicEccDegs(2)+0.5*mosaicSizeDegs(2);

    idx = find(...
        (optimizationPositionGrids(:,1) >= xMin) & ...
        (optimizationPositionGrids(:,1) <= xMax) & ...
        (optimizationPositionGrids(:,2) >= yMin) & ...
        (optimizationPositionGrids(:,2) <= yMax) ...
        );
    optimizationPositionGrids = optimizationPositionGrids(idx,:);


    % Optimization grid positions
	optimizationPositionsAndSizesGrids(:,1:2) = optimizationPositionGrids;

	% Add the optimization grid sizes.
	optimizationPositionEccentricities = sqrt(sum(optimizationPositionGrids.^2,2));
	for iPos = 1:size(optimizationPositionGrids,1)
		r = optimizationPositionEccentricities(iPos)/max(optimizationPositionEccentricities);
		samplingRegionWidthDegs = minGridSize + (maxGridSize - minGridSize)*r;
		optimizationPositionsAndSizesGrids(iPos,3:4) = samplingRegionWidthDegs*[1 1];
	end

	if (addTheEightExtrementPositions)
		% Add eight exteme positions
		% Add the 8 extreme points
		lastSamplingGridSize = squeeze(optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1),3:4));

		% Top left
		optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1)+1,:) = ...
			[mosaicEccDegs(1)-0.5*mosaicSizeDegs(1) mosaicEccDegs(2)-0.5*mosaicSizeDegs(2) lastSamplingGridSize(1) lastSamplingGridSize(2)];
	    % Top middle
		optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1)+1,:) = ...
			[mosaicEccDegs(1) mosaicEccDegs(2)-0.5*mosaicSizeDegs(2) lastSamplingGridSize(1) lastSamplingGridSize(2)];
		% Top right
		optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1)+1,:) = ...
			[mosaicEccDegs(1)+0.5*mosaicSizeDegs(1) mosaicEccDegs(2)-0.5*mosaicSizeDegs(2) lastSamplingGridSize(1) lastSamplingGridSize(2)];

		% Left middle
		optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1)+1,:) = ...
			[mosaicEccDegs(1)-0.5*mosaicSizeDegs(1) mosaicEccDegs(2) lastSamplingGridSize(1) lastSamplingGridSize(2)];

		% Right middle
		optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1)+1,:) = ...
			[mosaicEccDegs(1)+0.5*mosaicSizeDegs(1) mosaicEccDegs(2) lastSamplingGridSize(1) lastSamplingGridSize(2)];

		% Bottom left
		optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1)+1,:) = ...
			[mosaicEccDegs(1)-0.5*mosaicSizeDegs(1) mosaicEccDegs(2)+0.5*mosaicSizeDegs(2) lastSamplingGridSize(1) lastSamplingGridSize(2)];

	    % Bottom middle
		optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1)+1,:) = ...
			[mosaicEccDegs(1) mosaicEccDegs(2)+0.5*mosaicSizeDegs(2) lastSamplingGridSize(1) lastSamplingGridSize(2)];
		% Bottom right
		optimizationPositionsAndSizesGrids(size(optimizationPositionsAndSizesGrids,1)+1,:) = ...
			[mosaicEccDegs(1)+0.5*mosaicSizeDegs(1) mosaicEccDegs(2)+0.5*mosaicSizeDegs(2) lastSamplingGridSize(1) lastSamplingGridSize(2)];
	end

end



function optimizationPositionGrids = generateForUniformSampling(mosaicEccDegs,  peripheralSampling)
    switch (peripheralSampling)
		case 'high res'
				opticalPositionIncrementDegs = 1.0;
		case 'medium res'
				opticalPositionIncrementDegs = 2;
		case 'low res'
				opticalPositionIncrementDegs = 3;
		case 'very low res'
				opticalPositionIncrementDegs = 4;
		otherwise
			error('Unknown peripheral sampling: ''%s''.', peripheralSampling)
    end

    xPos = mosaicEccDegs(1) + opticalPositionIncrementDegs*(-100:1:100);
    yPos = mosaicEccDegs(2) + opticalPositionIncrementDegs*(-100:1:100);
    [xPosGrid, yPosGrid] = meshgrid(xPos, yPos);
	optimizationPositionGrids = [xPosGrid(:) yPosGrid(:)];

end


function optimizationPositionGrids = generateForPeripheralSampling(peripheralSampling, ...
	opticsXpositions, opticsYpositions, coneDensityXYGrid,  neighboringSampleDistancesDegs, XLims, YLims)
	
	switch (peripheralSampling)
		case 'high res'
				coneDensityThresholdDistanceDegs = 1.0;
				opticalToConeDensityGridThresholdDistanceDegs = 0.5;
				opticalPositionIncrementDegs = 1.0;
		case 'medium res'
				coneDensityThresholdDistanceDegs = 2.0;
				opticalToConeDensityGridThresholdDistanceDegs = 1.0;
				opticalPositionIncrementDegs = 2;
		case 'low res'
				coneDensityThresholdDistanceDegs = 2.5;
				opticalToConeDensityGridThresholdDistanceDegs = 1.5;
				opticalPositionIncrementDegs = 3;
		case 'very low res'
				coneDensityThresholdDistanceDegs = 2.5;
				opticalToConeDensityGridThresholdDistanceDegs = 2.5;
				opticalPositionIncrementDegs = 4;
		otherwise
			error('Unknown peripheral sampling: ''%s''.', peripheralSampling)
	end

	% Generate opticalPositionXYGrid
	pos = 0:opticalPositionIncrementDegs:max(opticsXpositions);
	opticsXpositions = [-fliplr(pos) pos(2:end)];
	pos = 0:opticalPositionIncrementDegs:max(opticsYpositions);
	opticsYpositions = [-fliplr(pos) pos(2:end)];
	[opticalPositionXGrid, opticalPositionYGrid] = meshgrid(opticsXpositions, opticsYpositions);
	opticalPositionXYGrid = [opticalPositionXGrid(:) opticalPositionYGrid(:)];

	% Only keep cone density sampling points whose neigboring distance is less than coneDensityThresholdDistanceDegs away from each other
	idx = find(neighboringSampleDistancesDegs <= coneDensityThresholdDistanceDegs);
	coneDensityXYGrid = coneDensityXYGrid(idx,:);
	optimizationPositionGrids = coneDensityXYGrid;

	% Measure distances between coneDensityGrid and opticalPositionGrid
	[D,I] = pdist2(coneDensityXYGrid, opticalPositionXYGrid, 'euclidean','Smallest',1);
	% Only keep optical position grid points that are at least opticalToConeDensityGridThresholdDistanceDegs 
	% away from the cone density grid points
	idx = find(D>=opticalToConeDensityGridThresholdDistanceDegs);
	if (~isempty(idx))
		optimizationPositionGrids = cat(1,optimizationPositionGrids, opticalPositionXYGrid(idx,:));
	end

	if (1==2)
		figure(100); clf;
		plot(coneDensityXYGrid(:,1), coneDensityXYGrid(:,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5]);
		hold on
	    plot(opticalPositionXYGrid(:,1), opticalPositionXYGrid(:,2), 'bx', 'MarkerSize', 12, 'LineWidth', 1.5);
	    axis 'equal';
	    set(gca, 'XLim', XLims, 'YLim', YLims);
	end
end


function [allPoints, meanNeighborDistances] = generateSamplingGrid(RGCRFpositionsDegs, RGCRFspacingsDegs)

 	XLims = [-1 1] * (max(abs(squeeze(RGCRFpositionsDegs(:,1))))+0.5);
    YLims = [-1 1] * (max(abs(squeeze(RGCRFpositionsDegs(:,2))))+0.5);
   
    minSpacing = min(RGCRFspacingsDegs(:)) * 1.2;
    maxSpacing = max(RGCRFspacingsDegs(:)) * 0.95;

    spacingZones = logspace(log10(minSpacing), log10(maxSpacing), 30);
    zoneRelativeTolerance = 0.05;

    [theXcoords, theYcoords, theEllipseXaxes, theEllipseYaxes] = generateSamplingEllipses(...
    	RGCRFpositionsDegs(:,1),  ...
    	RGCRFpositionsDegs(:,2), ...
    	RGCRFspacingsDegs(:), ...
    	spacingZones, zoneRelativeTolerance);

    
    angles = 0:10:360;
    theCosines = cosd(angles);
    theSines = sind(angles);

    samplingAngles = 0:60:300;
    samplingPointsX = [];
    samplingPointsY = [];

    visualizeMeasurements = false;

    for iSamplingZone = 1:numel(spacingZones)
    	if (isnan(theEllipseXaxes(iSamplingZone)) || isnan(theEllipseYaxes(iSamplingZone)))
    		continue;
    	end
    	x = theEllipseXaxes(iSamplingZone) * theCosines;
	    y = theEllipseYaxes(iSamplingZone) * theSines;

	    idx = 1:numel(samplingAngles);
	    newSamplingXcoords = theEllipseXaxes(iSamplingZone) * cosd(samplingAngles(idx));
	    newSamplingYcoords = theEllipseYaxes(iSamplingZone) * sind(samplingAngles(idx));
	    if (any(isnan(newSamplingXcoords))) || (any(isnan(newSamplingYcoords)))
	    	continue;
	    end

	    samplingPointsX = cat(1, samplingPointsX, newSamplingXcoords(:));
	    samplingPointsY = cat(1, samplingPointsY, newSamplingYcoords(:));
	    samplingAngles = mod(samplingAngles + 30, 360);

	    if (visualizeMeasurements)
    		hFig2 = figure(1+iSamplingZone); clf;
    		ax = subplot(1,1,1);
    	
    		plot(ax,theXcoords{iSamplingZone}, theYcoords{iSamplingZone} , 'k.');
    		hold(ax,'on');
    		plot(ax,x,y, 'r-', 'LineWidth', 1.5);
    		
    		axis(ax, 'equal');
    		set(ax, 'XLim', XLims, 'YLim', YLims);
    		title(ax,sprintf('spacing:%2.2f (arc min)', spacingZones(iSamplingZone)*60));
    	end
    end

    allPoints(1,:) = [0 0];
    allPoints = cat(1, allPoints, [samplingPointsX(:) samplingPointsY(:)]);

    meanNeighborDistances = zeros(1, size(allPoints,1));
    for iPoint = 1:size(allPoints,1)
    	xo = allPoints(iPoint,:);
    	d = sort(sqrt(sum((bsxfun(@minus, allPoints, xo)).^2,2)), 'ascend');
    	meanNeighborDistances(iPoint) = min(d(2:7));
    end

end

function [theXcoords, theYcoords, theEllipseXaxes, theEllipseYaxes] = generateSamplingEllipses(...
	xCoords, yCoords, spacings, spacingZones, zoneRelativeTolerance)

    theEllipseXaxes = nan(1, numel(spacingZones));
    theEllipseYaxes = nan(1, numel(spacingZones));
    theXcoords = cell(1, numel(spacingZones));
    theYcoords = cell(1, numel(spacingZones));

    for i = 1:numel(spacingZones)
    	idx = find(abs(spacings-spacingZones(i)) < zoneRelativeTolerance * spacingZones(i));
    	if (idx > 1)
	    	theXcoords{i} = xCoords(idx);
	    	theYcoords{i} = yCoords(idx);
	    	theAngles = atan2d(yCoords(idx), xCoords(idx));
	    	theRadii = sqrt(xCoords(idx).^2 + yCoords(idx).^2);

	    	idx = find((abs(theAngles-0) < 45) | (abs(theAngles-180) < 45));
	    	theEllipseXaxes(i) = prctile(theRadii(idx), 50);

	    	idx = find((abs(theAngles-90) < 45) | (abs(theAngles+270) < 45));
	    	theEllipseYaxes(i) = prctile(theRadii(idx), 50);
	    end
    end
end


function demoAllResolutions(opticsXpositions, opticsYpositions,  coneDensityXYGrid, neighboringSampleDistancesDegs, XLims, YLims)
	peripheralSampling = 'high res';
	highResOptimizationPositionGrids = generateForPeripheralSampling(peripheralSampling, ...
		opticsXpositions, opticsYpositions,  coneDensityXYGrid, neighboringSampleDistancesDegs, XLims, YLims);

	peripheralSampling = 'medium res';
	mediumResOptimizationPositionGrids = generateForPeripheralSampling(peripheralSampling, ...
		opticsXpositions, opticsYpositions, coneDensityXYGrid, neighboringSampleDistancesDegs, XLims, YLims);

	peripheralSampling = 'low res';
	lowResOptimizationPositionGrids = generateForPeripheralSampling(peripheralSampling, ...
		opticsXpositions, opticsYpositions, coneDensityXYGrid, neighboringSampleDistancesDegs, XLims, YLims);

	peripheralSampling = 'very low res';
	veryLowResOptimizationPositionGrids = generateForPeripheralSampling(peripheralSampling, ...
		opticsXpositions, opticsYpositions, coneDensityXYGrid, neighboringSampleDistancesDegs, XLims, YLims);

	subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 2, ...
        'colsNum', 2, ...
        'heightMargin',  0.04, ...
        'widthMargin',    0.03, ...
        'leftMargin',     0.01, ...
        'rightMargin',    0.0, ...
        'bottomMargin',   0.04, ...
        'topMargin',      0.02);

	hFig = figure(100); clf;
	set(hFig, 'Position', [10 10 1800 1200]);

	ax = subplot('Position', subplotPosVectors(1,1).v);
    plot(ax,highResOptimizationPositionGrids(:,1), highResOptimizationPositionGrids(:,2), 'k.', 'MarkerSize', 10, 'LineWidth', 1.5);
    axis(ax,'equal');
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);
    title(ax, 'high-res');

    ax = subplot('Position', subplotPosVectors(1,2).v);
    plot(ax,mediumResOptimizationPositionGrids(:,1), mediumResOptimizationPositionGrids(:,2), 'k.', 'MarkerSize', 10, 'LineWidth', 1.5);
    axis(ax,'equal');
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);
    title(ax, 'medium-res');

    ax = subplot('Position', subplotPosVectors(2,1).v);
    plot(ax,lowResOptimizationPositionGrids(:,1), lowResOptimizationPositionGrids(:,2), 'k.', 'MarkerSize', 10, 'LineWidth', 1.5);
    axis(ax,'equal');
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);
    title(ax, 'low-res');

	ax = subplot('Position', subplotPosVectors(2,2).v);
    plot(ax,veryLowResOptimizationPositionGrids(:,1), veryLowResOptimizationPositionGrids(:,2), 'k.', 'MarkerSize', 10, 'LineWidth', 1.5);
    axis(ax,'equal');
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);
    title(ax, 'very-low res');
end