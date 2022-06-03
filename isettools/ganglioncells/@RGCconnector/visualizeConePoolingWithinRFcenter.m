function hFig = visualizeConePoolingWithinRFcenter(obj, iRGC, varargin)
    
     % Parse input
    p = inputParser;
    p.addParameter('visualizedFieldOfViewMicrons', [], @(x)((isempty(x))||(isscalar(x))));
    p.parse(varargin{:});
    visualizedFieldOfViewMicrons = p.Results.visualizedFieldOfViewMicrons;
    
    % Set-up figure
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', 3, ...
       'rowsNum', 3, ...
       'heightMargin',  0.09, ...
       'widthMargin',    0.06, ...
       'leftMargin',     0.06, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);

    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 1100 1000], 'Color', [1 1 1]);

    axConeWiringMain = subplot('Position', subplotPosVectors(1,1).v);
    axConeAperturesMain = subplot('Position', subplotPosVectors(2,1).v);

    axConeWiringNearby = subplot('Position', subplotPosVectors(1,2).v);
    axConeAperturesNearby = subplot('Position', subplotPosVectors(2,2).v);

    axRFOverlapPair2D = subplot('Position', subplotPosVectors(1,3).v);
    axRFOverlapPair1Dx = subplot('Position', subplotPosVectors(2,3).v);
    axRFOverlapPair1Dy = subplot('Position', subplotPosVectors(3,3).v);

    % Plot the main RGC wiring
    [xSupport, ySupport, rfProfile2DmainRGC, xTicks, yTicks] = ...
        plotWiring(obj, iRGC, axConeWiringMain, axConeAperturesMain, ...
        [], [], [], [], visualizedFieldOfViewMicrons, 'black');

    % Find the closest neighboring RGC
    nearbyRGCindices = obj.neihboringRGCindices(iRGC);
    nearbyRGCindices = nearbyRGCindices(1:min([3 numel(nearbyRGCindices)]));
    if (isempty(nearbyRGCindices))
        return;
    end
    
    theNearbyRGC = nearbyRGCindices(1);

    % Compute the overlap coefficient
    weightsRGC = abs(full(obj.coneConnectivityMatrix(:, iRGC)));
    weightsNearbyRGC = abs(full(obj.coneConnectivityMatrix(:, theNearbyRGC)));
    overlapOfWeights = RGCconnector.overlap(weightsRGC, weightsNearbyRGC);
    
    
    % Plot the nearby RGC wiring
    [~, ~, rfProfile2DnearbyRGC, ~, ~] = ...
        plotWiring(obj, theNearbyRGC, axConeWiringNearby, axConeAperturesNearby,  ...
        xSupport, ySupport, xTicks, yTicks, visualizedFieldOfViewMicrons, 'red');

    % Plot overlap between main and nearby RGC
    plotRFoverlaps2D(obj, iRGC, theNearbyRGC, rfProfile2DmainRGC, rfProfile2DnearbyRGC, ... 
        xSupport, ySupport, xTicks, yTicks, axRFOverlapPair2D);


    % The horizontal and vertical line spread functions (integrated over
    % corresponding axis)

    theLegends = {sprintf('RGC #%d', iRGC), sprintf('RGC #%d',theNearbyRGC)};

    dX = xSupport(2)-xSupport(1);
    rfProfile1DXmainRGC = sum(rfProfile2DmainRGC,1)*dX;
    rfProfile1DXnearbyRGC = sum(rfProfile2DnearbyRGC,1)*dX;
    dY = ySupport(2)-ySupport(1);
    rfProfile1DYmainRGC = sum(rfProfile2DmainRGC,2)*dY;
    rfProfile1DYnearbyRGC = sum(rfProfile2DnearbyRGC,2)*dY;

    maxY = max([max(rfProfile1DXmainRGC(:)) max(rfProfile1DXnearbyRGC(:)) max(rfProfile1DYmainRGC(:)) max(rfProfile1DYnearbyRGC(:))]);

    renderLineSpreadDiagrams(axRFOverlapPair1Dx, xSupport, rfProfile1DXmainRGC, rfProfile1DXnearbyRGC, theLegends, xTicks, maxY, 'space, x (microns)');
    renderLineSpreadDiagrams(axRFOverlapPair1Dy, ySupport, rfProfile1DYmainRGC, rfProfile1DYnearbyRGC, theLegends, yTicks, maxY, 'space, y (microns)');
    

end

function renderLineSpreadDiagrams(ax, spatialSupport, rfProfile1DmainRGC, rfProfile1DnearbyRGC, theLegends, xyTicks, maxY, xLabelString)
    hold(ax, 'on');
    shadedAreaPlot(ax,spatialSupport, rfProfile1DmainRGC,   0, [0.3 0.3 0.3], [0. 0 0], 0.8, 1.5, '-'); 
    shadedAreaPlot(ax,spatialSupport, rfProfile1DnearbyRGC, 0, [1 0 0], [1 0 0]*0.5, 0.5, 1.5, '-'); 
    legend(ax,theLegends);

    overlapCoeff = RGCconnector.overlap(rfProfile1DmainRGC, rfProfile1DnearbyRGC);

    % Finalize plot
    axis(ax,'square');
    set(ax, 'XLim', [spatialSupport(1) spatialSupport(end)], 'YLim', [0 maxY], 'FontSize', 15);
    set(ax, 'XTick', xyTicks);
    if (contains(xLabelString, 'x'))
        title(ax, sprintf('1D-overlap across y-integrated RF: %2.0f%%', overlapCoeff*100));
    else
        title(ax, sprintf('1D-overlap across x-integrated RF: %2.0f%%', overlapCoeff*100));

    end
    box(ax, 'on'); grid(ax, 'on'); 
    xlabel(ax, xLabelString);
    ylabel(ax, 'sensitivity');
end


function plotRFoverlaps2D(obj, iRGC, nearbyRGCindex, ...
    rfProfile2DmainRGC, rfProfile2DnearbyRGC, ...
    xSupport, ySupport, xTicks, yTicks, axRFOverlap2D)

    zLevels = [0.07 0.999];

    cmapMain = brewermap(1024, 'greys');
    alphaMain = 0.95;
    contourLineColorMain = [0 0 0];

    cmapNearby = brewermap(1024, 'reds');
    alphaNearby = 0.4;
    contourLineColorNearby = [1 0 0];
    hold(axRFOverlap2D, 'on');


    zData1 = (rfProfile2DmainRGC/max(rfProfile2DmainRGC(:))).^0.5;
    cMosaic.semiTransparentContourPlot(axRFOverlap2D, xSupport, ySupport, ...
        zData1, zLevels, cmapMain, alphaMain, contourLineColorMain, ...
        'lineWidth', 2.0);
    

    zData2 = (rfProfile2DnearbyRGC/max(rfProfile2DnearbyRGC(:))).^0.5;
    cMosaic.semiTransparentContourPlot(axRFOverlap2D, xSupport, ySupport, ...
              zData2, zLevels, cmapNearby, alphaNearby, contourLineColorNearby, ...
              'lineWidth', 2.0);

    % Compute the overlap coefficient
    overlapCoeff = RGCconnector.overlap(zData1(:), zData2(:));
    
    
    % Finalize plot
    axis(axRFOverlap2D,'equal');
    set(axRFOverlap2D, 'XLim', [xSupport(1) xSupport(end)], 'YLim', [ySupport(1) ySupport(end)], 'FontSize', 15);
    set(axRFOverlap2D, 'XTick', xTicks, 'YTick', yTicks);
    box(axRFOverlap2D, 'on'); grid(axRFOverlap2D, 'on');
    title(axRFOverlap2D, sprintf('Rc/RGCsep: %1.2f, 2D-overlap: %2.0f%%', ...
        obj.wiringParams.RcToRGCseparationRatio, 100*overlapCoeff));
    xlabel(axRFOverlap2D, 'space, x (microns)');
    ylabel(axRFOverlap2D, 'space, y (microns)');

end



function [xSupport, ySupport, rfProfile2D, xTicks, yTicks] = ...
    plotWiring(obj, iRGC, axConeWiring, axConeApertures, ...
    xSupport, ySupport, xTicks, yTicks, visualizedFieldOfViewMicrons, colorString)

    % Indices and weights of non-overlapping & overlapping input cones
    connectedNonOverlappingConeIndices = find(squeeze(obj.coneConnectivityMatrix(:, iRGC))>0);
    connectedOverlappingConeIndices = find(squeeze(obj.coneConnectivityMatrix(:, iRGC))<0);
    connectedNonOverlappingConeWeights = full(obj.coneConnectivityMatrix(connectedNonOverlappingConeIndices, iRGC));
    connectedOverlappingConeWeights = full(obj.coneConnectivityMatrix(connectedOverlappingConeIndices, iRGC));

    % The centroid of the RGC
    theRGCCentroidMicrons = obj.RGCRFcentroidsFromInputs(iRGC,:);


    % Find 100 neigboring cones that are not-connected to this RGC
    [~, nearbyConeIndices] = RGCconnector.pdist2(...
            obj.inputConeMosaic.coneRFpositionsMicrons, ...
            theRGCCentroidMicrons, '', ...
            'smallest', 100);

    allConnectedConeIndices = [ ...
        connectedNonOverlappingConeIndices(:); ...
        connectedOverlappingConeIndices(:)
        ];

    neighboringNonConnectedConeIndices = setdiff(nearbyConeIndices, allConnectedConeIndices);

    % The cone outline
    coneOutline(:,1) = cosd(0:10:360);
    coneOutline(:,2) = sind(0:10:360);
    

    % The cone wiring diagram with weights
    hold(axConeWiring, 'on');
    [xSupport, ySupport] = renderPooledConesDiagram(axConeWiring, ...
        theRGCCentroidMicrons, ...
        obj.inputConeMosaic.coneRFpositionsMicrons(neighboringNonConnectedConeIndices,:), ...
        obj.inputConeMosaic.coneRFspacingsMicrons(neighboringNonConnectedConeIndices), ...
        obj.inputConeMosaic.coneTypes(neighboringNonConnectedConeIndices), ...
        obj.inputConeMosaic.coneRFpositionsMicrons(connectedNonOverlappingConeIndices, :),...
        obj.inputConeMosaic.coneRFpositionsMicrons(connectedOverlappingConeIndices, :),...
        obj.inputConeMosaic.coneRFspacingsMicrons(connectedNonOverlappingConeIndices), ...
        obj.inputConeMosaic.coneRFspacingsMicrons(connectedOverlappingConeIndices), ...
        obj.inputConeMosaic.coneTypes(connectedNonOverlappingConeIndices), ...
        obj.inputConeMosaic.coneTypes(connectedOverlappingConeIndices), ...
        connectedNonOverlappingConeWeights, ...
       -connectedOverlappingConeWeights, ...
        coneOutline, ...
        xSupport, ySupport, ...
        visualizedFieldOfViewMicrons ...
        );

    % Compute the 2D profile
    [X,Y] = meshgrid(xSupport, ySupport);
    rfProfile2D = computeRFprofile(...
        obj.inputConeMosaic.coneRFpositionsMicrons(connectedNonOverlappingConeIndices, :),...
        obj.inputConeMosaic.coneRFpositionsMicrons(connectedOverlappingConeIndices, :),...
        obj.inputConeMosaic.coneRFspacingsMicrons(connectedNonOverlappingConeIndices), ...
        obj.inputConeMosaic.coneRFspacingsMicrons(connectedOverlappingConeIndices), ...
        connectedNonOverlappingConeWeights, ...
       -connectedOverlappingConeWeights, ...
        X,Y);

    % Superimpose contour of 2Dprofile
    zLevels = [0.07 0.999];
    cmapMain = brewermap(1024, 'greys');
    alphaMain = 0.25;
    contourLineColorMain = [0 0 0];

    zData = (rfProfile2D/max(rfProfile2D(:))).^0.5;
    cMosaic.semiTransparentContourPlot(axConeWiring, xSupport, ySupport, ...
        zData, zLevels, cmapMain, alphaMain, contourLineColorMain, ...
        'lineWidth', 2.0);


    % Finalize plot
    if (isempty(xTicks))
        if (visualizedFieldOfViewMicrons < 30)
            deltaTick = 5;
        elseif (visualizedFieldOfViewMicrons < 60)
            deltaTick = 10;
        elseif (visualizedFieldOfViewMicrons < 100)
            deltaTick = 20;
        else
            deltaTick = 40;
        end
        
        xTicks = round(theRGCCentroidMicrons(1)/10)*10+(-visualizedFieldOfViewMicrons:deltaTick:visualizedFieldOfViewMicrons);
        yTicks = round(theRGCCentroidMicrons(2)/10)*10+(-visualizedFieldOfViewMicrons:deltaTick:visualizedFieldOfViewMicrons);
    end

    axis(axConeWiring,'equal');
    set(axConeWiring, 'XLim', [xSupport(1) xSupport(end)], 'YLim', [ySupport(1) ySupport(end)], 'FontSize', 15);
    set(axConeWiring, 'XTick', xTicks, ...
            'YTick', yTicks);
    box(axConeWiring, 'on'); grid(axConeWiring, 'on');
    title(axConeWiring, sprintf('\\color{%s} RGC #%d (connected cones)', colorString, iRGC));   
    xlabel(axConeWiring, 'space, x (microns)');
    ylabel(axConeWiring, 'space, y (microns)');
    
    % The weighted cone apertures
    hold(axConeApertures, 'on');
    renderWeightedConeAperturesDiagram(axConeApertures, ...
        obj.inputConeMosaic.coneRFpositionsMicrons(connectedNonOverlappingConeIndices, :),...
        obj.inputConeMosaic.coneRFpositionsMicrons(connectedOverlappingConeIndices, :),...
        obj.inputConeMosaic.coneRFspacingsMicrons(connectedNonOverlappingConeIndices), ...
        obj.inputConeMosaic.coneRFspacingsMicrons(connectedOverlappingConeIndices), ...
        obj.inputConeMosaic.coneTypes(connectedNonOverlappingConeIndices), ...
        obj.inputConeMosaic.coneTypes(connectedOverlappingConeIndices), ...
        connectedNonOverlappingConeWeights, ...
       -connectedOverlappingConeWeights, ...
        xSupport, ySupport);

    % Finalize plot
    axis(axConeApertures,'square');
    set(axConeApertures, 'XLim', [xSupport(1) xSupport(end)], 'YLim', [0 1.02], 'FontSize', 15);
    set(axConeApertures, 'XTick', xTicks, ...
            'YTick', 0:0.1:1.0);
    box(axConeApertures, 'on'); grid(axConeApertures, 'on');
    title(axConeApertures, sprintf('\\color{%s} RGC #%d (weighted cone apertures)', colorString, iRGC)); 
    xlabel(axConeApertures, 'space, x (microns)');
    ylabel(axConeApertures, 'pooling weight');
end


function rfProfile2D = computeRFprofile(...
        connectedNonOverlappingConeRFpositionsMicrons,...
        connectedOverlappingConeRFpositionsMicrons,...
        connectedNonOverlappingConeRFspacingsMicrons,...
        connectedOverlappingConeRFspacingsMicrons,...
        connectedNonOverlappingConeWeights, ...
        connectedOverlappingConeWeights, ...
        X,Y)

    allInputConePositions = [...
        connectedNonOverlappingConeRFpositionsMicrons;
        connectedOverlappingConeRFpositionsMicrons];

    allInputConeSpacings = [...
        connectedNonOverlappingConeRFspacingsMicrons(:); ...
        connectedOverlappingConeRFspacingsMicrons(:)];

    allConeWeights = [ ...
        connectedNonOverlappingConeWeights(:); ...
        connectedOverlappingConeWeights(:)
        ];


    rfProfile2D = [];
    for iCone = 1:size(allInputConePositions,1)
        coneRc = 0.204*sqrt(2.0)*allInputConeSpacings(iCone);
        gaussianProfile2D = exp(-((X-allInputConePositions(iCone,1))/coneRc).^2) .* ...
                            exp(-((Y-allInputConePositions(iCone,2))/coneRc).^2);

        gaussianProfile2D(gaussianProfile2D<0.001) = 0;
        if (isempty(rfProfile2D))
            rfProfile2D = abs(allConeWeights(iCone))*gaussianProfile2D;
        else
            rfProfile2D = rfProfile2D + abs(allConeWeights(iCone))*gaussianProfile2D;
        end
    end
end


function renderWeightedConeAperturesDiagram(ax, ...
        connectedNonOverlappingConeRFpositionsMicrons,...
        connectedOverlappingConeRFpositionsMicrons,...
        connectedNonOverlappingConeRFspacingsMicrons,...
        connectedOverlappingConeRFspacingsMicrons,...
        connectedNonOverlappingConeTypes, ...
        connectedOverlappingConeTypes, ...
        connectedNonOverlappingConeWeights, ...
        connectedOverlappingConeWeights, ...
        xSupport, ySupport)


    % The non-overlapping cone apertures
    renderWeightedConeApertures(ax, xSupport, ...
        connectedNonOverlappingConeRFpositionsMicrons, ...
        connectedNonOverlappingConeRFspacingsMicrons, ...
        connectedNonOverlappingConeTypes, ...
        connectedNonOverlappingConeWeights, ...
        '-');

    % The overlapping cone apertures
    renderWeightedConeApertures(ax, xSupport, ...
        connectedOverlappingConeRFpositionsMicrons, ...
        connectedOverlappingConeRFspacingsMicrons, ...
        connectedOverlappingConeTypes, ...
        connectedOverlappingConeWeights, ...
        '--');

end

function renderWeightedConeApertures(ax, xSupport, coneRFpositions, ...
    coneRFspacings, coneTypes, coneWeights, lineStyle)

    for iCone = 1:numel(coneRFspacings)
        switch (coneTypes(iCone))
            case cMosaic.LCONE_ID
                coneColor = [1 0.1000 0.5000];
            case cMosaic.MCONE_ID
                coneColor = [0.1000 1 0.5000];
            case cMosaic.SCONE_ID
                coneColor = [0.6000 0.1000 1];
        end
        coneRc = 0.204*sqrt(2.0)*coneRFspacings(iCone);
        gaussianProfile = exp(-((xSupport-coneRFpositions(iCone,1))/coneRc).^2);
        shadedAreaPlot(ax, xSupport, coneWeights(iCone)*gaussianProfile, 0, ...
            coneColor*0.5, [0 0 0], 0.5, 1.0, lineStyle); 
    end

end


function [xSupport, ySupport] = renderPooledConesDiagram(ax, ...
        theRGCCentroidMicrons, ...
        neighboringNonConnectedConeRFpositionsMicrons, ...
        neighboringNonConnectedConeRFspacingsMicrons, ...
        neighboringNonConnectedConeTypes, ...
        connectedNonOverlappingConeRFpositionsMicrons,...
        connectedOverlappingConeRFpositionsMicrons,...
        connectedNonOverlappingConeRFspacingsMicrons,...
        connectedOverlappingConeRFspacingsMicrons,...
        connectedNonOverlappingConeTypes, ...
        connectedOverlappingConeTypes, ...
        connectedNonOverlappingConeWeights, ...
        connectedOverlappingConeWeights, ...
        coneOutline, ...
        xSupport, ySupport, ...
        visualizedFieldOfViewMicrons ...
        )
    
    allInputConePositions = [...
        connectedNonOverlappingConeRFpositionsMicrons;
        connectedOverlappingConeRFpositionsMicrons];
    
    meanConeSpacing = mean([connectedNonOverlappingConeRFspacingsMicrons connectedOverlappingConeRFspacingsMicrons]);
    minXY = min(allInputConePositions,[],1);
    maxXY = max(allInputConePositions,[],1);
    rangeXY = maxXY-minXY;
    if (isempty(visualizedFieldOfViewMicrons))
        halfWidthMicrons = 1.25*(round(max(rangeXY(:))/2+meanConeSpacing));
    else
        halfWidthMicrons = 0.5*visualizedFieldOfViewMicrons;
    end
    
    % Spatial support
    if (isempty(xSupport))
        xSupport = linspace(theRGCCentroidMicrons(1)-halfWidthMicrons, theRGCCentroidMicrons(1)+halfWidthMicrons, 201);
        ySupport = linspace(theRGCCentroidMicrons(2)-halfWidthMicrons, theRGCCentroidMicrons(2)+halfWidthMicrons, 201);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
   
    % The non-overlapping cones
    renderCones(ax, connectedNonOverlappingConeRFpositionsMicrons, ...
        connectedNonOverlappingConeRFspacingsMicrons, ...
        connectedNonOverlappingConeTypes, ...
        coneOutline, '-', 0.5);

    % The overlapping cones
    renderCones(ax, connectedOverlappingConeRFpositionsMicrons, ...
        connectedOverlappingConeRFspacingsMicrons, ...
        connectedOverlappingConeTypes, ...
        coneOutline, '--', 0.5);

    % The neighboring non-connected cones
    renderCones(ax, neighboringNonConnectedConeRFpositionsMicrons, ...
        neighboringNonConnectedConeRFspacingsMicrons, ...
        neighboringNonConnectedConeTypes, ...
        coneOutline, 'none', 0.2);

    % The pooling weights (lines)
    if (1==2)
        renderPoolingWeightLines(ax, theRGCCentroidMicrons, ...
            connectedNonOverlappingConeRFpositionsMicrons, ...
            connectedNonOverlappingConeWeights);
        renderPoolingWeightLines(ax, theRGCCentroidMicrons, ...
            connectedOverlappingConeRFpositionsMicrons, ...
            connectedOverlappingConeWeights);
    end


    % The pooling weights (text)
    renderPoolingWeightText(ax, theRGCCentroidMicrons, ...
        connectedNonOverlappingConeRFpositionsMicrons, ...
        connectedNonOverlappingConeWeights, meanConeSpacing, ...
        xSupport, ySupport);
    renderPoolingWeightText(ax, theRGCCentroidMicrons, ...
        connectedOverlappingConeRFpositionsMicrons, ...
        connectedOverlappingConeWeights, meanConeSpacing, ...
        xSupport, ySupport);

    drawnow;
end

function renderPoolingWeightText(ax, centroidPosition, coneRFpositions, coneWeights, meanConeSpacing, xSupport, ySupport)
    
    
    for iCone = 1:size(coneRFpositions,1)
        dy = (centroidPosition(2)-coneRFpositions(iCone,2));
        dx = (centroidPosition(1)-coneRFpositions(iCone,1));
        theta = atan2d(dy,dx);
        
        if (coneWeights(iCone) > 0.999)
          theText = sprintf('1');
          dd = meanConeSpacing/6;
        elseif (coneWeights(iCone)>0.1)
          theText = sprintf('%.1f', coneWeights(iCone));
          dd = meanConeSpacing/4;
        else
           dd = meanConeSpacing/6;
           theText =  sprintf('%.2f', coneWeights(iCone));
        end
        xx = coneRFpositions(iCone,1)-dd;
        yy = coneRFpositions(iCone,2);
        
        if ( ...
            (xx > xSupport(1)) && ...
            (xx < xSupport(end)) && ...
            (yy > ySupport(1)) && ...
            (yy < ySupport(end)))
            text(ax, xx, yy, theText, 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0 0 0]); %, 'BackgroundColor', [0.2 0.2 0.2]);
        end
        
    end
end

function renderPoolingWeightLines(ax, centroidPosition, coneRFpositions, coneWeights)
    for iCone = 1:size(coneRFpositions,1)
        xx = [centroidPosition(1) coneRFpositions(iCone,1)];
        yy = [centroidPosition(2) coneRFpositions(iCone,2)];
        plot(ax, xx, yy, 'k-', 'LineWidth', max([0.5 coneWeights(iCone)*4]));
    end
end

function renderCones(ax, coneRFpositionsMicrons, coneRFspacingsMicrons, coneTypes, diskOutline, lineStyle, faceAlpha)

    for iCone = 1:numel(coneTypes)
        xx = coneRFpositionsMicrons(iCone,1) + 0.45*coneRFspacingsMicrons(iCone)*diskOutline(:,1);
        yy = coneRFpositionsMicrons(iCone,2) + 0.45*coneRFspacingsMicrons(iCone)*diskOutline(:,2);
        switch (coneTypes(iCone))
            case cMosaic.LCONE_ID
                coneColor = [1 0.1000 0.5000];
            case cMosaic.MCONE_ID
                coneColor = [0.1000 1 0.5000];
            case cMosaic.SCONE_ID
                coneColor = [0.6000 0.1000 1];
        end
        pz = -10*eps*ones(size(yy)); 
        patch(ax, xx, yy, pz, coneColor*0.3, ...
            'FaceColor', coneColor*0.5, 'EdgeColor', [0 0 0], ...
            'FaceAlpha', faceAlpha, 'LineWidth', 1.5, 'LineStyle', lineStyle); 
    end
end



function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x fliplr(x)];
    y = [y y*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end

