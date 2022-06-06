function hFig = visualizeConePoolingWithinRFcenter(obj, iRGC, varargin)
    
     % Parse input
    p = inputParser;
    p.addParameter('visualizedFieldOfViewMicrons', [], @(x)((isempty(x))||(isscalar(x))));
    p.addParameter('visualizedConesNum', [], @(x)((isempty(x))||(isscalar(x))));
    p.addParameter('visualizedNeighbor', [], @(x)((isempty(x))||(isscalar(x)&&(x>0))));
    p.parse(varargin{:});
    visualizedFieldOfViewMicrons = p.Results.visualizedFieldOfViewMicrons;
    visualizedConesNum = p.Results.visualizedConesNum;
    visualizedNeighbor = p.Results.visualizedNeighbor;

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

    hFig = figure(); clf;
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
        obj.visualizeConeAperturePooling(iRGC, axConeWiringMain, axConeAperturesMain, ...
        [], [], [], [], visualizedFieldOfViewMicrons, visualizedConesNum, 'black');

    % Find the closest neighboring RGC
    nearbyRGCindices = obj.neihboringRGCindices(iRGC);
    if (isempty(nearbyRGCindices))
        return;
    end
    
    theNearbyRGC = nearbyRGCindices(min([visualizedNeighbor  numel(nearbyRGCindices)]));
    % Plot the nearby RGC wiring
    [~, ~, rfProfile2DnearbyRGC, ~, ~] = ...
        obj.visualizeConeAperturePooling(theNearbyRGC, axConeWiringNearby, axConeAperturesNearby,  ...
        xSupport, ySupport, xTicks, yTicks, visualizedFieldOfViewMicrons, visualizedConesNum, 'red');

    % Plot overlap between main and nearby RGC
    obj.visualizeRFoverlap2D(iRGC, theNearbyRGC, rfProfile2DmainRGC, rfProfile2DnearbyRGC, ... 
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
    RGCconnector.shadedAreaPlot(ax,spatialSupport, rfProfile1DmainRGC,   0, [0.3 0.3 0.3], [0. 0 0], 0.8, 1.5, '-'); 
    RGCconnector.shadedAreaPlot(ax,spatialSupport, rfProfile1DnearbyRGC, 0, [1 0 0], [1 0 0]*0.5, 0.5, 1.5, '-'); 
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