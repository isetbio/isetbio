function visualizeRFs(obj, rgcIndices, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandles', [], @(x)(isempty(x)||(iscell(x)&&(numel(x)==6))));
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    axesHandles = p.Results.axesHandles;
    fontSize = 16;

    for iRGC = 1:numel(rgcIndices)

        theTargetRGC = rgcIndices(iRGC);

        if (isempty(axesHandles))
            hFigFinal = figure(1000+iRGC); clf;
            set(hFigFinal, 'Color', [1 1 1], 'Position', [10 10 1200 800], 'Name', sprintf('RGC #%d', theTargetRGC));
            
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 2, ...
                'colsNum', 3, ...
                'heightMargin',  0.12, ...
                'widthMargin',    0.05, ...
                'leftMargin',     0.05, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.08, ...
                'topMargin',      0.02);
            axesHandlesFinal{1} = subplot('Position', subplotPosVectors(1,1).v);
            axesHandlesFinal{2} = subplot('Position', subplotPosVectors(2,2).v);
            axesHandlesFinal{3} = subplot('Position', subplotPosVectors(2,1).v);
            axesHandlesFinal{4} = subplot('Position', subplotPosVectors(1,2).v);
            
        else
            hFigFinal = hFig;
            axesHandlesFinal = axesHandles;
        end


        % The center cone connection data
        theRFcenterConeData = subregionConeData(obj, theTargetRGC,  'center');
        % The surround cone connection data
        theRFsurroundConeData = subregionConeData(obj, theTargetRGC, 'surround');

        [theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, ...
         theRetinalLineWeightingFunctions, spatialSupportDegs] = generateSpatialRFmaps(obj, theRFcenterConeData, theRFsurroundConeData);

        renderSubregionMap(axesHandlesFinal{1}, theRetinalRFcenterConeMap, theRFcenterConeData, spatialSupportDegs, fontSize);
        renderSubregionMap(axesHandlesFinal{2}, -theRetinalRFsurroundConeMap, theRFsurroundConeData, spatialSupportDegs, fontSize);
        renderLineWeightingFunction(axesHandlesFinal{3}, theRetinalLineWeightingFunctions, spatialSupportDegs, 'x', fontSize);
        renderLineWeightingFunction(axesHandlesFinal{4}, theRetinalLineWeightingFunctions, spatialSupportDegs, 'y', fontSize);
    end

end


function renderLineWeightingFunction(ax, theRetinalLineWeightingFunctions, spatialSupportDegs, whichAxis, fontSize)
    
    xSupportDegs = spatialSupportDegs(:,1);
    ySupportDegs = spatialSupportDegs(:,2);
    
    cMapEntries = 1024;
    cMap = (brewermap(cMapEntries, '*RdBu')).^1.0;

    mxMin = min(xSupportDegs);
    mxMax = max(xSupportDegs);
    myMin = min(ySupportDegs);
    myMax = max(ySupportDegs);

    x1 = ceil(mxMin/0.05)*0.05;
    x2 = floor(mxMax/0.05)*0.05;
    y1 = ceil(myMin/0.05)*0.05;
    y2 = floor(myMax/0.05)*0.05;

    sensitivityRange(1) = -1.05*max([max(theRetinalLineWeightingFunctions.surroundX(:)) max(theRetinalLineWeightingFunctions.surroundY(:))]);
    sensitivityRange(2) = 1.05*max([max(theRetinalLineWeightingFunctions.centerX(:)) max(theRetinalLineWeightingFunctions.centerY(:))]);
    sensitivityTicks = -1:0.05:1;

    switch (whichAxis)
        case 'x'
            shadedAreaBetweenTwoLines(ax, xSupportDegs', theRetinalLineWeightingFunctions.centerX, ...
                theRetinalLineWeightingFunctions.centerX*0, cMap(512+256,:), 'none', 0.3, 1.5, '-');
            hold(ax, 'on');
            shadedAreaBetweenTwoLines(ax, xSupportDegs', -theRetinalLineWeightingFunctions.surroundX, ...
                -theRetinalLineWeightingFunctions.surroundX*0, cMap(512-256,:), 'none', 0.3, 1.5, '-');
            plot(ax,xSupportDegs, theRetinalLineWeightingFunctions.centerX-theRetinalLineWeightingFunctions.surroundX, 'k-', 'LineWidth', 1.0);
            plot(ax,xSupportDegs, ySupportDegs*0, 'k--', 'LineWidth', 0.5);
            set(ax, 'XTick', x1: 0.05: x2, 'XLim', [mxMin mxMax], 'YLim', sensitivityRange, 'YTick', sensitivityTicks);
            xlabel(ax, 'space, x (degs)');
            ylabel(ax, 'gain');

        case 'y'
            shadedAreaBetweenTwoLines(ax, theRetinalLineWeightingFunctions.centerY', ySupportDegs', ...
                ySupportDegs'*0, cMap(512+256,:), 'none', 0.3, 1.5, '-');
            hold(ax, 'on');
            shadedAreaBetweenTwoLines(ax, -theRetinalLineWeightingFunctions.surroundY', ySupportDegs', ...
                ySupportDegs'*0, cMap(512-256,:), 'none', 0.3, 1.5, '-');
            plot(ax,theRetinalLineWeightingFunctions.centerY-theRetinalLineWeightingFunctions.surroundY, ...
                ySupportDegs, 'k-', 'LineWidth', 1.0);
            plot(ax,xSupportDegs*0, ySupportDegs, 'k--', 'LineWidth', 0.5);
            
            set(ax, 'Xdir', 'reverse', 'YAxisLocation', 'right');
            set(ax, 'YTick', y1:0.05:y2, 'YLim', [myMin myMax], 'XLim', sensitivityRange, 'XTick', sensitivityTicks);
            xlabel(ax, 'gain');


        otherwise
            error('Unknown axis, ''%s''.', whichAxis);
    end

    set(ax, 'FontSize', fontSize);
    grid(ax, 'on');
    box(ax, 'on');
    axis(ax, 'square');
    
    xtickangle(ax, 0);
end


function renderSubregionMap(ax, theSubregionConeMap, theConeData, spatialSupportDegs, fontSize)
    cMapEntries = 1024;
    cMap = (brewermap(cMapEntries, '*greys')).^1.0;
    xSupportDegs = spatialSupportDegs(:,1);
    ySupportDegs = spatialSupportDegs(:,2);

    mxMin = min(xSupportDegs);
    mxMax = max(xSupportDegs);
    myMin = min(ySupportDegs);
    myMax = max(ySupportDegs);

    x1 = ceil(mxMin/0.05)*0.05;
    x2 = floor(mxMax/0.05)*0.05;
    y1 = ceil(myMin/0.05)*0.05;
    y2 = floor(myMax/0.05)*0.05;

    imagesc(ax,xSupportDegs, ySupportDegs, theSubregionConeMap);
    hold(ax, 'on');

    for iInputCone = 1:size(theConeData.positionsDegs,1)
        switch theConeData.types(iInputCone)
            case cMosaic.LCONE_ID
                coneColor = [1 0 0];
            case cMosaic.MCONE_ID
                coneColor = [0 1 0];
            case cMosaic.SCONE_ID
                coneColor = [0 0 1];
        end
        plot(ax,theConeData.positionsDegs(iInputCone,1), theConeData.positionsDegs(iInputCone,2), '.', 'MarkerSize', 12, 'Color', coneColor);
    end

    hold(ax, 'off');
    axis(ax, 'xy');
    axis(ax, 'image'); 
    colormap(ax, cMap);
    set(ax, 'Color', cMap(1,:));
    set(ax, 'CLim',  0.8*max(abs(theSubregionConeMap(:)))*[-1 1]);
    set(ax, 'XTick', x1: 0.05: x2, 'YTick', y1:0.05:y2, 'XLim', [mxMin mxMax], 'YLim', [myMin myMax]);
    xlabel(ax, 'space, x (degs)');
    ylabel(ax, 'space, y (degs)');
    grid(ax, 'on');
    xtickangle(ax, 00);
    set(ax, 'FontSize', fontSize);

end


function [theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, ...
          theRetinalLineWeightingFunctions, spatialSupportDegs] = ...
                generateSpatialRFmaps(obj, theRFcenterConeData, theRFsurroundConeData)
    
    spatialSupportDegs(:,1) = -0.1:0.001:0.1;
    spatialSupportDegs(:,2) = spatialSupportDegs(:,1);

    rfCenterPositionDegs = mean(theRFcenterConeData.positionsDegs,1);
    spatialSupportDegs(:,1) = spatialSupportDegs(:,1) + rfCenterPositionDegs(1);
    spatialSupportDegs(:,2) = spatialSupportDegs(:,2) + rfCenterPositionDegs(2);

    theRetinalRFcenterConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
            obj.inputConeMosaic,...
            theRFcenterConeData.positionsDegs , ...
            theRFcenterConeData.poolingWeigts, ...
            spatialSupportDegs);

    theRetinalRFsurroundConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
            obj.inputConeMosaic,...
            theRFsurroundConeData.positionsDegs , ...
            theRFsurroundConeData.poolingWeigts, ...
            spatialSupportDegs);

    % Spatial profiles along x, y
    theRetinalLineWeightingFunctions.centerX = sum(theRetinalRFcenterConeMap,1);
    theRetinalLineWeightingFunctions.centerY = sum(theRetinalRFcenterConeMap,2);
    theRetinalLineWeightingFunctions.surroundX = sum(theRetinalRFsurroundConeMap,1);
    theRetinalLineWeightingFunctions.surroundY = sum(theRetinalRFsurroundConeMap,2);

    
end


function theConeData = subregionConeData(obj, theRGCindex, subregionName)

    switch (subregionName)
        case 'center'
            connectivityVector = full(squeeze(obj.centerConePoolingMatrix(:, theRGCindex)));
        case 'surround'
            connectivityVector = full(squeeze(obj.surroundConePoolingMatrix(:, theRGCindex)));
        otherwise
            error('Incorrect subregion: ''%s''.', subregionName);
    end

   minPoolingWeight = max(connectivityVector) * 0.001;
   subregionConeIndices = find(connectivityVector >= minPoolingWeight);

   theConeData.poolingWeigts = connectivityVector(subregionConeIndices);
   theConeData.positionsDegs = obj.inputConeMosaic.coneRFpositionsDegs(subregionConeIndices,:);
   theConeData.types = obj.inputConeMosaic.coneTypes(subregionConeIndices);
%    theConeData.Rc = ...
%         obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
%         obj.inputConeMosaic.coneApertureDiametersDegs(subregionConeIndices);

end

   
function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end
