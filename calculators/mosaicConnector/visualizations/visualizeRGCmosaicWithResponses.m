function hFig = visualizeRGCmosaicWithResponses(figNo,theConeMosaic, xAxisScaling, ...
    xAxisData, theMidgetRGCmosaicResponses, eccentricityMicrons, sizeMicrons, ...
    theMidgetRGCmosaic,  zLevels, subregions)

    % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center so as to align with ecc-varying full mRGC mosaic
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, eccentricityMicrons);
   
    % Cone types
    coneTypes = cmStruct.coneTypes;
    
    % Cone apertures and spacing
    coneDiameterMicrons = cmStruct.coneApertures * theConeMosaic.micronsPerDegree;
    coneSpacingsMicrons = 1.0/0.7 * coneDiameterMicrons;
    
    % Sampling for RF center contours
    deltaX = 0.25;
    extraMicrons = 1.1*theMidgetRGCmosaic.extraMicronsForSurroundCones;
    xAxis = (eccentricityMicrons(1)-sizeMicrons(1)/2-extraMicrons): deltaX: (eccentricityMicrons(1)+sizeMicrons(1)/2+extraMicrons);
    yAxis = (eccentricityMicrons(2)-sizeMicrons(2)/2-extraMicrons): deltaX: (eccentricityMicrons(2)+sizeMicrons(2)/2+extraMicrons);
    
   
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [1 1 2040 950], 'Color', [1 1 1]);
        
%     % Full view, no cone labeling
%     renderPlot(theAxesGrid{1,1}, theMidgetRGCmosaicResponses, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
%         theMidgetRGCmosaic, eccentricityMicrons, xAxis, yAxis, zLevels, subregions, false);
%     xLims = [xAxis(1) xAxis(end)];
%     yLims = [yAxis(1) yAxis(end)];
%     set(theAxesGrid{1,1}, 'XLim', xLims, 'YLim', yLims);
    
    % Zoomedin view, labeled cone types
    theAxes = axes('Position', [0.02 0.03 0.45 0.96]);
    renderPlot(theAxes, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
        theMidgetRGCmosaic, eccentricityMicrons, xAxis, yAxis, zLevels, subregions, true);
    xLims = eccentricityMicrons(1)+1.1*sizeMicrons(1)/2.0*[-1 1];
    yLims = eccentricityMicrons(2)+1.1*sizeMicrons(2)/2.0*[-1 1];
    axis(theAxes, 'equal');
    box(theAxes, 'on');
    set(theAxes, 'XLim', xLims, 'YLim', yLims, 'FontSize', 12);

    
    rgcsNum = size(theMidgetRGCmosaic.centerWeights,2);
    RGCpositions = determineResponsePositions(theConeMosaic, eccentricityMicrons, theMidgetRGCmosaic.centerWeights);
    
    w = 0.1/2;
    h = 0.1/2;
    gw = 0.5-w;
    gh = 0.96-h;
    m = 0.015/2;

    for iRGC = 1:rgcsNum
            ax = axes('Position', [0.48+gw*RGCpositions(iRGC,1)+m 0.03+gh*RGCpositions(iRGC,2)+m w h]);
            line(ax, xAxisData, squeeze(theMidgetRGCmosaicResponses(iRGC,:)), 'Color', [0 0 0], 'LineWidth', 1.5);
            set(ax, 'XTickLabel', {}, 'YTickLabel', {}, ...
                'XLim', [xAxisData(1) xAxisData(end)], 'YLim', [min(theMidgetRGCmosaicResponses(:)) max(theMidgetRGCmosaicResponses(:))], ...
                'XColor', 'none', 'YColor', 'none');
            set(ax, 'XScale', xAxisScaling);
            box(ax, 'off'); grid(ax, 'off');
            axis(ax, 'square')
            drawnow;
    end
        
end


function RGCpositions = determineResponsePositions(theConeMosaic, eccentricityMicrons, centerWeights)
     % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center so as to align with ecc-varying full mRGC mosaic
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, eccentricityMicrons);
    
    rgcsNum = size(centerWeights,2);
    RGCpositions = zeros(rgcsNum,2);
    for iRGC = 1:rgcsNum
        weights = full(squeeze(centerWeights(:, iRGC)));
        centerIndices = find(weights>0);
        RGCpositions(iRGC,:) = mean(conePositionsMicrons(centerIndices,:),1);
    end
    
    % Normalize to [0..1]
    mins = min(RGCpositions, [], 1);
    maxs = max(RGCpositions, [], 1);
    RGCpositions(:,1) = (RGCpositions(:,1)-mins(1))/(maxs(1)-mins(1));
    RGCpositions(:,2) = (RGCpositions(:,2)-mins(2))/(maxs(2)-mins(2));
end


function renderPlot(theAxes, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
    theMidgetRGCmosaic, eccentricityMicrons, xAxis, yAxis, zLevels, subregions, coneLabeling)

    global LCONE_ID
    global MCONE_ID
    global SCONE_ID

    [X,Y] = meshgrid(xAxis,yAxis);
    
%     % Render the optical image
%     spatialSupportMM = oiGet(theOI, 'spatial support', 'mm');
%     spatialSupportMicrons = 1e3 * squeeze(spatialSupportMM(1,:,1));
%     imagesc(theAxes, ...
%             spatialSupportMicrons+eccentricityMicrons(1), ...
%             spatialSupportMicrons+eccentricityMicrons(2), ...
%             oiGet(theOI, 'rgbimage'));
%     axis(theAxes, 'xy')
%     axis(theAxes, 'image')

    hold(theAxes, 'on');
    colormap(theAxes, brewermap(512, 'greys'));
    
    xxx = cosd(0:10:360);
    yyy = sind(0:10:360);
    
    if (coneLabeling)
        % Display cones
        LconeIndices = find(coneTypes == LCONE_ID);
        MconeIndices = find(coneTypes == MCONE_ID);
        SconeIndices = find(coneTypes == SCONE_ID);
    
        for k = 1:numel(LconeIndices)
            xx = conePositionsMicrons(LconeIndices(k),1) + 0.5*coneDiameterMicrons(LconeIndices(k))*xxx;
            yy = conePositionsMicrons(LconeIndices(k),2) + 0.5*coneDiameterMicrons(LconeIndices(k))*yyy;
            line(theAxes, xx, yy, 'Color', 'r', 'LineWidth', 1.5);
        end
        for k = 1:numel(MconeIndices)
            xx = conePositionsMicrons(MconeIndices(k),1) + 0.5*coneDiameterMicrons(MconeIndices(k))*xxx;
            yy = conePositionsMicrons(MconeIndices(k),2) + 0.5*coneDiameterMicrons(MconeIndices(k))*yyy;
            line(theAxes, xx, yy, 'Color', [0 0.8 0], 'LineWidth', 1.5);
        end
        for k = 1:numel(SconeIndices)
            xx = conePositionsMicrons(SconeIndices(k),1) + 0.5*coneDiameterMicrons(SconeIndices(k))*xxx;
            yy = conePositionsMicrons(SconeIndices(k),2) + 0.5*coneDiameterMicrons(SconeIndices(k))*yyy;
            line(theAxes, xx, yy, 'Color', 'b', 'LineWidth', 1.5);
        end
    else
        for k = 1:size(conePositionsMicrons,1)
            xx = conePositionsMicrons(k,1) + 0.5*coneDiameterMicrons(k)*cosd(0:30:360);
            yy = conePositionsMicrons(k,2) + 0.5*coneDiameterMicrons(k)*sind(0:30:360);
            line(theAxes, xx, yy, 'Color', 'k');
        end
    end
    
    rgcsNum = size(theMidgetRGCmosaic.centerWeights,2);
    for mRGCindex = 1:rgcsNum
        centerWeights = full(squeeze(theMidgetRGCmosaic.centerWeights(:, mRGCindex)));
        centerIndices = find(centerWeights>0);
        centerWeights(centerWeights>0) = 1;
                
        % Generate RF centers of RGCs based on cone positions and connection matrix
        switch subregions
            case 'centers'        
                theRF = generateRGCRFcenterSubregionFromConnectivityMatrix(...
                    centerWeights, conePositionsMicrons, coneSpacingsMicrons, X,Y);
            case 'surrounds'
                surroundWeights = full(theMidgetRGCmosaic.surroundWeights(:, mRGCindex));
                theRF = generateRGCRFcenterSubregionFromConnectivityMatrix(...
                    surroundWeights, conePositionsMicrons, coneSpacingsMicrons, X,Y);
            otherwise
                error('Unknown subregion: ''%s''.', subregions)
        end
        
        whichLevelsToContour = 1;
        fitEllipse = false;
        
        C = contourc(xAxis, yAxis, theRF, zLevels);
        faceAlpha = 0.2;
        edgeAlpha = 0.8;
        fillRFoutline(theAxes, C, zLevels, whichLevelsToContour, fitEllipse, faceAlpha, edgeAlpha);
%         if (strcmp(subregions, 'centers'))
%             displayConnectedConesPolygon(theAxes, centerIndices, conePositionsMicrons);
%         end
        
    end % mRGCindex
    
    set(theAxes, 'CLim', [0 1])
end

