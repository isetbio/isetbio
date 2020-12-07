function hFig = visualizeConeAndRGCmosaicsWithRetinalImage(theConeMosaic, eccentricityMicrons, sizeMicrons, ...
    theMidgetRGCmosaic,  zLevels, subregions, theOI)

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
    
    
    hFig = figure(99); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 1, ...
            'colsNum', 2, ...
            'leftMargin', 0.05, ...
            'bottomMargin', 0.05, ...
            'rightMargin', 0.03, ...
            'topMargin', 0.1);
        
    % Full view, no cone labeling
    renderPlot(theAxesGrid{1,1}, theOI, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
        theMidgetRGCmosaic, eccentricityMicrons, xAxis, yAxis, zLevels, subregions, false);
    xLims = [xAxis(1) xAxis(end)];
    yLims = [yAxis(1) yAxis(end)];
    set(theAxesGrid{1,1}, 'XLim', xLims, 'YLim', yLims);
    
    % Zoomedin view, labeled cone types
    renderPlot(theAxesGrid{1,2}, theOI, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
        theMidgetRGCmosaic, eccentricityMicrons, xAxis, yAxis, zLevels, subregions, true);
    xLims = eccentricityMicrons(1)-10+sizeMicrons(1)/2.0*[-1 1];
    yLims = eccentricityMicrons(2)-10+sizeMicrons(2)/2.0*[-1 1];
    set(theAxesGrid{1,2}, 'XLim', xLims, 'YLim', yLims);
end

function renderPlot(theAxes, theOI, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
    theMidgetRGCmosaic, eccentricityMicrons, xAxis, yAxis, zLevels, subregions, coneLabeling)

    global LCONE_ID
    global MCONE_ID
    global SCONE_ID

    [X,Y] = meshgrid(xAxis,yAxis);
    
    % Render the optical image
    spatialSupportMM = oiGet(theOI, 'spatial support', 'mm');
    spatialSupportMicrons = 1e3 * squeeze(spatialSupportMM(1,:,1));
    imagesc(theAxes, ...
            spatialSupportMicrons+eccentricityMicrons(1), ...
            spatialSupportMicrons+eccentricityMicrons(2), ...
            oiGet(theOI, 'rgbimage'));
    axis(theAxes, 'xy')
    axis(theAxes, 'image')
    hold(theAxes, 'on');
    colormap(theAxes, brewermap(512, 'greys'));
    
    if (coneLabeling)
        % Display cones
        LconeIndices = find(coneTypes == LCONE_ID);
        MconeIndices = find(coneTypes == MCONE_ID);
        SconeIndices = find(coneTypes == SCONE_ID);
    
        for k = 1:numel(LconeIndices)
            xx = conePositionsMicrons(LconeIndices(k),1) + 0.5*coneDiameterMicrons(LconeIndices(k))*cosd(0:30:360);
            yy = conePositionsMicrons(LconeIndices(k),2) + 0.5*coneDiameterMicrons(LconeIndices(k))*sind(0:30:360);
            line(theAxes, xx, yy, 'Color', 'r');
        end
        for k = 1:numel(MconeIndices)
            xx = conePositionsMicrons(MconeIndices(k),1) + 0.5*coneDiameterMicrons(MconeIndices(k))*cosd(0:30:360);
            yy = conePositionsMicrons(MconeIndices(k),2) + 0.5*coneDiameterMicrons(MconeIndices(k))*sind(0:30:360);
            line(theAxes, xx, yy, 'Color', [0 0.8 0]);
        end
        for k = 1:numel(SconeIndices)
            xx = conePositionsMicrons(SconeIndices(k),1) + 0.5*coneDiameterMicrons(SconeIndices(k))*cosd(0:30:360);
            yy = conePositionsMicrons(SconeIndices(k),2) + 0.5*coneDiameterMicrons(SconeIndices(k))*sind(0:30:360);
            line(theAxes, xx, yy, 'Color', 'b');
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

