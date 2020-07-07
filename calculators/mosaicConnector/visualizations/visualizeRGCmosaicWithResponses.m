function hFig = visualizeRGCmosaicWithResponses(figNo,theConeMosaic, xAxisScaling, ...
    xAxisData, theMidgetRGCmosaicResponses, ...
    xAxisDataFit, theMidgetRGCmosaicResponsesFit, ...
    eccentricityMicrons, sizeMicrons, ...
    theMidgetRGCmosaic,  zLevels, subregions, maxSpikeRate)

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

    % Zoomedin view, labeled cone types
    theAxes = axes('Position', [0.02 0.03 0.45 0.96]);
    renderPlot(theAxes, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
        theMidgetRGCmosaic, xAxis, yAxis, zLevels, subregions);
    xLims = eccentricityMicrons(1)+1.2*sizeMicrons(1)/2.0*[-1 1];
    yLims = eccentricityMicrons(2)+1.2*sizeMicrons(2)/2.0*[-1 1];
    axis(theAxes, 'equal');
    box(theAxes, 'on');
    set(theAxes, 'XLim', xLims, 'YLim', yLims, 'FontSize', 14);

    
    rgcsNum = size(theMidgetRGCmosaic.centerWeights,2);
    [~,RGCpositionsNormalized] = determineRGCPositionsFromCenterInputs(theConeMosaic, eccentricityMicrons, theMidgetRGCmosaic.centerWeights);
    
    w = 0.1/2;
    h = 0.1/2;
    gw = 0.5-w;
    gh = 0.96-h;
    m = 0.015/2;

    for iRGC = 1:rgcsNum
            ax = axes('Position', [0.48+gw*RGCpositionsNormalized(iRGC,1)+m 0.03+gh*RGCpositionsNormalized(iRGC,2)+m w h]);
            if (strcmp(xAxisScaling, 'log'))
                markerSize = 169;
                lineColor = [1 0 0];
            else
                markerSize = 100;
                lineColor = [0 1 1];
            end
            
            scatter(ax, xAxisData, squeeze(theMidgetRGCmosaicResponses(iRGC,:)), markerSize, '.', ...
                    'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1]);
            if (~isempty(xAxisDataFit))
                hold(ax, 'on')
                line(ax, xAxisDataFit, squeeze(theMidgetRGCmosaicResponsesFit(iRGC,:)), 'Color', lineColor, 'LineWidth', 1.5);
                hold(ax, 'off');
            end
            set(ax, 'XTickLabel', {}, 'YTickLabel', {}, ...
                'XLim', [xAxisData(1) xAxisData(end)], ...
                'XColor', 'none', 'YColor', 'none');
            set(ax, 'XScale', xAxisScaling);
            if (strcmp(xAxisScaling, 'log'))
                set(ax, 'XTick', [0.1 0.3 1 3 10 30], 'YTick', (0:0.25:1)*max(theMidgetRGCmosaicResponses(:)), ...
                    'YLim', [0 2*maxSpikeRate]);
            else
                set(ax, 'XTick',linspace(xAxisData(1), xAxisData(end), 5), 'YTick', max(abs(theMidgetRGCmosaicResponses(:)))*(-1:0.25:1), ...
                    'YLim', maxSpikeRate*[-1.0 1.0]);
            end
            box(ax, 'on'); grid(ax, 'on');
            axis(ax, 'square')
            drawnow;
    end
        
end

function renderPlot(theAxes, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
    theMidgetRGCmosaic,  xAxis, yAxis, zLevels, subregions)

    hold(theAxes, 'on');
    colormap(theAxes, brewermap(512, 'greys'));
    renderConeConnections(theAxes, theMidgetRGCmosaic, conePositionsMicrons);
    renderCones(theAxes, coneTypes, conePositionsMicrons, coneDiameterMicrons);
    renderRGCoutlines(theAxes,theMidgetRGCmosaic, conePositionsMicrons, coneSpacingsMicrons, subregions, zLevels,xAxis, yAxis);
   
    
    set(theAxes, 'CLim', [0 1]);
end

function renderRGCoutlines(theAxes, theMidgetRGCmosaic, conePositionsMicrons, coneSpacingsMicrons, subregions, zLevels, xAxis, yAxis)
    [X,Y] = meshgrid(xAxis,yAxis);
    rgcsNum = size(theMidgetRGCmosaic.centerWeights,2);
    
    for mRGCindex = 1:rgcsNum
        centerWeights = full(squeeze(theMidgetRGCmosaic.centerWeights(:, mRGCindex)));
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
        faceAlpha = 0.35;
        edgeAlpha = 0.5;
        fillRFoutline(theAxes, C, zLevels, whichLevelsToContour, fitEllipse, faceAlpha, edgeAlpha);
        
    end % mRGCindex
end

function renderConeConnections(theAxes, theMidgetRGCmosaic, conePositionsMicrons)
    rgcsNum = size(theMidgetRGCmosaic.centerWeights,2);
    for mRGCindex = 1:rgcsNum
        centerWeights = full(squeeze(theMidgetRGCmosaic.centerWeights(:, mRGCindex)));
        centerIndices = find(centerWeights>0);

        displayConnectedConesPolygon(theAxes, centerIndices, conePositionsMicrons);
    end % mRGCindex
end


function renderCones(theAxes, coneTypes, conePositionsMicrons, coneDiameterMicrons)

    global LCONE_ID
    global MCONE_ID
    global SCONE_ID

    xxx = cosd(0:10:360);
    yyy = sind(0:10:360);
    coneLabeling = true;
    
    if (coneLabeling)
        % Display cones
        LconeIndices = find(coneTypes == LCONE_ID);
        MconeIndices = find(coneTypes == MCONE_ID);
        SconeIndices = find(coneTypes == SCONE_ID);
    
        for k = 1:numel(LconeIndices)
            xx = conePositionsMicrons(LconeIndices(k),1) + 0.5*coneDiameterMicrons(LconeIndices(k))*xxx;
            yy = conePositionsMicrons(LconeIndices(k),2) + 0.5*coneDiameterMicrons(LconeIndices(k))*yyy;
            patch(theAxes, 'Faces', 1:numel(xx),'Vertices',[xx(:) yy(:)],'FaceColor',[1 0.6 0.6], 'FaceAlpha', 1, 'EdgeColor', [1 0 0], 'LineWidth', 1.5);
        end
        for k = 1:numel(MconeIndices)
            xx = conePositionsMicrons(MconeIndices(k),1) + 0.5*coneDiameterMicrons(MconeIndices(k))*xxx;
            yy = conePositionsMicrons(MconeIndices(k),2) + 0.5*coneDiameterMicrons(MconeIndices(k))*yyy;
            patch(theAxes, 'Faces', 1:numel(xx),'Vertices',[xx(:) yy(:)],'FaceColor',[0.7 1 0.7], 'FaceAlpha', 1, 'EdgeColor', [0 1 0], 'LineWidth', 1.5);
        end
        for k = 1:numel(SconeIndices)
            xx = conePositionsMicrons(SconeIndices(k),1) + 0.5*coneDiameterMicrons(SconeIndices(k))*xxx;
            yy = conePositionsMicrons(SconeIndices(k),2) + 0.5*coneDiameterMicrons(SconeIndices(k))*yyy;
            patch(theAxes, 'Faces', 1:numel(xx),'Vertices',[xx(:) yy(:)],'FaceColor',[0.7 0.7 1], 'FaceAlpha', 1, 'EdgeColor', [0 0 1], 'LineWidth', 1.5);
        end
    else
        for k = 1:size(conePositionsMicrons,1)
            xx = conePositionsMicrons(k,1) + 0.5*coneDiameterMicrons(k)*cosd(0:30:360);
            yy = conePositionsMicrons(k,2) + 0.5*coneDiameterMicrons(k)*sind(0:30:360);
            patch(theAxes, 'Faces', 1:numel(xx),'Vertices',[xx(:) yy(:)],'FaceColor',[0.9 0.9 0.9], 'FaceAlpha', 1, 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
        end
    end
    
end
