function hFig = visualizeRGCmosaicWithResponses(figNo,theConeMosaic, xAxisScaling, ...
    xAxisData, theMidgetRGCmosaicResponses, ...
    xAxisDataFit, theMidgetRGCmosaicResponsesFit, ...
    eccentricityMicrons, sizeMicrons, ...
    theMidgetRGCmosaic,  zLevels, subregions, maxSpikeRate, ...
    exportFig, condition, LMScontrast, targetRGC, labelCells)

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
    
    figWidthInches = 28; figHeightInches = 14;
    plotlabOBJ = setupPlotLab(0, figWidthInches, figHeightInches);
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1]);
    
    theAxes = axes('Position', [0.02*figWidthInches 0.02*figHeightInches 0.49*figWidthInches 0.98*figHeightInches]);
    renderMosaicConnectivityPlot(theAxes, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
        theMidgetRGCmosaic, xAxis, yAxis, zLevels, subregions, targetRGC);
    
    xLims = eccentricityMicrons(1)+1.2*sizeMicrons(1)/2.0*[-1 1];
    yLims = eccentricityMicrons(2)+1.2*sizeMicrons(2)/2.0*[-1 1];
    axis(theAxes, 'equal');
    box(theAxes, 'on');
    set(theAxes, 'XLim', xLims, 'YLim', yLims, 'FontSize', 14);

    if (~isempty(targetRGC))
        w0 = 0.55;
        h0 = 0.15;
        w = 0.35;
        h = 0.75;
        axesPosition = [...
                w0*figWidthInches, ...
                h0*figHeightInches, ...
                w*figWidthInches, ...
                h*figHeightInches];
            
        ax = axes('Position', axesPosition, 'Color', [1 1 1]);
          renderResponsePlot(ax, xAxisScaling, xAxisData, squeeze(theMidgetRGCmosaicResponses(targetRGC,:)), ...
                xAxisDataFit, squeeze(theMidgetRGCmosaicResponsesFit(targetRGC,:)), maxSpikeRate,  targetRGC, false);
    else
        rgcsNum = size(theMidgetRGCmosaic.centerWeights,2);
        [~,RGCpositionsNormalized] = determineRGCPositionsFromCenterInputs(theConeMosaic, eccentricityMicrons, theMidgetRGCmosaic.centerWeights);
        w = 0.1/2;
        h = 0.1/2;
        gw = 0.5-w;
        w0 = 0.50;
        gh = 0.83-h;
        h0 = 0.09;
        m = 0.015/2;
        for iRGC = 1:rgcsNum
            axesPosition = [...
                (w0+gw*RGCpositionsNormalized(iRGC,1)+m)*figWidthInches, ...
                (h0+gh*RGCpositionsNormalized(iRGC,2)+m)*figHeightInches, ...
                w*figWidthInches, ...
                h*figHeightInches];
            
            ax = axes('Position', axesPosition, 'Color', [1 1 1]);
            renderResponsePlot(ax, xAxisScaling, xAxisData, squeeze(theMidgetRGCmosaicResponses(iRGC,:)), ...
                xAxisDataFit, squeeze(theMidgetRGCmosaicResponsesFit(iRGC,:)), maxSpikeRate, iRGC, labelCells);
        end
    end
    
     if (exportFig)
          pdfFileName = sprintf('Ensemble_%s_LMS_%0.2f_%0.2f_%0.2f', condition,LMScontrast(1), LMScontrast(2), LMScontrast(3));
         plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, pwd());
    end
            
    setupPlotLab(-1);
end

function  renderResponsePlot(ax, xAxisScaling, xAxisData, yAxisData, xAxisDataFit, yAxisDataFit, maxSpikeRate, iRGC, labelCells)
    
    if (strcmp(xAxisScaling, 'log'))
        markerSize = 169;
        lineColor = [1 0 0];
    else
        markerSize = 100;
        lineColor = [0 1 1];
    end

    scatter(ax, xAxisData, yAxisData, markerSize, '.', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
    if (~isempty(xAxisDataFit))
        hold(ax, 'on')
        line(ax, xAxisDataFit, yAxisDataFit, 'Color', lineColor, 'LineWidth', 1.5);
        hold(ax, 'off');
    end
    
    set(ax, 'XTickLabel', {}, 'YTickLabel', {}, ...
        'XLim', [xAxisData(1) xAxisData(end)], ...
        'XColor', 'none', 'YColor', 'none');
    set(ax, 'XScale', xAxisScaling);
    set(ax, 'YTick', (0:0.25:1)*maxSpikeRate, 'YLim', [0 maxSpikeRate]);
    
    if (strcmp(xAxisScaling, 'log'))
        set(ax, 'XTick', [0.1 0.3 1 3 10 30])
    else
        set(ax, 'XTick',linspace(xAxisData(1), xAxisData(end), 5));
    end
    box(ax, 'on'); grid(ax, 'on');
    axis(ax, 'square')
    if (labelCells)
        if (strcmp(xAxisScaling, 'log'))
            xo = 0.12;
            yo = maxSpikeRate*0.85;
        else
            xo = xAxisData(1);
            yo = maxSpikeRate*0.85;
        end
        text(ax, xo,yo, sprintf('%d', iRGC), 'FontSize',10);
    end
    drawnow;
end

function renderMosaicConnectivityPlot(theAxes, conePositionsMicrons, coneDiameterMicrons, coneSpacingsMicrons, coneTypes, ...
    theMidgetRGCmosaic,  xAxis, yAxis, zLevels, subregions, targetRGC)

    hold(theAxes, 'on');
    colormap(theAxes, brewermap(512, 'greys'));
    renderConeConnections(theAxes, theMidgetRGCmosaic, conePositionsMicrons);
    renderCones(theAxes, coneTypes, conePositionsMicrons, coneDiameterMicrons);
    renderRGCoutlines(theAxes,targetRGC, coneTypes, theMidgetRGCmosaic, conePositionsMicrons, coneSpacingsMicrons, subregions, zLevels,xAxis, yAxis);  
    set(theAxes, 'CLim', [0 1]);
end

function renderRGCoutlines(theAxes, targetRGC, coneTypes, theMidgetRGCmosaic, conePositionsMicrons, coneSpacingsMicrons, subregions, zLevels, xAxis, yAxis)
    [X,Y] = meshgrid(xAxis,yAxis);
    rgcsNum = size(theMidgetRGCmosaic.centerWeights,2);
    
    for mRGCindex = 1:rgcsNum
        
        if ((~isempty(targetRGC)) && (mRGCindex ~= targetRGC))
            continue;
        end
        centerWeights = full(squeeze(theMidgetRGCmosaic.centerWeights(:, mRGCindex)));
        %centerWeights(centerWeights>0) = 1;
        
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

function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'in', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 14, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 

