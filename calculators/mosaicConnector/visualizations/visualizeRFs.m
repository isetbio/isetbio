function visualizeRFs(patchEccDegs, zLevels, whichLevelsToContour, connectivityMatrix, RGCRFPositionsMicrons, ...
    conePositionsMicrons, coneSpacingsMicrons, coneTypes, roi, fitEllipse,  showConnectedConePolygon, plotlabOBJ, exportFileName, exportsDir)

    % Define constants
    global LCONE_ID
    global MCONE_ID
    global SCONE_ID 
    
    
    % Sampling for contours
    deltaX = 0.2;
    xAxis = (roi.center(1)-roi.size(1)/2): deltaX: (roi.center(1)+roi.size(1)/2);
    yAxis = (roi.center(2)-roi.size(2)/2): deltaX: (roi.center(2)+roi.size(2)/2);
    [X,Y] = meshgrid(xAxis,yAxis);
    
    hFig = figure(99); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.05, ...
            'bottomMargin', 0.05, ...
            'rightMargin', 0.03, ...
            'topMargin', 0.1);
        
    theAxesGrid = theAxesGrid{1,1};
    set(theAxesGrid, 'XLim', roi.center(1)+roi.size(1)/2*[-1 1], ...
        'YLim', roi.center(2)+roi.size(2)/2*[-1 1]);
    hold(theAxesGrid, 'on');
    
    rgcsNum = size(RGCRFPositionsMicrons,1);
    mixedInputRGCs = 0;
    
    coneInputsPerRGC = zeros(1,100);
    rgcsNumWithNonZeroInputs = 0;
    
    for RGCindex = 1:rgcsNum
        
        connectivityVector = full(squeeze(connectivityMatrix(:, RGCindex)));
        inputIDs = find(connectivityVector > 0.01);
        inputsNum = numel(inputIDs);
        
        if (inputsNum == 0)
            fprintf(2, 'RGC %d has zero inputs!!!\n', RGCindex);
            continue;
        end
        
        rgcsNumWithNonZeroInputs = rgcsNumWithNonZeroInputs + 1;
        coneInputsPerRGC(inputsNum) = coneInputsPerRGC(inputsNum) + 1;
  
        inputTypesNum = numel(unique(coneTypes(inputIDs)));
        if (inputTypesNum>1)
            mixedInputRGCs = mixedInputRGCs + 1;
        end
        
        % Generate RF centers of RGCs based on cone positions and connection matrix
        theRF = generateRGCRFcenterSubregionFromConnectivityMatrix(...
            connectivityVector/max(connectivityVector), conePositionsMicrons, coneSpacingsMicrons, X,Y);

        if (isempty(theRF))
            fprintf(2,'No cone inputs to this RF -> No visualization\n');
            continue;
        end
        
        C = contourc(xAxis, yAxis,theRF, zLevels);
        faceAlpha = 0.4;
        edgeAlpha = 0.8;
        fillRFoutline(theAxesGrid, C, zLevels, whichLevelsToContour, fitEllipse, faceAlpha, edgeAlpha);
    
        indicesOfConeInputsToThisRGC = find(connectivityVector>0);
        
        if (showConnectedConePolygon)
            % Connected cones
            displayConnectedConesPolygon(theAxesGrid, indicesOfConeInputsToThisRGC, conePositionsMicrons);
        end
        
        if (mod(RGCindex-1,10) == 9)
            drawnow;
        end
    end % RGCindex
       
    % Only show cones within ROI
    xPos = conePositionsMicrons(:,1);
    yPos = conePositionsMicrons(:,2);
    indicesOfConnectedCones = find((xPos>xAxis(1)) & (xPos<xAxis(end)) & (yPos>yAxis(1)) & (yPos<yAxis(end)));
    conePositionsMicrons = conePositionsMicrons(indicesOfConnectedCones,:);
    coneTypes = coneTypes(indicesOfConnectedCones);
    

    % Display cones
    LconeIndices = find(coneTypes == LCONE_ID);
    MconeIndices = find(coneTypes == MCONE_ID);
    SconeIndices = find(coneTypes == SCONE_ID);
    scatter(theAxesGrid,conePositionsMicrons(LconeIndices,1), conePositionsMicrons(LconeIndices,2), 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]);
    scatter(theAxesGrid,conePositionsMicrons(MconeIndices,1), conePositionsMicrons(MconeIndices,2), 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0.5 0.9 0.5]);
    scatter(theAxesGrid,conePositionsMicrons(SconeIndices,1), conePositionsMicrons(SconeIndices,2), 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0.5 0.5 1.0]);
    
    
    conesNum = numel(LconeIndices)+numel(MconeIndices);
    
    ratio = conesNum/rgcsNum;
    title(theAxesGrid,...
        sprintf('\\color[rgb]{0.3 0.3 0.3} eccentricity = %2.1f degs, \\color[rgb]{1.0 0.1 0.1} LMcones-to-mRGCs ratio = %2.2f\n\\color[rgb]{0.4 0.4 0.6}RGCs with 1 input: %2.1f%%, 2 inputs: %2.1f%%, 3 inputs: %2.1f%%, 4+ inputs: %2.1f%%\n RGCs with mixed cone inputs: %2.1f%%', ...
        patchEccDegs(1), ratio, ...
        100*coneInputsPerRGC(1)/rgcsNumWithNonZeroInputs, ...
        100*coneInputsPerRGC(2)/rgcsNumWithNonZeroInputs, ...
        100*coneInputsPerRGC(3)/rgcsNumWithNonZeroInputs, ...
        100*sum(coneInputsPerRGC(4:end))/rgcsNumWithNonZeroInputs, ...
        100*mixedInputRGCs/rgcsNumWithNonZeroInputs));
     
    xLims = [xAxis(1) xAxis(end)];
    yLims = [yAxis(1) yAxis(end)];
    ylabel(theAxesGrid, 'microns');
    set(theAxesGrid, 'CLim', [0 1], 'XLim', xLims, 'YLim', yLims);
    
    colormap(brewermap(512, 'greys'));
 
    
    % Export the figure to the gallery directory in PNG format
    exportFig = true;
    if (exportFig)
        rfEccDegs = WatsonRGCModel.rhoMMsToDegs(roi.center/1000);
        fName = sprintf('%s__RFecc_x=%2.2f_y=%2.2fdegs', exportFileName, rfEccDegs(1), rfEccDegs(2));
        plotlabOBJ.exportFig(hFig, 'png', fName, exportsDir);
    end
end

