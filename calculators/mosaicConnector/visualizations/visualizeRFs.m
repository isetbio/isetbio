function visualizeRFs(patchEccDegs, zLevels, whichLevelsToContour, connectivityMatrix, RGCRFPositionsMicrons, ...
    conePositionsMicrons, coneSpacingsMicrons, coneTypes, roi, fitEllipse,  plotlabOBJ, exportFileName, exportsDir)

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
            'topMargin', 0.07);
        
    theAxesGrid = theAxesGrid{1,1};
    set(theAxesGrid, 'XLim', roi.center(1)+roi.size(1)/2*[-1 1], ...
        'YLim', roi.center(2)+roi.size(2)/2*[-1 1]);
    hold(theAxesGrid, 'on');
    
    rgcsNum = size(RGCRFPositionsMicrons,1);
    multiInputRGCs = 0;
    mixedInputRGCs = 0;
    
    for RGCindex = 1:rgcsNum
        
        connectivityVector = full(squeeze(connectivityMatrix(:, RGCindex)));
        inputIDs = find(connectivityVector>0);
        inputsNum = numel(inputIDs);
        if (inputsNum > 1)
            multiInputRGCs = multiInputRGCs + 1;
            inputTypesNum = numel(unique(coneTypes(inputIDs)));
            if (inputTypesNum>1)
                mixedInputRGCs = mixedInputRGCs + 1;
            end
        end
        
        % Generate RF centers of RGCs based on cone positions and connection matrix
        theRF = generateRGCRFcenterSubregionFromConnectivityMatrix(...
            connectivityVector, conePositionsMicrons, coneSpacingsMicrons, X,Y);

        C = contourc(xAxis, yAxis,theRF, zLevels);
        renderContourPlot(theAxesGrid, C, zLevels, whichLevelsToContour, fitEllipse);
        
        indicesOfConeInputsToThisRGC = find(connectivityVector>0);
        
        showConnectedConePolygon = true;
        if (showConnectedConePolygon)
            % Connected cones
            displayConnectedConesPolygon(indicesOfConeInputsToThisRGC, conePositionsMicrons);
        end
    end
       
    % Only show cones within ROI
    xPos = conePositionsMicrons(:,1);
    yPos = conePositionsMicrons(:,2);
    indicesOfConnectedCones = find((xPos>xAxis(1)) & (xPos<xAxis(end)) & (yPos>yAxis(1)) & (yPos<yAxis(end)));
    conePositionsMicrons = conePositionsMicrons(indicesOfConnectedCones,:);
    coneTypes = coneTypes(indicesOfConnectedCones);
    

    % Display cones
    LconeIndices = find(coneTypes == 2);
    MconeIndices = find(coneTypes == 3);
    SconeIndices = find(coneTypes == 4);
    scatter(conePositionsMicrons(LconeIndices,1), conePositionsMicrons(LconeIndices,2), 'r');
    scatter(conePositionsMicrons(MconeIndices,1), conePositionsMicrons(MconeIndices,2), 'g');
    scatter(conePositionsMicrons(SconeIndices,1), conePositionsMicrons(SconeIndices,2), 'b');
    
    
    conesNum = numel(LconeIndices)+numel(MconeIndices);
    
    ratio = conesNum/rgcsNum;
    title(theAxesGrid,sprintf('\\color[rgb]{0.3 0.3 0.3} eccentricity = %2.1f degs, \\color[rgb]{1.0 0.1 0.1} LMcones-to-mRGCs ratio = %2.2f\n\\color[rgb]{0.4 0.4 0.6}RGCs with > 1 inputs = %d/%d (%2.1f%%), RGCs with mixed inputs = %d/%d (%2.1f%%)', patchEccDegs(1), ratio, multiInputRGCs,rgcsNum, 100*multiInputRGCs/rgcsNum, mixedInputRGCs, rgcsNum, 100*mixedInputRGCs/rgcsNum));
     
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

function  [semiAxes, rfCenter] = renderContourPlot(theAxes, C, zLevels, whichLevelsToContour, fitEllipse )
    k = 1;
    contoursNum = 0;
    while k < size(C,2)
        level = C(1,k);
        points = C(2,k);
        if (level == zLevels(whichLevelsToContour(1)))
            edgeAlpha = 0.7;
            faceAlpha = 0.2;
        elseif (ismember(level, zLevels(whichLevelsToContour)))
            edgeAlpha = 0;
            faceAlpha = 0.2+(level)*0.05;
        else
            % skip this contour
            k = k+points+1;
            continue;
        end
        
        xRGCEnsembleOutline = C(1,k+(1:points));
        yRGCEnsembleOutline = C(2,k+(1:points));
       
        if (fitEllipse)
            [xRGCEnsembleOutline,  yRGCEnsembleOutline, ...
                semiAxes, rfCenter, noFit] = fitEllipseToContour(xRGCEnsembleOutline,  yRGCEnsembleOutline);
        else
            semiAxes = [nan nan];
            rfCenter = [nan nan];
        end
        
        faceColor = [0.5 0.5 0.5]-level*0.05;
        edgeColor = [0.2 0.2 0.2];
        patchContour(theAxes, xRGCEnsembleOutline, yRGCEnsembleOutline, faceColor, edgeColor, faceAlpha, edgeAlpha);

        k = k+points+1;
        contoursNum = contoursNum + 1;
    end
end

function patchContour(theAxes, xRGCEnsembleOutline, yRGCEnsembleOutline, faceColor, edgeColor, faceAlpha, edgeAlpha)
    v = [xRGCEnsembleOutline(:) yRGCEnsembleOutline(:)];
    f = 1:numel(xRGCEnsembleOutline);
    patch(theAxes, 'Faces', f, 'Vertices', v, 'FaceColor', faceColor, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', edgeColor, ... 
           'EdgeAlpha', edgeAlpha, 'LineWidth', 1.0);
end

function displayConnectedConesPolygon(indicesOfConeInputs, conePositionsMicrons)
    % Polygon connecting input cones
    xx = conePositionsMicrons(indicesOfConeInputs,1);
    yy = conePositionsMicrons(indicesOfConeInputs,2);
    xo = mean(xx);
    yo = mean(yy);
    dx = xx-xo;
    dy = yy-yo;
    [~,idx] = sort(unwrap(atan2(dy,dx)));
    xx = xx(idx);
    yy = yy(idx);
    xx(end+1) = xx(1);
    yy(end+1) = yy(1);
    plot(xx,yy, 'k:', 'LineWidth', 1.0);
end
