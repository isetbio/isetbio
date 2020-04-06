function visualizeRFs(connectivityMatrix, conePositionsMicrons, RGCRFPositionsMicrons, coneSpacingsMicrons, roi)

    % Sampling for contours
    deltaX = 0.2;
    xAxis = (roi.center(1)-roi.size(1)/2): deltaX: (roi.center(1)+roi.size(1)/2);
    yAxis = (roi.center(2)-roi.size(2)/2): deltaX: (roi.center(2)+roi.size(2)/2);
    [X,Y] = meshgrid(xAxis,yAxis);
    
   
    zLevels = 0.2:0.2:0.7;
    
    whichLevelsToContour = 1:numel(zLevels);
    
    hFig = figure(1); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
        'leftMargin', 0.05, 'bottomMargin', 0.15, 'topMargin', 0.05);
    theAxesGrid = theAxesGrid{1,1};
    set(theAxesGrid, 'XLim', roi.center(1)+roi.size(1)/2*[-1 1], 'YLim', roi.center(2)+roi.size(2)/2*[-1 1]);
    hold(theAxesGrid, 'on');
    
    rgcsNum = size(RGCRFPositionsMicrons,1);
    for RGCindex = 1:rgcsNum
        
        % Generate RFs of RGCs based on cone psotions and connection matrix
        theRF = generateRGCRFsFromConnectivityMatrix(...
            squeeze(connectivityMatrix(:, RGCindex)), conePositionsMicrons, coneSpacingsMicrons, X,Y);
    
        if (all(isnan(theRF(:))))
            fprintf('Zero RF. skipping\n');
            continue;
        end
        %contourf(theAxesGrid, xAxis, yAxis, theRF, zLevels);
        C = contourc(xAxis, yAxis,theRF, zLevels);
        renderContourPlot(theAxesGrid, C, zLevels, whichLevelsToContour );
        
    end
        
    scatter(conePositionsMicrons(:,1), conePositionsMicrons(:,2), 'r');
    
    colormap(brewermap(512, 'greys'))
    set(theAxesGrid, 'CLim', [0 1]);
end

function renderContourPlot(theAxes, C, zLevels, whichLevelsToContour )
    k = 1;
    contoursNum = 0;
    while k < size(C,2)
        level = C(1,k);
        points = C(2,k);
        if (level == zLevels(whichLevelsToContour(1)))
            edgeAlpha = 0.5;
            faceAlpha = 0.2;
        else
            edgeAlpha = 0;
            faceAlpha = 0.2+(level)*0.05;
        end
        
        xRGCEnsembleOutline = C(1,k+(1:points));
        yRGCEnsembleOutline = C(2,k+(1:points));
        v = [xRGCEnsembleOutline(:) yRGCEnsembleOutline(:)];
        f = 1:numel(xRGCEnsembleOutline);
        patch(theAxes, 'Faces', f, 'Vertices', v, 'FaceColor', [0.7 0.7 0.7]-level*0.05, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', [0 0 0], 'EdgeAlpha', edgeAlpha, 'LineWidth', 1.0);
        
        k = k+points+1;
        contoursNum = contoursNum + 1;
    end


end


function theRF = generateRGCRFsFromConnectivityMatrix(connectivityVectorForRGC, conePositions, coneSpacings,  X,Y)
    
    theRF = [];
    connectedConeIDs = find(connectivityVectorForRGC>0);
    
    for k = 1:numel(connectedConeIDs)
        coneIndex = connectedConeIDs(k);
        cP = squeeze(conePositions(coneIndex,:));
        coneSigma = coneSpacings(coneIndex)/6;
        coneProfile = exp(-0.5*((X-cP(1))/coneSigma).^2) .* exp(-0.5*((Y-cP(2))/coneSigma).^2);
        coneProfile = coneProfile.^0.25;
        if (isempty(theRF))
            theRF = coneProfile * connectivityVectorForRGC(coneIndex);
        else
            theRF = theRF + coneProfile  * connectivityVectorForRGC(coneIndex);
        end
    end
    theRF(theRF>0.6) = 1;
    
end




