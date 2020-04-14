function [semiAxes,rfCenters] = computeRFsizes(zLevels, whichLevelsToContour, connectivityMatrix, conePositionsMicrons, RGCRFPositionsMicrons, coneSpacingsMicrons, coneTypes, roi)

    % Sampling for contours
    deltaX = 0.2;
    xAxis = (roi.center(1)-roi.size(1)/2): deltaX: (roi.center(1)+roi.size(1)/2);
    yAxis = (roi.center(2)-roi.size(2)/2): deltaX: (roi.center(2)+roi.size(2)/2);
    [X,Y] = meshgrid(xAxis,yAxis);
    
    rgcsNum = size(RGCRFPositionsMicrons,1);    
    semiAxes = zeros(rgcsNum, 2);
    rfCenters = zeros(rgcsNum,2);

    parfor RGCindex = 1:rgcsNum
        fprintf('Fitting cell %d of %d\n', RGCindex, rgcsNum);
        connectivityVector = full(squeeze(connectivityMatrix(:, RGCindex)));
        
        % Generate RFs of RGCs based on cone positions and connection matrix
        theRF = generateRGCRFsFromConnectivityMatrix(...
            connectivityVector, conePositionsMicrons, coneSpacingsMicrons, X,Y);

        C = contourc(xAxis, yAxis,theRF, zLevels);
        [semiAxes(RGCindex,:), rfCenters(RGCindex,:)] = computeContourPlot(C, zLevels, whichLevelsToContour);
        
    end
end

function  [semiAxes, rfCenter] = computeContourPlot(C, zLevels, whichLevelsToContour)
    k = 1;
    contoursNum = 0;
    while k < size(C,2)
        level = C(1,k);
        points = C(2,k);
        if (level == zLevels(whichLevelsToContour(1)))
        elseif (ismember(level, zLevels(whichLevelsToContour)))
        else
            % skip this contour
            k = k+points+1;
            continue;
        end
        
        xRGCEnsembleOutline = C(1,k+(1:points));
        yRGCEnsembleOutline = C(2,k+(1:points));
        
        [~,~,  semiAxes, rfCenter, noFit] = fitEllipseToContour(xRGCEnsembleOutline,  yRGCEnsembleOutline);
        
        k = k+points+1;
        contoursNum = contoursNum + 1;
    end
end



