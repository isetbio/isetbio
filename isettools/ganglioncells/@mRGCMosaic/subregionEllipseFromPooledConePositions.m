function [contourData, theCenter, xAlpha, yAlpha, theRotationRadians] = subregionEllipseFromPooledConePositions(...
                         theConePositions, ellipseContourPoints, ellipseContourAngles, pLevel, maxNumberOfConesOutsideContour)

    dataDimensionality = size(theConePositions,2);
    conesNum = size(theConePositions,1);
    unitCirclePoints = [cosd(ellipseContourAngles(:)) sind(ellipseContourAngles(:))];
    theCenter = mean(theConePositions,1);

    if (~isempty(pLevel))
        pLevel = min([1 max([0 pLevel])]);
        k = finv(pLevel,dataDimensionality,conesNum-dataDimensionality) * dataDimensionality * (conesNum-1)/(conesNum-dataDimensionality);
        [principalComponentVectors,scores,principalComponentVariances] = pca(theConePositions);
        ab = diag(sqrt(k * principalComponentVariances));
    else
        startingPLevel = 0.6;
        k = finv(startingPLevel,dataDimensionality,conesNum-dataDimensionality) * dataDimensionality * (conesNum-1)/(conesNum-dataDimensionality);
        [principalComponentVectors,scores,principalComponentVariances] = pca(theConePositions);
        ab = diag(sqrt(k * principalComponentVariances));

        numberOfConesOutsideContour = conesNum;
        coneXcoords = theConePositions(:,1);
        coneYcoords = theConePositions(:,2);
        while (numberOfConesOutsideContour > maxNumberOfConesOutsideContour)
            contourData.vertices = bsxfun(@plus, unitCirclePoints * ab * principalComponentVectors', theCenter);
             
            indicesOfConePositionsInsideContour = inpolygon(coneXcoords, coneYcoords, contourData.vertices(:,1), contourData.vertices(:,2));
            numberOfConesOutsideContour = numel(coneXcoords(~indicesOfConePositionsInsideContour));
            % Increase alphas so we can include more cones within the contour
            ab = ab * 1.01;
        end
    end

    contourData.vertices = bsxfun(@plus, unitCirclePoints * ab * principalComponentVectors', theCenter);
    contourData.faces = 1:size(contourData.vertices,1);


    xAlpha = ab(1,1); yAlpha = ab(2,2);
    theRotationRadians = atan2(principalComponentVectors(2,1), principalComponentVectors(1,1));
end
