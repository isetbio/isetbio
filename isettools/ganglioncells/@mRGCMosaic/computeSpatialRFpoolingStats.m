function [centerRadiusDegsX , centerRadiusDegsY, surroundRadiusDegs, ...
    centerRadiusMicronsX, centerRadiusMicronsY, surroundRadiusMicrons, ...
    eccDegs, eccMicrons] = computeSpatialRFpoolingStats(obj, theRGCindex, ...
    minConeWeightForRFcenterPooling, ...
    minConeWeightForRFsurroundPooling)

    connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
    pooledConeIndicesAndWeights.centerConeIndices = find(connectivityVector > 0);
    pooledConeIndicesAndWeights.centerConeWeights = connectivityVector(pooledConeIndicesAndWeights.centerConeIndices);
    
    connectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex)));
    pooledConeIndicesAndWeights.surroundConeIndices = find(connectivityVector > 0);
    pooledConeIndicesAndWeights.surroundConeWeights = connectivityVector(pooledConeIndicesAndWeights.surroundConeIndices);
    

    maxValue = max(pooledConeIndicesAndWeights.centerConeWeights);

    idx = find(pooledConeIndicesAndWeights.centerConeWeights >=  maxValue * minConeWeightForRFcenterPooling);
    identifiedCenterConeIndices = pooledConeIndicesAndWeights.centerConeIndices(idx);
    identifiedCenterConeWeights = pooledConeIndicesAndWeights.centerConeWeights(idx);

    idx = find(pooledConeIndicesAndWeights.surroundConeWeights >= maxValue * minConeWeightForRFsurroundPooling);
    identifiedSurroundConeIndices = pooledConeIndicesAndWeights.surroundConeIndices(idx); 
    identifiedSurroundConeWeights = pooledConeIndicesAndWeights.surroundConeWeights(idx); 

    centerConePositionsDegs = obj.inputConeMosaic.coneRFpositionsDegs(identifiedCenterConeIndices,:);
    surroundConePositionsDegs = obj.inputConeMosaic.coneRFpositionsDegs(identifiedSurroundConeIndices,:);

    centerConePositionsMicrons = obj.inputConeMosaic.coneRFpositionsMicrons(identifiedCenterConeIndices,:);
    surroundConePositionsMicrons = obj.inputConeMosaic.coneRFpositionsMicrons(identifiedSurroundConeIndices,:);

    centerCenterMicrons = mean(centerConePositionsMicrons,1);
    centerCenterDegs = mean(centerConePositionsDegs,1);
    eccMicrons = sqrt(sum(centerCenterMicrons.^2,2));
    eccDegs = sqrt(sum(centerCenterDegs.^2,2));


    if (numel(centerConePositionsDegs) == 0)
        fprintf(2,'Zero center cones. Returning nan.\n')
        centerRadiusDegsX = nan;
        centerRadiusDegsY = nan;
        surroundRadiusDegs = nan;
        centerRadiusMicronsX = nan;
        centerRadiusMicronsY = nan;
        surroundRadiusMicrons = nan;
        eccMicrons = nan;
        eccDegs = nan;
        return;
    end
 

    ellipseContourAngles = 0:10:350;
    maxNumberOfConesOutsideContour = 1;
    [~, ~, centerRadiusDegsX , centerRadiusDegsY, theRotationRadians] = ...
        mRGCMosaic.subregionEllipseFromPooledConePositions(...
             centerConePositionsDegs, [], ellipseContourAngles, [], maxNumberOfConesOutsideContour);

     [~, ~, centerRadiusMicronsX , centerRadiusMicronsY, theRotationRadians] = ...
        mRGCMosaic.subregionEllipseFromPooledConePositions(...
             centerConePositionsMicrons, [], ellipseContourAngles, [], maxNumberOfConesOutsideContour);


    if (numel(surroundConePositionsDegs) == 0)
        fprintf(2,'Zero surround cones. Returning nan.\n');
        [max(pooledConeIndicesAndWeights.surroundConeWeights) maxValue * minConeWeightForRFsurroundPooling]
        surroundRadiusDegs = nan;
        surroundRadiusMicrons = nan;
        return;
    end
 
    surroundCenterDegs = mean(surroundConePositionsDegs,1);
    surroundRadiusDegs = max(sqrt(sum((bsxfun(@minus, surroundConePositionsDegs, surroundCenterDegs)).^2,2)));

    surroundCenterMicrons = mean(surroundConePositionsMicrons,1);
    surroundRadiusMicrons = max(sqrt(sum((bsxfun(@minus, surroundConePositionsMicrons, surroundCenterMicrons)).^2,2)));


    debug = false;

    if (debug)
        centerCenter = mean(centerConePositions,1);
        centerCircleXX = contourData.vertices(:,1);
        centerCircleYY = contourData.vertices(:,2);
    
        surroudCircleXX = surroundCenter(1) + cosd(0:10:360) * surroundRadiusDegs;
        surroudCircleYY = surroundCenter(2) + sind(0:10:360) * surroundRadiusDegs;
    
        xx = prctile(surroundConePositions(:,1), [0 50 100]);
        yy = prctile(surroundConePositions(:,2), [0 50 100]);
    
        xRange = xx(3)-xx(1);
        yRange = yy(3)-yy(1);
        maxRange = max([xRange yRange]);
        xLims = xx(2) + 0.5*maxRange*[-1 1];
        yLims = yy(2) + 0.5*maxRange*[-1 1];
    
        figure(111); clf;
        subplot(1,2,1);
        plot(centerConePositions(:,1), centerConePositions(:,2), 'ko');
        hold on;
        plot(centerCircleXX, centerCircleYY, 'r-');
        plot(centerCenter(1) + centerRadiusDegsX *[-1 1], centerCenter(2)*[1 1], 'r-');
        plot(centerCenter(1)*[1 1], centerCenter(2) + centerRadiusDegsX  *[-1 1], 'b-');
        title('center cones');
        axis 'equal'
        set(gca, 'XLim', xLims, 'YLim', yLims);

        subplot(1,2,2);
        plot(surroundConePositions(:,1), surroundConePositions(:,2), 'ko');
        hold on;
        plot(surroudCircleXX, surroudCircleYY, 'b-');
        title('surround cones');
        axis 'equal'
        set(gca, 'XLim', xLims, 'YLim', yLims);
    end

end





