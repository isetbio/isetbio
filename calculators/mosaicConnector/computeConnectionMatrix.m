function [connectionMatrix, conePositionsMicrons, RGCRFPositionsMicrons, coneSpacingsMicrons] = computeConnectionMatrix(RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios, roi, thresholdSeparationMicronsForRemovingUnitsFromMosaic)

    % Find cones within the roi
    idxCones = positionsWithinROI(roi, conePositionsMicrons,  thresholdSeparationMicronsForRemovingUnitsFromMosaic);
    conePositionsMicrons = conePositionsMicrons(idxCones,:);
    coneSpacingsMicrons = coneStats(conePositionsMicrons);
      
    % Find RGCs within the roi
    idxRGC = positionsWithinROI(roi, RGCRFPositionsMicrons,  thresholdSeparationMicronsForRemovingUnitsFromMosaic );
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(idxRGC,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(idxRGC);
    desiredConesToRGCratios = desiredConesToRGCratios(idxRGC);
    
    % Step1. Align each RGC with its nearest cone. This ensure all RGC's
    % are connected to at least one cone. Since cones are more numerous
    % than RGCs some cones will not connect to an RGC at this step. This 
    % step occurs only for RGCs for which the cone-to-RGC ratio is [1..2]
    
    visualizeProcess = true;
    RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
        conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        desiredConesToRGCratios, visualizeProcess);
    
    % Step 2. Connect each cone to its neigboring RGC. Connection weights
    % depend on 2 factors: cone-to-RGC ration and proximity
    connectionMatrix = connectConesToRGC(conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        desiredConesToRGCratios, visualizeProcess);
end

    
function indices = positionsWithinROI(roi, positions, thresholdSeparation)
    d = bsxfun(@minus,positions, roi.center);
    ecc = sqrt(sum(positions.^2,2));
    indices = find((abs(d(:,1)) <= 0.5*roi.size(1)) & (abs(d(:,2)) <= 0.5*roi.size(2)));
    
    % Re-order according to increasing eccentricity
    [~,sortedIdx] = sort(ecc(indices), 'ascend');
    indices = indices(sortedIdx);
    
    % mosaic correction: Remove units that are too close to each other
    % Compute all distances in included positions
    posIncluded = positions(indices,:);
    indicesToBeRemoved = indicesOfUnitsWithNearestSeparationLessThanThresholdSeparation(posIncluded, thresholdSeparation);
    % Remove units that are too close to each other
    indices = setdiff(indices, indices(indicesToBeRemoved));
end

function indices = indicesOfUnitsWithNearestSeparationLessThanThresholdSeparation(positions, thresholdSeparation)
    % Find the closest distance from each cone to all other cones
    allPairwiseDistances = pdist(positions);
    %thresholdSeparation = prctile(allPairwiseDistances(:), 0.02);
    idx = find(allPairwiseDistances < thresholdSeparation);
    [unit1Indices,unit2Indices] = returnRowColFromLowerTriIndex(idx, size(positions,1));
    indices = unit1Indices;
    for k = 1:numel(indices)
        distEst1 = sqrt(sum((positions(unit1Indices(k),:)-positions(unit2Indices(k),:)).^2));
        distEst2 = allPairwiseDistances(indices(k));
        fprintf('Elements %d and %d are too close. (distance = %2.2f, %2.2f). Removed %d.\n', ...
            unit1Indices(k), unit2Indices(k), distEst1, distEst2, unit1Indices(k));
    end
end

function [rows,cols] = returnRowColFromLowerTriIndex(indices, N)
    totalElementsNum = (N*(N-1))/2;
    reversedIndices = totalElementsNum-indices;
    k = floor( (sqrt(1+8*reversedIndices)-1)/2);
    j = reversedIndices - k.*(k+1)/2;
    rows = N-j;
    cols = N-(k+1);
end
function coneSpacings = coneStats(conePositions)
    p = pdist2(conePositions, conePositions, 'euclidean', 'Smallest', 3);
    p = p(2:end,:);
    coneSpacings = mean(p,1);
end