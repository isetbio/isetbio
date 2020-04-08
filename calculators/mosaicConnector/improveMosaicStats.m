function [conePositionsMicrons, coneSpacingsMicrons,...
          RGCRFPositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios] = ...
    improveMosaicStats(conePositionsMicrons, ...
                       RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
                       desiredConesToRGCratios, thresholdFraction, roi)

    % Find cones within the roi
    idxCones = positionsWithinROI(roi, conePositionsMicrons);
    conePositionsMicrons = conePositionsMicrons(idxCones,:);
    coneSpacingsMicrons = coneStats(conePositionsMicrons);
    
    % Remove cones that may be too close to each other
    idxCones = correctMosaicIncosistencies(...
        conePositionsMicrons, coneSpacingsMicrons, ...
        thresholdFraction, 'cones');
    conePositionsMicrons = conePositionsMicrons(idxCones,:);
    coneSpacingsMicrons = coneSpacingsMicrons(idxCones);
    
    % Find RGCs within the roi
    idxRGC = positionsWithinROI(roi, RGCRFPositionsMicrons);
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(idxRGC,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(idxRGC);
    desiredConesToRGCratios = desiredConesToRGCratios(idxRGC);
    
    % Remove RGCs that may be too close to each other
    idxRGC = correctMosaicIncosistencies(RGCRFPositionsMicrons, ...
        RGCRFSpacingsMicrons, thresholdFraction, 'RGCs');
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(idxRGC,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(idxRGC);
    desiredConesToRGCratios = desiredConesToRGCratios(idxRGC);
    
end

function indices = positionsWithinROI(roi, positions)
    d = bsxfun(@minus,positions, roi.center);
    ecc = sqrt(sum(positions.^2,2));
    indices = find((abs(d(:,1)) <= 0.5*roi.size(1)) & (abs(d(:,2)) <= 0.5*roi.size(2)));
    
    % Re-order according to increasing eccentricity
    [~,sortedIdx] = sort(ecc(indices), 'ascend');
    indices = indices(sortedIdx);
end


function indices = correctMosaicIncosistencies(positions, spacings, thresholdFraction, unitType)
    micronsPerDegree = 300;
    ecc = sqrt(sum(positions.^2,2))/micronsPerDegree;
    eccRanges = [0 logspace(log10(0.1), log10(30), 30)];
    indices  = [];
    for iecc = 2:numel(eccRanges)
        idx = find((ecc >= eccRanges(iecc-1)) & (ecc < eccRanges(iecc)));
        meanSpacing = median(spacings(idx));
        posInEccBand = positions(idx,:);
        thresholdSeparation = thresholdFraction * meanSpacing;
        indicesToBeRemoved = indicesOfUnitsWithNearestSeparationLessThanThresholdSeparation(posInEccBand, thresholdSeparation);
        if (~isempty(indicesToBeRemoved))
            fprintf('Removing %d %s in ecc range %2.2f = %2.2f degs\n', numel(indicesToBeRemoved), unitType, eccRanges(iecc-1), eccRanges(iecc));
        end
        idx = setdiff(idx, idx(indicesToBeRemoved));
        indices = cat(1, indices , idx);
    end
end

function indices = indicesOfUnitsWithNearestSeparationLessThanThresholdSeparation(positions, thresholdSeparation)
    % Find the closest distance from each cone to all other cones
    allPairwiseDistances = pdist(positions);
    idx = find(allPairwiseDistances < thresholdSeparation);
    [unit1Indices,unit2Indices] = returnRowColFromLowerTriIndex(idx, size(positions,1));
    indices = unit1Indices;
    for k = 1:numel(indices)
        distEst1 = sqrt(sum((positions(unit1Indices(k),:)-positions(unit2Indices(k),:)).^2));
        distEst2 = allPairwiseDistances(idx(k));
        fprintf('\tElements %d and %d are too close. (distance = %2.2f, %2.2f). Removed %d.\n', ...
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
    p = pdist2(conePositions, conePositions, 'euclidean', 'Smallest', 5);
    p = p(2:end,:);
    coneSpacings = median(p,1);
end
