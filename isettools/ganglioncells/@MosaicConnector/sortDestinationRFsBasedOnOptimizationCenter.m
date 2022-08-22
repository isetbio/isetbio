function sortedIndices = sortDestinationRFsBasedOnOptimizationCenter(obj,unsortedIndices)

    centroidsOfDestinationRFsInThisGroup = ...
        obj.destinationRFcentroidsFromInputs(unsortedIndices,:);

    switch (obj.wiringParams.optimizationCenter)
        case 'origin'
            ecc = sum(centroidsOfDestinationRFsInThisGroup.^2,2);
        case 'latticeCenter'
            if (isempty(obj.sourceLatticeCenter))
                obj.sourceLatticeCenter = mean(obj.sourceLattice.RFpositionsMicrons,1);
            end

            diff = bsxfun(@minus, centroidsOfDestinationRFsInThisGroup, obj.sourceLatticeCenter);
            ecc = sum(diff.^2,2);
    end % switch

   % Compute sorted indices of destination RFs in increasing eccentricity
   [~, sortedIndices] = sort(ecc, 'ascend');

end
