function [sourceRGCconeIndex, idx] = closestConePositionFromSourcePoolToDestinationPoolCones(...
        theSourceRGCinputConeIndices, theDestinationRGCinputConeIndices, allConePositions)

    sourceRGCconePositions = allConePositions(theSourceRGCinputConeIndices,:);
    destinationRGCconePositions = allConePositions(theDestinationRGCinputConeIndices,:);

    for iSourceConeIndex = 1:numel(theSourceRGCinputConeIndices)
        d = pdist2(destinationRGCconePositions, sourceRGCconePositions(iSourceConeIndex,:));
        [minDistance(iSourceConeIndex), bestTargetConeIndex] = min(d);
    end

   [~,idx] = min(minDistance);
   sourceRGCconeIndex = theSourceRGCinputConeIndices(idx);
end


    