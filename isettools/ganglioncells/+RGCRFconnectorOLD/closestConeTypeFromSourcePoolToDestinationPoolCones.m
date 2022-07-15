function [sourceRGCconeIndex, idx] = closestConeTypeFromSourcePoolToDestinationPoolCones(...
        theSourceRGCinputConeIndices, theDestinationRGCinputConeIndices, allConeTypes, allConePositions)

    sourceRGCconeTypes = allConeTypes(theSourceRGCinputConeIndices);
    destinationRGCconeTypes = allConeTypes(theDestinationRGCinputConeIndices);

    lConeIndicesInDestinationRGC = find(destinationRGCconeTypes == cMosaic.LCONE_ID);
    mConeIndicesInDestinationRGC = find(destinationRGCconeTypes == cMosaic.MCONE_ID);
    lConesNumInDestinationRGC = numel(lConeIndicesInDestinationRGC);
    mConesNumInDestinationRGC = numel(mConeIndicesInDestinationRGC);

    if (lConesNumInDestinationRGC >= mConesNumInDestinationRGC)
        candidateSourceIndices = find(sourceRGCconeTypes == cMosaic.LCONE_ID);
        majorityDestinationIndices = lConeIndicesInDestinationRGC;
    else
        candidateSourceIndices = find(sourceRGCconeTypes == cMosaic.MCONE_ID);
        majorityDestinationIndices = mConeIndicesInDestinationRGC;
    end

    % If more than one of the source cones match the majority type in the
    % destination RGC, choose the one that is closest to those
    if (numel(candidateSourceIndices)>1)
        [~,idx] = RGCRFconnector.closestConePositionFromSourcePoolToDestinationPoolCones(...
                theSourceRGCinputConeIndices(candidateSourceIndices), ...
                theDestinationRGCinputConeIndices(majorityDestinationIndices), ...
                allConePositions);
        candidateSourceIndices = candidateSourceIndices(idx);
    elseif (numel(candidateSourceIndices)==0)
        [~,idx] = RGCRFconnector.closestConePositionFromSourcePoolToDestinationPoolCones(...
                theSourceRGCinputConeIndices, ...
                theDestinationRGCinputConeIndices(majorityDestinationIndices), ...
                allConePositions);
        candidateSourceIndices = idx;
    end

    sourceRGCconeIndex = theSourceRGCinputConeIndices(candidateSourceIndices);
    idx = candidateSourceIndices;
end


    