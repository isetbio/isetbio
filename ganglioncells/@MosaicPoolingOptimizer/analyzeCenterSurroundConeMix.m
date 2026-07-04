function [theCenterMajorityConeType, netCenterLconeWeight, netCenterMconeWeight, ...
         netSurroundLconeWeight, netSurroundMconeWeight, surroundConeMix, centerConeMix] = analyzeCenterSurroundConeMix(...
         theMRGCmosaic, theRGCindex, performSurroundAnalysisForConesExclusiveToTheSurround)

    [theCenterConeTypeNetWeights, ~, theCenterMajorityConeType, ...
        theCenterConeTypes, theCenterConeIndices] = theMRGCmosaic.centerConeTypeWeights(theRGCindex);
    netCenterLconeWeight = theCenterConeTypeNetWeights(find(theCenterConeTypes == cMosaic.LCONE_ID));
    netCenterMconeWeight = theCenterConeTypeNetWeights(find(theCenterConeTypes == cMosaic.MCONE_ID));

    if (performSurroundAnalysisForConesExclusiveToTheSurround)
        [~, theSurroundConeTypeNetWeights] = theMRGCmosaic.surroundConeTypeWeights(theRGCindex, theCenterConeIndices);
    else
        theSurroundConeTypeNetWeights = theMRGCmosaic.surroundConeTypeWeights(theRGCindex, []);
    end

    netSurroundLconeWeight = theSurroundConeTypeNetWeights(cMosaic.LCONE_ID);
    netSurroundMconeWeight = theSurroundConeTypeNetWeights(cMosaic.MCONE_ID);

    % Surround cone mix = net surround cone weights for the
    % non-dominant center cone / total surround L+M cone weight
    if (theCenterMajorityConeType == cMosaic.LCONE_ID)
        surroundConeMix = netSurroundMconeWeight / (netSurroundMconeWeight+netSurroundLconeWeight);
        centerConeMix = netCenterLconeWeight/(netCenterLconeWeight+netCenterMconeWeight);
    else
        surroundConeMix = netSurroundLconeWeight / (netSurroundMconeWeight+netSurroundLconeWeight);
        centerConeMix = netCenterMconeWeight/(netCenterLconeWeight+netCenterMconeWeight);
    end

end