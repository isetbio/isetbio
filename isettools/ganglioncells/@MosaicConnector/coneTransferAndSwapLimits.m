% Method to determine the max number of cones that can be
% transferred or swapped between nearby RF centers, depending on
% the mosaic eccentricity

function [maxConeInputsPerRGCToConsiderTransferToNearbyRGCs, maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs] = ...
    coneTransferAndSwapLimits(mosaicEccDegs)

    xEccDegs = abs(round(mosaicEccDegs(1)));
    if (xEccDegs <= 2)
        maxConeInputsPerRGCToConsiderTransferToNearbyRGCs = 4; 
        maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs = 4;
    elseif (xEccDegs <= 4)
        maxConeInputsPerRGCToConsiderTransferToNearbyRGCs = 10;
        maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs = 10;
    elseif (xEccDegs <= 10)
        maxConeInputsPerRGCToConsiderTransferToNearbyRGCs = 10;
        maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs = 10;
    elseif (xEccDegs <= 14)
        maxConeInputsPerRGCToConsiderTransferToNearbyRGCs = 15;
        maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs = 15;
    else
        maxConeInputsPerRGCToConsiderTransferToNearbyRGCs = 30;
        maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs = 30;
    end

end
