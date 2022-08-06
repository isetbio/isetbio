function extraMicronsForSurroundCones = estimateMaxSurroundCharacteristicRadiusMicrons(mosaicEccMicrons, mosaicSizeMicrons)

    mosaicEccDegs = WatsonRGCModel.rhoMMsToDegs(mosaicEccMicrons/1000);
    mosaicEccDegs = sqrt(sum(mosaicEccDegs.^2,2));
    surroundCharacteristicRadiusDegs(1) = max([0.1 0.203 * mosaicEccDegs.^0.472]);
    
    minMosaicEccDegs = WatsonRGCModel.rhoMMsToDegs((mosaicEccMicrons-0.5*max(mosaicSizeMicrons))/1000);
    minMosaicEccDegs = sqrt(sum(minMosaicEccDegs.^2,2));
    surroundCharacteristicRadiusDegs(2) = max([0.1 0.203 * minMosaicEccDegs .^0.472]);
    
    maxMosaicEccDegs = WatsonRGCModel.rhoMMsToDegs((mosaicEccMicrons+0.5*max(mosaicSizeMicrons))/1000);
    maxMosaicEccDegs = sqrt(sum(maxMosaicEccDegs.^2,2));
    surroundCharacteristicRadiusDegs(3) = max([0.1 0.203 * maxMosaicEccDegs .^0.472]);
    
    extraMicronsForSurroundCones = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(max(surroundCharacteristicRadiusDegs), mosaicEccDegs);
end

