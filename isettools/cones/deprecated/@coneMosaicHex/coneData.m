function mosaicMetaData = coneData(obj)
        
    mosaicMetaData.positionUnits = 'microns';
    mosaicMetaData.blurApertureDiameterMicrons = 2*sqrt(obj.pigment.pdArea / pi)*1e6;

    sampledHexMosaicXaxis = obj.patternSupport(1, :, 1);
    sampledHexMosaicYaxis = obj.patternSupport(:, 1, 2);

    sampledHexMosaicXaxis = sampledHexMosaicXaxis - mean(sampledHexMosaicXaxis);
    sampledHexMosaicYaxis = sampledHexMosaicYaxis - mean(sampledHexMosaicYaxis);

    % Compute correction factors
    [correctionFactors, outerSegmentLengthAttenationFactors, innerSegmentDiameterBoostFactors] = ...
        coneMosaicHex.computeConeEfficiencyCorrectionFactors(obj, 'coneData');
    
    % L-cones
    idx = find(obj.pattern == 2);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(end-iRows+1);
    lConeCoords = [coneXcoords(:) coneYcoords(:)];
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        lConeApertures = 2*sqrt(obj.pigment.pdArea / pi * innerSegmentDiameterBoostFactors(idx));
    else
        lConeApertures = 2*sqrt(obj.pigment.pdArea / pi) * ones(1, size(lConeCoords,1)); 
    end

    if (obj.eccBasedConeQuantalEfficiency)
        lConeOuterSegmentLengthAttenationFactors = outerSegmentLengthAttenationFactors(idx);
        lConeTotalCorrectionFactors = correctionFactors(idx);
    else
        lConeOuterSegmentLengthAttenationFactors = ones(1, size(lConeCoords,1));
        lConeTotalCorrectionFactors = ones(1, size(lConeCoords,1));
    end
    
    % M-cones
    idx = find(obj.pattern == 3);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(end-iRows+1);
    mConeCoords = [coneXcoords(:) coneYcoords(:)];
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        mConeApertures = 2*sqrt(obj.pigment.pdArea / pi * innerSegmentDiameterBoostFactors(idx));
    else
        mConeApertures = 2*sqrt(obj.pigment.pdArea / pi) * ones(1, size(mConeCoords,1)); 
    end
    if (obj.eccBasedConeQuantalEfficiency)
        mConeOuterSegmentLengthAttenationFactors = outerSegmentLengthAttenationFactors(idx);
        mConeTotalCorrectionFactors = correctionFactors(idx);
    else
        mConeOuterSegmentLengthAttenationFactors = ones(1, size(mConeCoords,1));
        mConeTotalCorrectionFactors = ones(1, size(mConeCoords,1));
    end
    
    % S-cones
    idx = find(obj.pattern == 4);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(end-iRows+1);
    sConeCoords = [coneXcoords(:) coneYcoords(:)];
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        sConeApertures = 2*sqrt(obj.pigment.pdArea / pi * innerSegmentDiameterBoostFactors(idx));
    else
        sConeApertures = 2*sqrt(obj.pigment.pdArea / pi) * ones(1, size(sConeCoords,1)); 
    end
    if (obj.eccBasedConeQuantalEfficiency)
        sConeOuterSegmentLengthAttenationFactors = outerSegmentLengthAttenationFactors(idx);
        sConeTotalCorrectionFactors = correctionFactors(idx);
    else
        sConeOuterSegmentLengthAttenationFactors = ones(1, size(sConeCoords,1));
        sConeTotalCorrectionFactors = ones(1, size(sConeCoords,1));
    end
    
    
    % Assemble data from all cones
    % L-cones
    mosaicMetaData.positions = lConeCoords * 1e6;
    mosaicMetaData.types = ones(size(lConeCoords,1),1);
    mosaicMetaData.lightGatheringApertureDiameters = lConeApertures(:)*1e6;
    mosaicMetaData.outerSegmentLengthAttenationFactors = lConeOuterSegmentLengthAttenationFactors(:);
    mosaicMetaData.totalAttenationFactors = lConeTotalCorrectionFactors(:);
 
    % M-cones
    mosaicMetaData.positions = cat(1, mosaicMetaData.positions, mConeCoords * 1e6);
    mosaicMetaData.types = cat(1, mosaicMetaData.types, 2*ones(size(mConeCoords,1),1));
    mosaicMetaData.lightGatheringApertureDiameters = cat(1, mosaicMetaData.lightGatheringApertureDiameters, mConeApertures(:)*1e6);
    mosaicMetaData.outerSegmentLengthAttenationFactors = cat(1, mosaicMetaData.outerSegmentLengthAttenationFactors, mConeOuterSegmentLengthAttenationFactors(:));
    mosaicMetaData.totalAttenationFactors = cat(1, mosaicMetaData.totalAttenationFactors, mConeTotalCorrectionFactors(:));
    
    % S-cones
    mosaicMetaData.positions = cat(1, mosaicMetaData.positions, sConeCoords * 1e6);
    mosaicMetaData.types = cat(1, mosaicMetaData.types, 3*ones(size(sConeCoords,1),1));
    mosaicMetaData.lightGatheringApertureDiameters = cat(1, mosaicMetaData.lightGatheringApertureDiameters, sConeApertures(:)*1e6);
    mosaicMetaData.outerSegmentLengthAttenationFactors = cat(1, mosaicMetaData.outerSegmentLengthAttenationFactors, sConeOuterSegmentLengthAttenationFactors(:));
    mosaicMetaData.totalAttenationFactors = cat(1, mosaicMetaData.totalAttenationFactors, sConeTotalCorrectionFactors(:));
end

