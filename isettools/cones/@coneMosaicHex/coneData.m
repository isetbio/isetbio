function mosaicMetaData = coneData(obj)

     % Account for the fact mismatch between obj.pigment.pdArea and the
     % actual area at 0,0 degs
    [~,apertureMetersAtZeroEcc, ~] = coneSizeReadData(...
            'eccentricity', 0.0,...
            'angle', 0.0);
        
    mosaicMetaData.positionUnits = 'microns';
    mosaicMetaData.boostDueToMismatchBetweenConeAreaAtZeroEccAndPigmentPDarea = obj.pigment.pdArea / (pi*(apertureMetersAtZeroEcc/2).^2);
       
    sampledHexMosaicXaxis = obj.patternSupport(1, :, 1);
    sampledHexMosaicYaxis = obj.patternSupport(:, 1, 2);


    [correctionFactors, outerSegmentLengthAttenationFactors] = coneMosaicHex.computeConeEfficiencyCorrectionFactors(obj, 'coneData');
    
    % L-cones
    idx = find(obj.pattern == 2);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(iRows);
    lConeCoords = [coneXcoords(:) coneYcoords(:)];
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [~, ~, ~, lConeApertures] = coneMosaicHex.computeApertureSizes(...
            [],[], ...
            [],[], ...
            coneXcoords, coneYcoords ...
        );
    else
        [~, apertureMeters, ~] = coneSizeReadData(...
        'eccentricity',0, 'angle',0);
        lConeApertures = apertureMeters * ones(1, size(lConeCoords,1)); 
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
    coneYcoords = sampledHexMosaicYaxis(iRows);
    mConeCoords = [coneXcoords(:) coneYcoords(:)];
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [~, ~, ~, mConeApertures] = coneMosaicHex.computeApertureSizes(...
            [],[], ...
            [],[], ...
            coneXcoords, coneYcoords ...
        );
    else
        [~, apertureMeters, ~] = coneSizeReadData(...
        'eccentricity',0, 'angle',0);
        mConeApertures = apertureMeters * ones(1, size(mConeCoords,1)); 
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
    coneYcoords = sampledHexMosaicYaxis(iRows);
    sConeCoords = [coneXcoords(:) coneYcoords(:)];
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [~, ~, ~, sConeApertures] = coneMosaicHex.computeApertureSizes(...
            [],[], ...
            [],[], ...
            coneXcoords, coneYcoords ...
        );
    else
        [~, apertureMeters, ~] = coneSizeReadData(...
        'eccentricity',0, 'angle',0);
        sConeApertures = apertureMeters * ones(1, size(sConeCoords,1)); 
    end
    if (obj.eccBasedConeQuantalEfficiency)
        sConeOuterSegmentLengthAttenationFactors = outerSegmentLengthAttenationFactors(idx);
        sConeTotalCorrectionFactors = correctionFactors(idx);
    else
        sConeOuterSegmentLengthAttenationFactors = ones(1, size(sConeCoords,1));
        sConeTotalCorrectionFactors = ones(1, size(sConeCoords,1));
    end
    
    
    % Assemble all cones
    % add the L-cones
    mosaicMetaData.positions = lConeCoords * 1e6;
    mosaicMetaData.types = ones(size(lConeCoords,1),1);
    mosaicMetaData.apertures = lConeApertures'*1e6;
    mosaicMetaData.outerSegmentLengthAttenationFactors = lConeOuterSegmentLengthAttenationFactors(:);
    mosaicMetaData.totalAttenationFactors = lConeTotalCorrectionFactors(:);
 
    % add the M-cones
    mosaicMetaData.positions = cat(1, mosaicMetaData.positions, mConeCoords * 1e6);
    mosaicMetaData.types = cat(1, mosaicMetaData.types, 2*ones(size(mConeCoords,1),1));
    mosaicMetaData.apertures = cat(1, mosaicMetaData.apertures, mConeApertures'*1e6);
    mosaicMetaData.outerSegmentLengthAttenationFactors = cat(1, mosaicMetaData.outerSegmentLengthAttenationFactors, mConeOuterSegmentLengthAttenationFactors(:));
    mosaicMetaData.totalAttenationFactors = cat(1, mosaicMetaData.totalAttenationFactors, mConeTotalCorrectionFactors(:));
    
    % add the S-cones
    mosaicMetaData.positions = cat(1, mosaicMetaData.positions, sConeCoords * 1e6);
    mosaicMetaData.types = cat(1, mosaicMetaData.types, 3*ones(size(sConeCoords,1),1));
    mosaicMetaData.apertures = cat(1, mosaicMetaData.apertures, sConeApertures'*1e6);
    mosaicMetaData.outerSegmentLengthAttenationFactors = cat(1, mosaicMetaData.outerSegmentLengthAttenationFactors, sConeOuterSegmentLengthAttenationFactors(:));
    mosaicMetaData.totalAttenationFactors = cat(1, mosaicMetaData.totalAttenationFactors, sConeTotalCorrectionFactors(:));
end

