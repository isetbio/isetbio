function mosaicMetaData = coneData(obj)

     % Account for the fact mismatch between obj.pigment.pdArea and the
     % actual area at 0,0 degs
    [~,apertureMetersAtZeroEcc, ~] = coneSizeReadData(...
            'eccentricity', 0.0,...
            'angle', 0.0);
    mosaicMetaData.absorptionCorrectionFactor = obj.pigment.pdArea / (pi*(apertureMetersAtZeroEcc/2).^2);
       
    sampledHexMosaicXaxis = obj.patternSupport(1, :, 1);
    sampledHexMosaicYaxis = obj.patternSupport(:, 1, 2);

    
    % L-cones
    idx = find(obj.pattern == 2);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(iRows);
    mosaicMetaData.lConeCoords = [coneXcoords(:) coneYcoords(:)];
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [~, ~, ~, mosaicMetaData.lConeApertures] = coneMosaicHex.computeApertureSizes(...
            [],[], ...
            [],[], ...
            coneXcoords, coneYcoords ...
        );
    end
    
    
    % M-cones
    idx = find(obj.pattern == 3);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(iRows);
    mosaicMetaData.mConeCoords = [coneXcoords(:) coneYcoords(:)];
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [~, ~, ~, mosaicMetaData.mConeApertures] = coneMosaicHex.computeApertureSizes(...
           [],[], ...
            [],[], ...
            coneXcoords, coneYcoords ...
        );
    end
    
    
    % S-cones
    idx = find(obj.pattern == 4);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx);
    coneXcoords = sampledHexMosaicXaxis(iCols);
    coneYcoords = sampledHexMosaicYaxis(iRows);
    mosaicMetaData.sConeCoords = [coneXcoords(:) coneYcoords(:)];
    
    if (obj.shouldCorrectAbsorptionsWithEccentricity())
        [~, ~, ~, mosaicMetaData.sConeApertures] = coneMosaicHex.computeApertureSizes(...
           [],[], ...
            [],[], ...
            coneXcoords, coneYcoords ...
        );
    end
    
    
    
end

