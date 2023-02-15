function absorptionsRate = computeAbsorptionRate(obj, currentEMposMicrons, ...
    oiPositionsVectorsMicrons, absorptionsDensityImage,  ...
    oiResMicrons, lowOpticalImageResolutionWarning)
    
    if (~obj.employsImportedConeData)
        if (isempty(obj.oiResMicronsForZoning) || (oiResMicrons ~= obj.oiResMicronsForZoning))
            % Recompute cone apertures for current oiRes
            obj.oiResMicronsForZoning = oiResMicrons;
            obj.computeConeApertures(lowOpticalImageResolutionWarning);
        end
    end
    
    conesNum = numel(obj.coneApertureDiametersMicrons);
    coneTypesNum = size(absorptionsDensityImage,3);
    zonesNum = numel(obj.blurApertureDiameterMicronsZones);
    
    % Alocate memory
    absorptionsRateAllConeTypes = zeros(conesNum, size(absorptionsDensityImage,3));
    
    for zoneIndex = 1:zonesNum
        % Retrieve aperture size for this cone zone
        blurApertureDiameterMicrons = obj.blurApertureDiameterMicronsZones(zoneIndex);
        
        % Compute aperture kernel
        apertureKernel = obj.generateApertureKernel(blurApertureDiameterMicrons(1), oiResMicrons, lowOpticalImageResolutionWarning);

        % Determine which cones should receive this blur.
        coneIDsInZone = obj.coneIndicesInZones{zoneIndex};

        % Determine aperture areas in meters^2
        apertureAreasMetersSquared = obj.computeApertureAreasMetersSquared(coneIDsInZone);
        
        % Interpolate from oiPositions to conePositions at current emPos
        shiftedConePositions = bsxfun(@plus, obj.coneRFpositionsMicrons(coneIDsInZone,:),reshape(currentEMposMicrons, [1 2]));
        
        interpolationMethod = 'linear';
        extrapolationMethod = 'nearest';

        
        for coneTypeIndex = 1:coneTypesNum
            % Convolve with the cone aperture
            absorptionsDensityImageFiltered = conv2(squeeze(absorptionsDensityImage(:, :, coneTypeIndex)), apertureKernel, 'same');
           
            % Compute gridded interpolant for the original cone positions
            F = griddedInterpolant(oiPositionsVectorsMicrons, ...
                absorptionsDensityImageFiltered, interpolationMethod, extrapolationMethod);
            
            % Interpolate at current cone positions using gridded
            % interpolant, F, and multiply by aperture area
            absorptionsRateAllConeTypes(coneIDsInZone, coneTypeIndex) = ...
                F([shiftedConePositions(:,2) shiftedConePositions(:,1)]) .* apertureAreasMetersSquared;
        end
    end  % zoneIndex
    
    % Clear memory
    clear('absorptionsDensityImageFiltered')
    
    % Allocate new memory
    absorptionsRate = zeros(1,conesNum);
    
    for coneTypeIndex = 1:coneTypesNum
        % Get cone position indices for cones of this type
        switch (coneTypeIndex)
            case obj.LCONE_ID
                coneIndicesForSubmosaic = obj.lConeIndices;
            case obj.MCONE_ID
                coneIndicesForSubmosaic = obj.mConeIndices;
            case obj.SCONE_ID
                coneIndicesForSubmosaic = obj.sConeIndices;
            case obj.KCONE_ID
                coneIndicesForSubmosaic = obj.kConeIndices;
        end
        
        absorptionsRate(coneIndicesForSubmosaic) = absorptionsRateAllConeTypes(coneIndicesForSubmosaic, coneTypeIndex);
    end
    
    % Replace nan with zeros
    absorptionsRate(isnan(absorptionsRate)) = 0;
    
    % Account for decrease in outer segment length with eccentricity
    effectiveOSlengthAttenuationFactors = obj.computeEffectiveOSlengthAttenuationFactors(1:conesNum);
    if (~isempty(effectiveOSlengthAttenuationFactors))
        absorptionsRate = absorptionsRate .* effectiveOSlengthAttenuationFactors;
    end

end