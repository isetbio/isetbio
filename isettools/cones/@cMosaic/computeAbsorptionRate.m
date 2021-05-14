function absorptionsRate = computeAbsorptionRate(obj, currentEMposMicrons, oiPositionsVectorsMicrons, absorptionsDensityImage,  ...
    oiResMicrons, coneApertureDiametersMicrons, ...
    coneIndicesInZones)

    conesNum = numel(coneApertureDiametersMicrons);
    coneTypesNum = size(absorptionsDensityImage,3);
    zonesNum = numel(obj.blurApertureDiameterMicronsZones);
    
    % Alocate memory
    absorptionsRateAllConeTypes = zeros(conesNum, size(absorptionsDensityImage,3));
    
    for zoneIndex = 1:zonesNum
        % Retrieve aperture size for this cone zone
        blurApertureDiameterMicrons = obj.blurApertureDiameterMicronsZones(zoneIndex);
        
        % Compute aperture kernel
        if (numel(obj.blurApertureDiameterMicronsZones) == 1)
            % Do not include the effect of spatial filtering via
            % ecc-dependent cone aperture. Instead do a single concolution
            % with the median cone aperture.
            apertureKernel = cMosaic.generateApertureKernel(blurApertureDiameterMicrons(1), oiResMicrons);
            
        else
            % Include the effect of spatial filtering via
            % ecc-dependent cone aperture, different aperture for each zone
            apertureKernel = cMosaic.generateApertureKernel(blurApertureDiameterMicrons, oiResMicrons);
        end
        
        % Determine which cones should receive this blur.
        coneIDsInZone = coneIndicesInZones{zoneIndex};
        
        % Determine aperture area for these cones
        apertureAreasMetersSquared = ((pi * (0.5*coneApertureDiametersMicrons(coneIDsInZone)*1e-6).^2))';

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
    if (isempty(obj.importedOSLengthAttenuationFactors))
        if (obj.eccVaryingOuterSegmentLength)
            % Do it for each cone separately based  on its eccentricity
            absorptionsRate = absorptionsRate .* obj.outerSegmentLengthEccVariationAttenuationFactors;
        else
            % Do it for alls cones based on the median attenuation factor 
            absorptionsRate = absorptionsRate * median(obj.outerSegmentLengthEccVariationAttenuationFactors);
        end
    else
        if (obj.eccVaryingOuterSegmentLength)
            absorptionsRate = absorptionsRate .* obj.importedOSLengthAttenuationFactors;
        end
    end
    
end