function absorptionsRate = computeAbsorptionRate(obj, currentEMposMicrons, oiPositionsMicrons, absorptionsDensityImage,  ...
    oiResMicrons, coneApertureDiametersMicrons, ...
    coneIndicesInZones)

    conesNum = numel(coneApertureDiametersMicrons);
    coneTypesNum = size(absorptionsDensityImage,3);
    zonesNum = numel(obj.blurApertureDiameterMicronsZones);
    
    % Alocate memory
    absorptionsDensityImageFiltered = absorptionsDensityImage * 0;
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
        

        % Do the convolution
        for coneType = 1:coneTypesNum
            absorptionsDensityImageFiltered(:, :, coneType) = conv2(squeeze(absorptionsDensityImage(:, :, coneType)), apertureKernel, 'same');
        end
            
        % Determine which cones should receive this blur.
        coneIDsInZone = coneIndicesInZones{zoneIndex};
        
        % Determine aperture area for these cones
        apertureAreasMetersSquared = ((pi * (0.5*coneApertureDiametersMicrons(coneIDsInZone)*1e-6).^2))';

        % Interpolate from oiPositions to conePositions at current emPos
        shiftedConePositions = bsxfun(@plus, obj.coneRFpositionsMicrons(coneIDsInZone,:),reshape(currentEMposMicrons, [1 2]));
        
        for coneTypeIndex = 1:coneTypesNum
            fullImage = squeeze(absorptionsDensityImageFiltered(:,:,coneTypeIndex));
            F = griddedInterpolant({oiPositionsMicrons(:,1), oiPositionsMicrons(:,2)}, fullImage', 'linear');
            vq = F([shiftedConePositions(:,1) shiftedConePositions(:,2)]);
            absorptionsRateAllConeTypes(coneIDsInZone, coneTypeIndex) = vq .* apertureAreasMetersSquared;
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