function absorptionsRate = computeAbsorptionRate(obj, currentEMposMicrons, oiPositionsVectorsMicrons, absorptionsDensityImage,  ...
    oiResMicrons)
    
    if (isempty(obj.oiResMicronsForZoning) || (oiResMicrons ~= obj.oiResMicronsForZoning))
        % Recompute cone apertures for current oiRes
        obj.oiResMicronsForZoning = oiResMicrons;
        obj.computeConeApertures();
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
        apertureKernel = obj.generateApertureKernel(blurApertureDiameterMicrons(1), oiResMicrons);

        % Determine which cones should receive this blur.
        coneIDsInZone = obj.coneIndicesInZones{zoneIndex};

        
        % Determine aperture area
        coneApertureDiametersMicrons = obj.coneApertureDiametersMicrons(coneIDsInZone);
        if (~obj.eccVaryingConeAperture)
            coneApertureDiametersMicrons = coneApertureDiametersMicrons*0 + median(coneApertureDiametersMicrons);
        end

        if (obj.anchorAllEccVaryingParamsToTheirFovealValues)
            % Foveal value of cone aperture
            coneApertureDiametersMicrons = coneApertureDiametersMicrons*0 + min(obj.coneApertureDiametersMicrons(:));
        end

        if (isfield(obj.coneApertureModifiers, 'shape'))
            switch (obj.coneApertureModifiers.shape)
                case 'Gaussian'
                    gaussianSigmaMicrons =  obj.coneApertureModifiers.sigma * coneApertureDiametersMicrons;
                    characteristicRadiusMicrons = gaussianSigmaMicrons  * sqrt(2.0);
                    apertureAreasMetersSquared = ((pi * (characteristicRadiusMicrons*1e-6).^2))';
                case 'Pillbox'
                    apertureAreasMetersSquared = ((pi * (0.5*coneApertureDiametersMicrons*1e-6).^2))';
                otherwise
                    error('Do not know how to generate a ''%s'' aperture.', obj.coneApertureModifiers.apertureShape)
            end
        else
            % By default we are using the pillbox aperture area
            apertureAreasMetersSquared = ((pi * (0.5*coneApertureDiametersMicrons*1e-6).^2))';
        end

        
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
    if (~obj.anchorAllEccVaryingParamsToTheirFovealValues)
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

    
end