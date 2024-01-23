function apertureAreasMetersSquared = computeApertureAreasMetersSquared(obj, coneIndices)

    coneApertureDiametersMicrons = obj.coneApertureDiametersMicrons(coneIndices);
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
                characteristicRadiusMicrons =  obj.coneApertureToConeCharacteristicRadiusConversionFactor  * coneApertureDiametersMicrons;
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

end
