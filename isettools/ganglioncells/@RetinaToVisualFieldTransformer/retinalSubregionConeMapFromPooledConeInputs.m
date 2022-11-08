function retinalSubregionConeMap = retinalSubregionConeMapFromPooledConeInputs(coneRc, coneDiameter, conePos, coneWeights, spatialSupport)
    
    [X,Y] = meshgrid(spatialSupport(:,1), spatialSupport(:,2));

    conesNumPooled = size(conePos,1);
    if (conesNumPooled == 0)
        retinalSubregionConeMap = X*0;
        return;
    end

    
    for iCone = 1:conesNumPooled
        d2 = (X-conePos(iCone,1)).^2 +  (Y-conePos(iCone,2)).^2;
        coneRc2 = (coneRc(iCone))^2;
        theConeApertureRF = exp(-d2/coneRc2);

        % Aperture RF must be zero outside the cone aperture
        %radiusSquared = (0.5*coneDiameter(iCone))^2;
        %idx = (d2 >= radiusSquared);
        %theConeApertureRF(idx) = 0;

        if (iCone == 1)
            retinalSubregionConeMap = coneWeights(iCone) * theConeApertureRF;
        else
            retinalSubregionConeMap = retinalSubregionConeMap + coneWeights(iCone) * theConeApertureRF;
        end
    end

end