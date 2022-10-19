function retinalSubregionConeMap = retinalSubregionConeMapFromPooledConeInputs(coneRc, conePos, coneWeights, spatialSupport)
    
    [X,Y] = meshgrid(spatialSupport(:,1), spatialSupport(:,2));

    conesNumPooled = size(conePos,1);
    if (conesNumPooled == 0)
        retinalSubregionConeMap = X*0;
        return;
    end

    
    for iCone = 1:conesNumPooled
        theConeApertureRF = exp(-((X-conePos(iCone,1))/coneRc(iCone)).^2) .* ...
                            exp(-((Y-conePos(iCone,2))/coneRc(iCone)).^2);
        if (iCone == 1)
            retinalSubregionConeMap = coneWeights(iCone) * theConeApertureRF;
        else
            retinalSubregionConeMap = retinalSubregionConeMap + coneWeights(iCone) * theConeApertureRF;
        end
    end

end