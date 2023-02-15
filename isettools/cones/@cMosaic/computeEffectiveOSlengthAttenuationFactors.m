function osLengthAttenuationFactors = computeEffectiveOSlengthAttenuationFactors(obj, coneIndices)

    osLengthAttenuationFactors = [];

    if (~obj.anchorAllEccVaryingParamsToTheirFovealValues)
        if (isempty(obj.importedOSLengthAttenuationFactors))
            if (obj.eccVaryingOuterSegmentLength)
                % Do it for each cone separately based  on its eccentricity
                v = obj.outerSegmentLengthEccVariationAttenuationFactors;
                osLengthAttenuationFactors = v(coneIndices);
            else
                % Do it for alls cones based on the median attenuation factor 
                v = median(obj.outerSegmentLengthEccVariationAttenuationFactors);
                osLengthAttenuationFactors = v(coneIndices);
            end
        else
            if (obj.eccVaryingOuterSegmentLength)
                v = obj.importedOSLengthAttenuationFactors;
                osLengthAttenuationFactors = v(coneIndices);
            end
        end
    end

    
end
