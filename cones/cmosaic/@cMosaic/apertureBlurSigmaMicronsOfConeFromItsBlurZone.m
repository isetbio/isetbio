% Blur sigma (in microns) of a cone with index theConeIndex, from its blur zone
function apertureBlurSigmaMicrons = apertureBlurSigmaMicronsOfConeFromItsBlurZone(obj, theConeIndex)
    apertureBlurSigmaMicrons = [];
    if (isempty(obj.coneApertureModifiers))
        return;
    end
    foundZoneOfCone = false;
    for zoneIndex = 1:numel(obj.coneIndicesInZones)
        coneIndicesInThisZone = obj.coneIndicesInZones{zoneIndex};
        if (ismember(theConeIndex, coneIndicesInThisZone))
            apertureBlurSigmaMicrons = obj.blurApertureDiameterMicronsZones(zoneIndex)*obj.coneApertureModifiers.sigma;
            foundZoneOfCone = true;
        end  
        if (foundZoneOfCone)
            continue;
        end
    end
end
