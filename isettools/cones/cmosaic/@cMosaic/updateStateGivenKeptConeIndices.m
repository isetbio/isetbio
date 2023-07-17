function updateStateGivenKeptConeIndices(obj, keptConeIndices)
    
    obj.coneRFpositionsDegs = obj.coneRFpositionsDegs(keptConeIndices,:);
    obj.coneRFpositionsMicrons = obj.coneRFpositionsMicrons(keptConeIndices,:);
    
    obj.coneRFspacingsDegs = obj.coneRFspacingsDegs(keptConeIndices);
    obj.coneRFspacingsMicrons = obj.coneRFspacingsMicrons(keptConeIndices);
    
    obj.coneApertureDiametersMicrons = obj.coneApertureDiametersMicrons(keptConeIndices);
    obj.coneApertureDiametersDegs = obj.coneApertureDiametersDegs(keptConeIndices);
    
    if (~isempty(obj.coneApertureRodIntrusionInducedShrinkageFactors))
        obj.coneApertureRodIntrusionInducedShrinkageFactors = obj.coneApertureRodIntrusionInducedShrinkageFactors(keptConeIndices);
    end

    minRFpositionDegs = squeeze(min(obj.coneRFpositionsDegs,[],1));
    maxRFpositionDegs = squeeze(max(obj.coneRFpositionsDegs,[],1));
    minRFpositionMicrons = squeeze(min(obj.coneRFpositionsMicrons,[],1));
    maxRFpositionMicrons = squeeze(max(obj.coneRFpositionsMicrons,[],1));

    obj.eccentricityDegs = 0.5*(maxRFpositionDegs+minRFpositionDegs);
    obj.eccentricityMicrons = 0.5*(maxRFpositionMicrons+minRFpositionMicrons);

    obj.sizeDegs = maxRFpositionDegs-minRFpositionDegs;

    keptConeIndices = keptConeIndices(:);
    theConeIndices = obj.coneIndicesInZones;
    parfor zoneBandIndex = 1:numel(theConeIndices)
        coneIndices = theConeIndices{zoneBandIndex};
        newConeIndices = [];
        for k = 1:numel(coneIndices)
            newConeIndices = cat(1, newConeIndices, find(keptConeIndices == coneIndices(k)));
        end
        theConeIndices{zoneBandIndex} = newConeIndices;
    end
    obj.coneIndicesInZones = theConeIndices;
end