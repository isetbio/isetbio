function cropMosaicToIncluceConesWithIndices(obj,keptConeIndices)
    
    % Update state
    obj.updateStateGivenKeptConeIndices(keptConeIndices);

    % Update cone-specific cone indices
    obj.lConeIndices = [];
    obj.mConeIndices = [];
    obj.sConeIndices = [];
    obj.kConeIndices = [];

    for iCone = 1:numel(keptConeIndices)
        theConeIndex = keptConeIndices(iCone);
        switch (obj.coneTypes(theConeIndex))
            case obj.LCONE_ID
                obj.lConeIndices(numel(obj.lConeIndices)+1) = iCone;
            case obj.MCONE_ID
                obj.mConeIndices(numel(obj.mConeIndices)+1) = iCone;
            case obj.SCONE_ID
                obj.sConeIndices(numel(obj.sConeIndices)+1) = iCone;
            case obj.KCONE_ID
                obj.kConeIndices(numel(obj.kConeIndices)+1) = iCone;

        end
    end

    % Compute achieved cone densities
    conesNum = obj.conesNum;

    % Assign cone types 
    obj.coneTypes = zeros(conesNum,1);
    obj.coneTypes(obj.lConeIndices) = obj.LCONE_ID;
    obj.coneTypes(obj.mConeIndices) = obj.MCONE_ID;
    obj.coneTypes(obj.sConeIndices) = obj.SCONE_ID;
    obj.coneTypes(obj.kConeIndices) = obj.KCONE_ID;
    
    % Reshape indices
    obj.lConeIndices = reshape(obj.lConeIndices, [numel(obj.lConeIndices) 1]);
    obj.mConeIndices = reshape(obj.mConeIndices, [numel(obj.mConeIndices) 1]);
    obj.sConeIndices = reshape(obj.sConeIndices, [numel(obj.sConeIndices) 1]);
    obj.kConeIndices = reshape(obj.kConeIndices, [numel(obj.kConeIndices) 1]);

    achievedConeDensities = [...
        numel(obj.lConeIndices)/conesNum ...
        numel(obj.mConeIndices)/conesNum ...
        numel(obj.sConeIndices)/conesNum ...
        numel(obj.kConeIndices)/conesNum];
    
    % Update coneDensities
    obj.achievedConeDensities = achievedConeDensities;

    % Update min and max cone position
    obj.minRFpositionMicrons = squeeze(min(obj.coneRFpositionsMicrons,[],1));
    obj.maxRFpositionMicrons = squeeze(max(obj.coneRFpositionsMicrons,[],1));
            
    % Update photon absorption attenuation factors to account for
    % the decrease in outer segment legth with ecc.
    obj.outerSegmentLengthEccVariationAttenuationFactors = obj.outerSegmentLengthEccVariationAttenuationFactors(keptConeIndices);

end
