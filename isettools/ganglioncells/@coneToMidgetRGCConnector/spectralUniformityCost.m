function theConeInputUniformityCost = spectralUniformityCost(inputConeTypes, inputConeWeights)
    

    % Compute the chromatic purity cost based on the relative L:Mcone inputs
    lConeIndices = find(inputConeTypes == cMosaic.LCONE_ID);
    mConeIndices = find(inputConeTypes == cMosaic.MCONE_ID);

    if (isempty(inputConeWeights))
        lConeSignal = numel(lConeIndices);
        mConeSignal = numel(mConeIndices);
    else
        lConeSignal = sum(inputConeWeights(lConeIndices));
        mConeSignal = sum(inputConeWeights(mConeIndices));
    end

    totalLMConeSignal = lConeSignal + mConeSignal;
    if (lConeSignal <= mConeSignal)
        % M-cone dominated, so theConePurityCost is the relative L-cone strength.
        % The higher the relative L-cone signal, the higher the conePurityCost,
        % Max value = 0.5 (equal L and Mcone inputs)
        theConeInputUniformityCost = lConeSignal/totalLMConeSignal;
    else
        % L-cone dominated, so theConePurityCost is the relative M-cone strength.
        % The higher the relative M-cone signal, the higher the conePurityCost,
        % Max value = 0.5 (equal L and Mcone inputs)
        theConeInputUniformityCost = mConeSignal/totalLMConeSignal;
    end

    % Put it in [0..1] range. 
    % So when theConeInputUniformityCost is 0, we have full cone input uniformity (single cone type)
    % When theConeInputUniformityCost ia 1, we have zero uniformity (equal L and M cones)
    theConeInputUniformityCost = 2*theConeInputUniformityCost;
end