function cost = chromaticVarianceCost(inputWeights, inputConeTypes)

    lConeIndices = find(inputConeTypes == cMosaic.LCONE_ID);
    mConeIndices = find(inputConeTypes == cMosaic.MCONE_ID);
    lConeSignal = sum(inputWeights(lConeIndices));
    mConeSignal = sum(inputWeights(mConeIndices));
    
    totalLMConeSignal = lConeSignal + mConeSignal;
    if (lConeSignal <= mConeSignal)
        cost = lConeSignal/totalLMConeSignal;
    else
        cost = mConeSignal/totalLMConeSignal;
    end

end
