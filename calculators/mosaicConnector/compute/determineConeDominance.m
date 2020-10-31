function LMconeBalance = determineConeDominance(RGCconeInputInfo)
    global LCONE_ID
    global MCONE_ID
    global SCONE_ID
    
    if (isempty(RGCconeInputInfo))
        LMconeBalance = [];
        return;
    end
    
    rgcsNum = numel(RGCconeInputInfo);

    LMconeBalance.center = zeros(1, rgcsNum);
    LMconeBalance.surround = zeros(1, rgcsNum);
    
    for mRGCindex = 1:rgcsNum
        coneTypes = RGCconeInputInfo{mRGCindex}.center.types;
        coneWeights = RGCconeInputInfo{mRGCindex}.center.weights;
        
        lConeIndices = find(coneTypes == LCONE_ID);
        mConeIndices = find(coneTypes == MCONE_ID);
        
        lConeSignal = sum(coneWeights(lConeIndices));
        mConeSignal = sum(coneWeights(mConeIndices));
        LMconeBalance.center(mRGCindex) = lConeSignal / (lConeSignal+mConeSignal);
        
        coneTypes = RGCconeInputInfo{mRGCindex}.surround.types;
        coneWeights = RGCconeInputInfo{mRGCindex}.surround.weights;
        
        lConeIndices = find(coneTypes == LCONE_ID);
        lConeSignal = sum(coneWeights(lConeIndices));
        mConeIndices = find(coneTypes == MCONE_ID);
        mConeSignal = sum(coneWeights(mConeIndices));
        LMconeBalance.surround(mRGCindex) = lConeSignal / (lConeSignal+mConeSignal);
        
    end
    
end

