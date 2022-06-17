function [cost, spatialVarianceCost, chromaticVarianceCost] = costToMaintainInputs(...
    obj, inputConeIndices, inputConeWeights, localRGCRFSpacingMicrons)

    if (isempty(localRGCRFSpacingMicrons))
        error('LocalRGCRFSpacingsMicrons has not been computed yet.')
    end

    if (isempty(inputConeIndices))
        fprintf(2,'\t\t* * * * - - - > Cost to maintain 0 cone indices: -99\n');
        cost = -99;
        spatialVarianceCost = -99;
        chromaticVarianceCost = -99;
        return;
    elseif (numel(inputConeIndices) == 1)
        % zero cost to maintain 1 cone
        cost = 0;
        spatialVarianceCost = 0;
        chromaticVarianceCost = 0;
        return;
    else
        inputConePositions = obj.inputConeMosaic.coneRFpositionsMicrons(inputConeIndices,:);
        inputConeTypes = obj.inputConeMosaic.coneTypes(inputConeIndices);
   
        switch (obj.wiringParams.spatialVarianceMetric)
            case 'maximal interinput distance'
                spatialVarianceCost = RGCconnector.maximalInterInputDistance(inputConePositions);
            case 'spatial variance'
                varianceXY = var(inputConePositions,inputConeWeights,1);
                spatialVarianceXY = varianceXY(:);
                spatialVarianceCost = sqrt(spatialVarianceXY(1)+spatialVarianceXY(2));
            otherwise
                error('Unknown spatialVarianceMetric: ''%s''.', obj.wiringParams.spatialVarianceMetric);
        end

        spatialVarianceCost = spatialVarianceCost / localRGCRFSpacingMicrons;

        lConeIndices = find(inputConeTypes == cMosaic.LCONE_ID);
        mConeIndices = find(inputConeTypes == cMosaic.MCONE_ID);
        lConeSignal = sum(inputConeWeights(lConeIndices));
        mConeSignal = sum(inputConeWeights(mConeIndices));
        
        totalLMConeSignal = lConeSignal + mConeSignal;
        if (lConeSignal <= mConeSignal)
            chromaticVarianceCost = lConeSignal/totalLMConeSignal;
        else
            chromaticVarianceCost = mConeSignal/totalLMConeSignal;
        end


        % Total cost
        w = obj.wiringParams.chromaticSpatialVarianceTradeoff;
        cost = w * spatialVarianceCost + (1-w)*chromaticVarianceCost;
    end
end