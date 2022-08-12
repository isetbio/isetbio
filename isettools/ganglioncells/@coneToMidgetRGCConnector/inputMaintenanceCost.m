function theCostComponents = inputMaintenanceCost(obj, ...
    inputIndices, inputWeights, destinationRFspacing)
    
    if (isempty(inputIndices))
        fprintf(2,'\t\t* * * * - - - > Cost to maintain 0 inputs: -99\n');
        theCostComponents = [-99 -99 -99];
        return;
    end

    if (numel(inputIndices) == 1)
        % zero cost to maintain 1 cone
        theCostComponents = [0 0 0];
        return;
    else
        % Retrieve the input positions
        inputPositions = obj.sourceLattice.RFpositionsMicrons(inputIndices,:);

        % Retrieve the input cone types
        if (isfield(obj.sourceLattice, 'metaData')) && ...
           (isfield(obj.sourceLattice.metaData, 'coneTypes')) && ...
           (isfield(obj.sourceLattice.metaData, 'coneTypeIDs'))
                inputConeTypes = obj.sourceLattice.metaData.coneTypes(inputIndices);
        else
                %fprintf(2,'No metaData were passed for the cone mosaic. Assuming all L-cones.\n')
                inputConeTypes = cMosaic.LCONE_ID;
        end

        switch (obj.wiringParams.spatialVarianceMetric)
            case 'maximal interinput distance'
                spatialVarianceCost = MosaicConnector.maximalInterInputDistance(inputPositions);
            case 'spatial variance'
                varianceXY = var(inputPositions,inputWeights,1);
                spatialVarianceXY = varianceXY(:);
                spatialVarianceCost = sqrt(spatialVarianceXY(1)+spatialVarianceXY(2));
            otherwise
                error('Unknown spatialVarianceMetric: ''%s''.', obj.wiringParams.spatialVarianceMetric);
        end

        spatialVarianceCost = spatialVarianceCost / destinationRFspacing;

        lConeIndices = find(inputConeTypes == cMosaic.LCONE_ID);
        mConeIndices = find(inputConeTypes == cMosaic.MCONE_ID);
        lConeSignal = sum(inputWeights(lConeIndices));
        mConeSignal = sum(inputWeights(mConeIndices));
        
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

    % Form vector with cost components
    theCostComponents = [cost spatialVarianceCost chromaticVarianceCost];
end
