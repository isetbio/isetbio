function theCostComponents = inputMaintenanceCost(obj, ...
    inputIndices, inputWeights, destinationRFspacing)
    
    if (isempty(inputIndices))
        %fprintf(2,'\t\t* * * * - - - > Cost to maintain 0 inputs: -99\n');
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
        if (obj.coneTypeInfoIsAvailable)
           inputConeTypes = obj.sourceLattice.metaData.coneTypes(inputIndices);
        else
            %fprintf(2,'No metaData were passed for the cone mosaic. Assuming all L-cones.\n')
            inputConeTypes = cMosaic.LCONE_ID;
        end

        % Spatial variance cost
        spatialVarianceCost = coneToMidgetRGCConnector.spatialVarianceCost(...
            obj.wiringParams.spatialVarianceMetric, inputWeights, inputPositions, destinationRFspacing);

        % Chromatic variance cost
        chromaticVarianceCost = coneToMidgetRGCConnector.chromaticVarianceCost(inputWeights, inputConeTypes);

        % Total cost
        w = obj.wiringParams.chromaticSpatialVarianceTradeoff;
        cost = w * spatialVarianceCost + (1-w)*chromaticVarianceCost;
    end

    % Form vector with cost components
    theCostComponents = [cost spatialVarianceCost chromaticVarianceCost];
end
