function [RGCRFinputs, RGCRFweights, RGCRFcentroids, RGCRFspacings] = ...
            transferSpecificConeInputFromSourceRGCToNearbyRGC(...
                RGCRFinputs, RGCRFweights, RGCRFcentroids, ...
                theSourceRGCindex, coneIndexToTransfer, ...
                localConeToRGCDensityRatio, ...
                wiringParams, ...
                allConePositions, allConeSpacings, allConeTypes)
    
    % Nearby RGC indices
    [~, nearbyRGCindices] = pdist2(RGCRFcentroids, allConePositions(coneIndexToTransfer,:), '', ...
           'smallest', wiringParams.maxNearbyRGCsNum);

    % Exclude ourself
    nearbyRGCindices = setdiff(nearbyRGCindices, theSourceRGCindex);

    % We want to transfer coneIndexToTransfer to one of the
    % nearbyRGCindices

    for kk = 1:numel(nearbyRGCindices)
        % The current inputs and weights
        theExistingConeInputs = RGCRFinputs{nearbyRGCindices(kk)};
        theExistingConeWeights = RGCRFweights{nearbyRGCindices(kk)};

        % The current cost
        theOriginalCost = RGCRFconnector.costToMaintainInputCones( ...
            wiringParams.chromaticSpatialVarianceTradeoff, ...
            allConePositions(theExistingConeInputs,:), ...
            allConeSpacings(theExistingConeInputs), ...
            allConeTypes(theExistingConeInputs), ...
            theExistingConeWeights, ...
            localConeToRGCDensityRatio);

        % Add the coneIndexToTransfer
        theExistingConeInputs(numel(theExistingConeInputs)+1) = coneIndexToTransfer;
        theExistingConeWeights(numel(theExistingConeWeights)+1) = 1;

        % Compute projected cost
        theProjectedCost = RGCRFconnector.costToMaintainInputCones( ...
            wiringParams.chromaticSpatialVarianceTradeoff, ...
            allConePositions(theExistingConeInputs,:), ...
            allConeSpacings(theExistingConeInputs), ...
            allConeTypes(theExistingConeInputs), ...
            theExistingConeWeights, ...
            localConeToRGCDensityRatio);

        changeInCost(kk) = theOriginalCost - theProjectedCost;
    end % for kk

    % Find the RGC that has the minimum change in cost
    [~,idx] = min(changeInCost);
    destinationRGCindex = nearbyRGCindices(idx);

    % Do the transfer
    [RGCRFinputs{theSourceRGCindex}, RGCRFweights{theSourceRGCindex}, ...
        RGCRFinputs{destinationRGCindex}, RGCRFweights{destinationRGCindex}, ...
        RGCRFcentroids(theSourceRGCindex,:), coneInputsNum(theSourceRGCindex), ...
        RGCRFcentroids(destinationRGCindex,:), coneInputsNum(destinationRGCindex)] = ...
        RGCRFconnector.transferConeFromSourceRGCtoDestinationRGC(...
            RGCRFinputs{theSourceRGCindex}, RGCRFweights{theSourceRGCindex}, ...
            coneIndexToTransfer, 1.0, ...
            RGCRFinputs{destinationRGCindex}, RGCRFweights{destinationRGCindex}, ...
            allConePositions);

    % Compute new RGCRFspacings
    RGCRFspacings = RGCmodels.Watson.convert.positionsToSpacings(RGCRFcentroids);

end