function [actionBenefit, actionParamsStruct] = computeCostOfConeTransfer(...
    theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
    theDestinationRGCinputConeIndices, theDestinationRGCinputConeWeights, ...
    localConeToRGCDensityRatio, ...
    wiringParams, allConePositions, allConeSpacings, allConeTypes)
    
    % Original cost for source RGC to retain the input cones
    originalCostInSourceRGC = RGCRFconnector.costToMaintainInputCones( ...
                    wiringParams.chromaticSpatialVarianceTradeoff, ...
                    allConePositions(theSourceRGCinputConeIndices,:), ...
                    allConeSpacings(theSourceRGCinputConeIndices), ...
                    allConeTypes(theSourceRGCinputConeIndices), ...
                    theSourceRGCinputConeWeights, ...
                    localConeToRGCDensityRatio);

    % Compute projected cost for source RGC when transfering each of its
    % cone inputs 
    for iCone = 1:numel(theSourceRGCinputConeIndices)
        % The index of the cone to be transfered out of the source RGC
        theCandidateSourceRGCConeIndex = theSourceRGCinputConeIndices(iCone);
        
        % Remaining pool of cone inputs to the source RGC after the cone transfer
        [newSourceRGCinputConeIndices,ia] = setdiff(theSourceRGCinputConeIndices, theCandidateSourceRGCConeIndex);
        newSourceRGCinputConeWeights = theSourceRGCinputConeWeights(ia);

        % Projected cost in source RGC after the transfer.
        projectedCostInSourceRGC(iCone) = RGCRFconnector.costToMaintainInputCones( ...
            wiringParams.chromaticSpatialVarianceTradeoff, ...
            allConePositions(newSourceRGCinputConeIndices,:), ...
            allConeSpacings(newSourceRGCinputConeIndices), ...
            allConeTypes(newSourceRGCinputConeIndices), ...
            newSourceRGCinputConeWeights, ...
            localConeToRGCDensityRatio);

    end % iCone
   
    % Find the iCone that maximizes the benefit for the sourceRGC
    actionBenefitSourceRGCallCones = originalCostInSourceRGC - projectedCostInSourceRGC;
    actionBenefitSourceRGCMax = max(actionBenefitSourceRGCallCones);

    % If the maximal benefit is < 0, abort
    if (actionBenefitSourceRGCMax < 0)
        actionBenefit = -inf;
        actionParamsStruct = struct();
        return;
    end

    % Fraction used to determine if some cone inputs are near enough to the maxBenefit
    thresholdBenefitFraction = 0.3;

    % Check to see if there are more than one iCones that are close enough to the maximum sourceRGC benefit,
    thresholdBenefit = thresholdBenefitFraction*actionBenefitSourceRGCMax;
    [~, iConesBest] = find(abs(actionBenefitSourceRGCallCones - actionBenefitSourceRGCMax) <= thresholdBenefit);

    if (numel(iConesBest)==1)
        %fprintf('best benefit for a single cone: %d\n', iConesBest);
    else
        if (wiringParams.chromaticSpatialVarianceTradeoff >= 0.5)
            % Choose the cone that lies closest to the cones of the nearby RGC
            [~,idx] = RGCRFconnector.closestConePositionFromSourcePoolToDestinationPoolCones(...
                theSourceRGCinputConeIndices(iConesBest), theDestinationRGCinputConeIndices, ...
                allConePositions);
        else
            % Choose the cone whos type matches the majority of the cone
            % types of the nearby RGC
            [~,idx] = RGCRFconnector.closestConeTypeFromSourcePoolToDestinationPoolCones(...
                theSourceRGCinputConeIndices(iConesBest), theDestinationRGCinputConeIndices, ...
                allConeTypes, allConePositions);
        end

        iConesBest = iConesBest(idx);
    end

    % Output
    actionBenefit = actionBenefitSourceRGCallCones(iConesBest);

    % Form the actionParams struct
    actionParamsStruct = struct(...
        'type', 'transfer cone', ...
        'sourceRGCconeIndex', theSourceRGCinputConeIndices(iConesBest) ...
        );

end