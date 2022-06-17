function [status, RGCRFinputs, RGCRFweights, RGCRFcentroids, RGCRFspacings, coneInputsNum] = tryToReassignConeInput( ...
                        targetSourceRGCindex, nearbyRGCindices,  ...
                        RGCRFinputs, RGCRFweights, RGCRFcentroids, RGCRFspacings, coneInputsNum, ...
                        localConeToRGCDensityRatio, ...
                        wiringParams, allConePositions, allConeSpacings, allConeTypes)

    % Retrieve the target source RGC cone inputs
    theSourceRGCinputConeIndices = RGCRFinputs{targetSourceRGCindex};
    theSourceRGCinputConeWeights = RGCRFweights{targetSourceRGCindex};
    sourceRGCconeInputsNum = numel(theSourceRGCinputConeIndices);
     
    for iNearbyRGC = 1:numel(nearbyRGCindices)
        theCandidateNearbyRGCindex = nearbyRGCindices(iNearbyRGC);

        % How many inputs does the candidate nearby RGC has?
        theCandidateRGCinputConeIndices = RGCRFinputs{theCandidateNearbyRGCindex};
        theCandidateRGCinputConeWeights = RGCRFweights{theCandidateNearbyRGCindex};
        candidateRGCconeInputsNum = numel(theCandidateRGCinputConeIndices);
        
        if (candidateRGCconeInputsNum < sourceRGCconeInputsNum-1)
            % eg: source RGC has N cones and destination RGC has at most N-2 cones, so
            % after the transfer the source RGC will have N-1 cones and the destination 
            % RGC will have at most N-1 cones.
            % Check the benefits of TRANSFERING one cone from the source RGC to
            % any of the nearby RGCs, and find the most beneficial transfer
            % params,i.e.: (sourceConeIndex, destinationNearbyRGCindex)

            [actionBenefit(iNearbyRGC), actionParams{iNearbyRGC}] = ...
                RGCRFconnector.computeCostOfConeTransfer(...
                    theSourceRGCinputConeIndices, ...
                    theSourceRGCinputConeWeights, ...
                    theCandidateRGCinputConeIndices, ...
                    theCandidateRGCinputConeWeights, ...
                    localConeToRGCDensityRatio(targetSourceRGCindex), ...
                    wiringParams, ...
                    allConePositions, allConeSpacings, allConeTypes);

        elseif (candidateRGCconeInputsNum <= sourceRGCconeInputsNum)
           fprintf('Source RGC has %d cones, candidate destination RGC (%d) has %d cones. WAITING TO IMPLEMENT CONE EXCHANGE.\n', ...
                sourceRGCconeInputsNum, iNearbyRGC, candidateRGCconeInputsNum);
            % eg: source RGC has N cone and destination RGC has N-1 or N cones.
            % check the benefits of EXCHANGING one cone from the source RGC
            % with a cone of any of the nearby RGCs, and find
            % the most beneficial exchange params, i.e: (sourceConeIndex,
            % destinationNearbyRGCindex, destinationNearbyRGCconeInputIndex)

            %[actionBenefit(iNearbyRGC), actionParams{iNearbyRGC}] = ...
            %    RGCRFconnector.computeCostOfConeExcange();

            actionBenefit(iNearbyRGC) = -inf;
            actionParams{iNearbyRGC} = struct();

        else
            % source RGC has N cones and destination RGC has N+1 or more
            % cones. Do not do anything
            actionBenefit(iNearbyRGC) = -inf;
            actionParams{iNearbyRGC} = struct();
        end

    end %iNearbyRGC

    % Find the max benefit
    [maxBenefit,iNearbyRGCbest] = max(actionBenefit);
    
    if (maxBenefit < 0)
        status = 'failed';
        return;
    end
    
    % Select the best actionParams
    status = 'success';
    actionParams = actionParams{iNearbyRGCbest};
    
    switch (actionParams.type)
        case 'transfer cone'
            % Do the cone transfer. 

            destinationRGCindex = nearbyRGCindices(iNearbyRGCbest);
            sourceRGCconeIndex = actionParams.sourceRGCconeIndex;

            [RGCRFinputs{targetSourceRGCindex}, RGCRFweights{targetSourceRGCindex}, ...
             RGCRFinputs{destinationRGCindex}, RGCRFweights{destinationRGCindex}, ...
             RGCRFcentroids(targetSourceRGCindex,:), coneInputsNum(targetSourceRGCindex), ...
             RGCRFcentroids(destinationRGCindex,:), coneInputsNum(destinationRGCindex)] = ...
                    RGCRFconnector.transferConeFromSourceRGCtoDestinationRGC(...
                        RGCRFinputs{targetSourceRGCindex}, RGCRFweights{targetSourceRGCindex}, ...
                        sourceRGCconeIndex, 1.0, ...
                        RGCRFinputs{destinationRGCindex}, RGCRFweights{destinationRGCindex}, ...
                        allConePositions);

            % Compute new RGCRFspacings
            RGCRFspacings = RGCmodels.Watson.convert.positionsToSpacings(RGCRFcentroids);

        case 'exchange cone'
            % Do the cone exchange

            destinationRGCindex = nearbyRGCindices(iNearbyRGCbest);
            sourceRGCconeIndex = actionParams.sourceRGCconeIndex;
            destinationRGCconeIndex = actionParams.destinationRGCconeIndex;

            % Source RGC ----> destination RGC
            [RGCRFinputs{targetSourceRGCindex}, RGCRFweights{targetSourceRGCindex}, ...
             RGCRFinputs{destinationRGCindex}, RGCRFweights{destinationRGCindex}, ...
             RGCRFcentroids(targetSourceRGCindex,:), coneInputsNum(targetSourceRGCindex), ...
             RGCRFcentroids(destinationRGCindex,:), coneInputsNum(destinationRGCindex)] = ...
                    RGCRFconnector.transferConeFromSourceRGCtoDestinationRGC(...
                        RGCRFinputs{targetSourceRGCindex}, RGCRFweights{targetSourceRGCindex}, ...
                        sourceRGCconeIndex, 1.0, ...
                        RGCRFinputs{destinationRGCindex}, RGCRFweights{destinationRGCindex}, ...
                        allConePositions);

            % Source RGC <--- destination RGC
            [RGCRFinputs{destinationRGCindex}, RGCRFweights{destinationRGCindex}, ...
             RGCRFinputs{targetSourceRGCindex}, RGCRFweights{targetSourceRGCindex}, ...
             RGCRFcentroids(destinationRGCindex,:), coneInputsNum(destinationRGCindex), ...
             RGCRFcentroids(targetSourceRGCindex,:), coneInputsNum(targetSourceRGCindex)] = ...
                    RGCRFconnector.transferConeFromSourceRGCtoDestinationRGC(...
                        RGCRFinputs{destinationSourceRGCindex}, RGCRFweights{destinationSourceRGCindex}, ...
                        destinationRGCconeIndex, 1.0, ...
                        RGCRFinputs{sourceRGCindex}, RGCRFweights{sourceRGCindex}, ...
                        allConePositions);

            % Compute new RGCRFspacings
            RGCRFspacings = RGCmodels.Watson.convert.positionsToSpacings(RGCRFcentroids);

        otherwise
            error('Unknown action.type: ''%s''.', actionParams.type);
    end

    

end