function theSTFdata = currentSTF(pooledConeIndicesAndWeights, inputConeMosaicSTFresponses, ...
    visualRcDegs, fixedOptimalOrientation, targetVisualSTFparams, targetVisualSTFparamFractionalTolerances,  multiStartsNumDoGFit,  deltaThresholdToDeclareLocalMinInSTF, stimParams,axDoGFitToCompositeSTFparams)

    theOptimalOrientationIndex = find(stimParams.orientationDegs == fixedOptimalOrientation);

    % Compute responses of all center cones
    conesNum = numel(pooledConeIndicesAndWeights.centerConeIndices);    
    centerConeWeights = reshape(pooledConeIndicesAndWeights.centerConeWeights(:), [1 1 conesNum]);
    centerConeResponses = squeeze(inputConeMosaicSTFresponses(theOptimalOrientationIndex,:,:, pooledConeIndicesAndWeights.centerConeIndices));
    netCenterResponses = sum(bsxfun(@times, centerConeResponses, centerConeWeights), 3);
    

    % Compute responses of all surronud cones
    conesNum = numel(pooledConeIndicesAndWeights.surroundConeWeights);
    surroundConeWeights = reshape(pooledConeIndicesAndWeights.surroundConeWeights(:), [1 1 conesNum]);
    surroundConeResponses = squeeze(inputConeMosaicSTFresponses(theOptimalOrientationIndex,:,:, pooledConeIndicesAndWeights.surroundConeIndices));
    netSurroundResponses = sum(bsxfun(@times, surroundConeResponses, surroundConeWeights), 3);

   	% Composite responses (center-surround) of the current modelRGC (defined by the current pooledConeIndicesAndWeights)
    theCompositeSTFresponses = netCenterResponses - netSurroundResponses;
        
    theSTIMamplitudeSpectrum = RGCMosaicConstructor.helper.simulateExperiment.stfFromResponseTimeSeries(...
            stimParams.spatialFrequencyCPD, theCompositeSTFresponses, ...
            stimParams.spatialPhasesDegs, stimParams.temporalSupportSeconds);
        
	% Fit the visual STF with a DoG model
    [theFittedDoGModelParams, theFittedDoGModelToTheVisualSTF] = ...
        RGCMosaicConstructor.helper.fit.DifferenceOfGaussiansToCompositeSTF(...
            stimParams.spatialFrequencyCPD, theSTIMamplitudeSpectrum, ...
            visualRcDegs, targetVisualSTFparams, targetVisualSTFparamFractionalTolerances, ...
            multiStartsNumDoGFit, deltaThresholdToDeclareLocalMinInSTF, ...
            axDoGFitToCompositeSTFparams);

    theSTFdata = struct;
    theSTFdata.spatialFrequencySupport = stimParams.spatialFrequencyCPD;
    theSTFdata.visualSTF = theSTIMamplitudeSpectrum;
    theSTFdata.fittedDoGModelParams = theFittedDoGModelParams;
    theSTFdata.fittedDoGModelToVisualSTF = theFittedDoGModelToTheVisualSTF;
    theSTFdata.residual = theFittedDoGModelToTheVisualSTF.RMSE;

    if (1==2)
        % Penalize strongly if the surround becomes stronger than the center at the lowest SF
        lowestSFcenterResponseTimeSeries = netCenterResponses(1,:);
        lowestSFsurroundResponseTimeSeries = netSurroundResponses(1,:);
        [maxCenter,idx] = max(lowestSFcenterResponseTimeSeries);
        maxSurround = (lowestSFsurroundResponseTimeSeries(idx));
        if (maxSurround > maxCenter)
            amplifier = 10.0*double(1+(maxSurround - 2*maxCenter)/maxCenter);
            theSTFdata.residual = theSTFdata.residual*amplifier;
        end
    end

    theSTFdata.fittedRcDegs = theFittedDoGModelParams.finalValues(4);
    theSTFdata.fittedDoGModelRsRcRatio = theFittedDoGModelParams.finalValues(3);
    theSTFdata.fittedDoGModelSCIntSensRatio = theFittedDoGModelParams.finalValues(2);

    theSTFdata.spatialFrequencySupportHR = theFittedDoGModelToTheVisualSTF.sfHiRes;
    theSTFdata.visualCompositeSTFHR = theFittedDoGModelToTheVisualSTF.compositeSTFHiRes;
    theSTFdata.visualCenterSTFHR = theFittedDoGModelToTheVisualSTF.centerSTFHiRes;
    theSTFdata.visualSurroundSTFHR = theFittedDoGModelToTheVisualSTF.surroundSTFHiRes;

end
