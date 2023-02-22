function theSTFdata = visualSTFfromCronerKaplanAnalysisOfconeMosaicSTFresponses(obj, pooledConeIndicesAndWeights)

    % Compute the pooled center signals for all sfs, all orientations, all
    % spatial phases
    conesNum = numel(pooledConeIndicesAndWeights.centerConeWeights);
    centerConeWeights = reshape(pooledConeIndicesAndWeights.centerConeWeights(:), [1 1 1 conesNum]);
    centerConeResponses = obj.stfComputeMethodResources.theConeMosaicResponses(:,:,:, pooledConeIndicesAndWeights.centerConeIndices);
    netCenterResponses = sum(bsxfun(@times, centerConeResponses, centerConeWeights), 4);

    % Compute the pooled surround signals
    conesNum = numel(pooledConeIndicesAndWeights.surroundConeWeights);
    surroundConeWeights = reshape(pooledConeIndicesAndWeights.surroundConeWeights(:), [1 1 1 conesNum]);
    surroundConeResponses = obj.stfComputeMethodResources.theConeMosaicResponses(:,:,:, pooledConeIndicesAndWeights.surroundConeIndices);
    netSurroundResponses = sum(bsxfun(@times, surroundConeResponses, surroundConeWeights), 4);

    % Composite (center-surround) signals
    compositeResponses = netCenterResponses - netSurroundResponses;

    % Compute STFs for all orientations
    theMeasuredSTFs = zeros(...
        numel(obj.stfComputeMethodResources.orientationsTested), ...
        numel(obj.stfComputeMethodResources.spatialFrequenciesTested) ...
        );

    for iOri = 1:numel(obj.stfComputeMethodResources.orientationsTested)
        for iSF = 1:numel(obj.stfComputeMethodResources.spatialFrequenciesTested)
            % Retrieve cone mosaic responses for all frames of this stimulus
            theRGCresponseTimeSeries = squeeze(compositeResponses(iOri, iSF,:));
            theMeasuredSTFs(iOri, iSF) = max(theRGCresponseTimeSeries) - min(theRGCresponseTimeSeries);
        end
    end

    % Pick the highest extension STF as the visual STF for this cell
    theVisualSTF = highestExtensionSTF(obj.stfComputeMethodResources.orientationsTested, obj.stfComputeMethodResources.spatialFrequenciesTested, theMeasuredSTFs);

    % Fit the visual STF with a DoG model
    [theFittedDoGModelParams, theFittedDoGModelToTheVisualSTF] = RTVF.fitDoGmodelToMeasuredSTF(...
                      obj.stfComputeMethodResources.spatialFrequenciesTested, ...
                      theVisualSTF, ...
                      obj.visualRFcenterRcDegs, ...
                      obj.visualRFcenterRcDegs*[0.5 1 1.5], ...
                      obj.multiStartsNumDoGFit);

    theSTFdata = struct;
    theSTFdata.spatialFrequencySupport = obj.stfComputeMethodResources.spatialFrequenciesTested;
    theSTFdata.visualSTF = theVisualSTF;
    theSTFdata.fittedDoGModelParams = theFittedDoGModelParams;
    theSTFdata.fittedDoGModelToVisualSTF = theFittedDoGModelToTheVisualSTF;
    theSTFdata.fittedRcDegs = theFittedDoGModelParams.finalValues(4);
    theSTFdata.fittedDoGModelRsRcRatio = theFittedDoGModelParams.finalValues(3);
    theSTFdata.fittedDoGModelSCIntSensRatio = theFittedDoGModelParams.finalValues(2) * (theFittedDoGModelParams.finalValues(3))^2;
    
end


function theHighestExtensionSTF = highestExtensionSTF(orientationsTested, spatialFrequenciesTested, theMeasuredSTFs)

    theMeasuredSTFs = theMeasuredSTFs / max(theMeasuredSTFs(:));

    % Determine the orientation that maximizes the STF extension to high spatial frequencies
    maxSF = nan(1,numel(orientationsTested));
    spatialFrequenciesInterpolated = linspace(spatialFrequenciesTested(1),spatialFrequenciesTested(end), 50);

    for iOri = 1:numel(orientationsTested)
        % Find spatial frequency at which STF drops to 20% of max
        theSTFatThisOri = squeeze(theMeasuredSTFs(iOri,:));
        theSTFatThisOriInterpolated = interp1(spatialFrequenciesTested, theSTFatThisOri, spatialFrequenciesInterpolated);
        [mag, iSFpeak] = max(theSTFatThisOri);
        thresholdSTF = mag * 0.2;

        ii = iSFpeak;
        keepGoing = true; iStop = [];
        while (ii < numel(spatialFrequenciesInterpolated)-1)&&(keepGoing)
            ii = ii + 1;
            if (theSTFatThisOriInterpolated(ii)>=thresholdSTF) && (theSTFatThisOriInterpolated(ii+1)<thresholdSTF)
                keepGoing = false;
                iStop = ii;
            end
        end % while
        if (~isempty(iStop))
            maxSF(iOri) = spatialFrequenciesInterpolated(iStop);
        end
    end % iOri

    % Best orientation
    if (any(isnan(maxSF)))
        theSTFatTheHighestSF = squeeze(theMeasuredSTFs(:,end));
        [~, iBestOri] = max(theSTFatTheHighestSF(:));
    else
        [~, iBestOri] = max(maxSF);
    end
    theHighestExtensionSTF = squeeze(theMeasuredSTFs(iBestOri,:));
end

