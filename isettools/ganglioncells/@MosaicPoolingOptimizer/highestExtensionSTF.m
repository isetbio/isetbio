function theHighestExtensionSTF = highestExtensionSTF(obj,theMeasuredSTFs)

    orientationsTested = obj.inputConeMosaicVisualSTFdata.orientationsTested;
    spatialFrequenciesTested = obj.inputConeMosaicVisualSTFdata.spatialFrequenciesTested;

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