function [theOptimalNormalizedSTF, theNormalizedSTFsAcrossAllOrientations, ...
    theHighestExtensionOrientation, theOptimalSTFMagnitudeSpectrum, theOptimalSTFphaseSpectrum] = optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies(...
        orientationsTested, spatialFrequenciesTested, ...
        theResponsesAcrossAllOrientationsAndSpatialFrequencies)

    orientationsNum = size(theResponsesAcrossAllOrientationsAndSpatialFrequencies,1);
    sfsNum = size(theResponsesAcrossAllOrientationsAndSpatialFrequencies,2);

    % Compute modelRGC STFs for all orientations
    theSTFmagnitudeSpectraAcrossAllOrientations = zeros(orientationsNum, sfsNum);
    theSTFPhaseSpectraAcrossAllOrientations = zeros(orientationsNum, sfsNum);

    timeSupport = 0:(size(theResponsesAcrossAllOrientationsAndSpatialFrequencies,3)-1);
    timeSupport = timeSupport/(numel(timeSupport)-1);
    timeSupportHR = timeSupport(1):0.01:timeSupport(end);

    temporalFrequency = 1.0;

    parfor iOri = 1:orientationsNum
        for iSF = 1:sfsNum
            % Retrieve cone mosaic responses for all frames of this stimulus
            theResponseTimeSeries = squeeze(theResponsesAcrossAllOrientationsAndSpatialFrequencies(iOri, iSF,:));

            % Fit a sinusoid to extract the magnitude and phase
            [~, fittedParams] = MosaicPoolingOptimizer.fitSinusoidToResponseTimeSeries(...
                timeSupport, theResponseTimeSeries, temporalFrequency, timeSupportHR);
            theSTFmagnitudeSpectraAcrossAllOrientations(iOri, iSF) = fittedParams(1);
            theSTFPhaseSpectraAcrossAllOrientations(iOri, iSF) = fittedParams(2);

        end
    end

    % Pick the highest extension STF as the visual STF for this cell
    [theOptimalNormalizedSTF,theNormalizedSTFsAcrossAllOrientations, theHighestExtensionOrientation, ...
        theOptimalSTFMagnitudeSpectrum, theOptimalSTFphaseSpectrum] = highestExtensionSTF(orientationsTested, spatialFrequenciesTested, ...
        theSTFmagnitudeSpectraAcrossAllOrientations, theSTFPhaseSpectraAcrossAllOrientations);
end

        
function [theHighestExtensionSTF, theNormalizedSTFs, theHighestExtensionOrientation, ...
    theOptimalSTFMagnitudeSpectrum, theOptimalSTFphaseSpectrum] = highestExtensionSTF(...
        orientationsTested, spatialFrequenciesTested, theMeasuredSTFs, theMeasuredSTFphases)

    % Normalize to max
    theNormalizedSTFs = theMeasuredSTFs / max(theMeasuredSTFs(:));

    % Determine the orientation that maximizes the STF extension to high spatial frequencies
    maxSF = nan(1,numel(orientationsTested));
    theSTFtoAnalyze = cell(1, numel(orientationsTested));

    for iOri = 1:numel(orientationsTested)
        % Find spatial frequency at which STF drops to  max x
        % MosaicPoolingOptimizer.highSFAttenuationFactorForOptimalOrientation
        theSTFatThisOrientation = squeeze(theNormalizedSTFs(iOri,:));
        
        % Only fit the part of the STF that contains the main peak
        [spatialFrequenciesToAnalyze, theSTFtoAnalyze{iOri}] = ...
            MosaicPoolingOptimizer.stfPortionToAnalyze(spatialFrequenciesTested, theSTFatThisOrientation);

        maxSF(iOri) = max(spatialFrequenciesToAnalyze);
    end

    idx = find(maxSF == max(maxSF));
    STFmagAtHighestSF = zeros(1, numel(idx));
    for i = 1:numel(idx)
        iOri = idx(i);
        theSTF = theSTFtoAnalyze{iOri};
        STFmagAtHighestSF(i) = theSTF(end);
    end

    [~,idx2] = max(STFmagAtHighestSF);

    bestOri = idx(idx2);
    theHighestExtensionSTF = squeeze(theNormalizedSTFs(bestOri,:));
    theHighestExtensionOrientation = orientationsTested(bestOri);

    theOptimalSTFMagnitudeSpectrum = theMeasuredSTFs(bestOri,:);
    theOptimalSTFphaseSpectrum = theMeasuredSTFphases(bestOri,:);
end


