function [theOptimalSTF, theSTFsAcrossAllOrientations, theHighestExtensionOrientation] = optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies(...
    orientationsTested, spatialFrequenciesTested, ...
    theResponsesAcrossAllOrientationsAndSpatialFrequencies)

    orientationsNum = size(theResponsesAcrossAllOrientationsAndSpatialFrequencies,1);
    sfsNum = size(theResponsesAcrossAllOrientationsAndSpatialFrequencies,2);

    % Compute modelRGC STFs for all orientations
    theSTFsAcrossAllOrientations = zeros(orientationsNum, sfsNum);

    for iOri = 1:orientationsNum
        for iSF = 1:sfsNum
            % Retrieve cone mosaic responses for all frames of this stimulus
            theResponseTimeSeries = squeeze(theResponsesAcrossAllOrientationsAndSpatialFrequencies(iOri, iSF,:));
            % STF as the amplitude of modulation of the response (max-min)
            theSTFsAcrossAllOrientations(iOri, iSF) = max(theResponseTimeSeries(:)) - min(theResponseTimeSeries(:));
        end
    end

    % Pick the highest extension STF as the visual STF for this cell
    [theOptimalSTF,theSTFsAcrossAllOrientations, theHighestExtensionOrientation] = ...
        highestExtensionSTF(orientationsTested, spatialFrequenciesTested, ...
        theSTFsAcrossAllOrientations);
end

        
function [theHighestExtensionSTF, theMeasuredSTFs, theHighestExtensionOrientation] = highestExtensionSTF(...
    orientationsTested, spatialFrequenciesTested, theMeasuredSTFs)

    theMeasuredSTFs = theMeasuredSTFs / max(theMeasuredSTFs(:));

    % Determine the orientation that maximizes the STF extension to high spatial frequencies
    maxSF = nan(1,numel(orientationsTested));
    
    for iOri = 1:numel(orientationsTested)
        % Find spatial frequency at which STF drops to  max x
        % MosaicPoolingOptimizer.highSFAttenuationFactorForOptimalOrientation
        theSTFatThisOrientation = squeeze(theMeasuredSTFs(iOri,:));
        
        % Only fit the part of the STF that contains the main peak
        [spatialFrequenciesToAnalyze, theSTFtoAnalyze{iOri}] = ...
            MosaicPoolingOptimizer.stfPortionToAnalyze(spatialFrequenciesTested, theSTFatThisOrientation);

        maxSF(iOri) = max(spatialFrequenciesToAnalyze);
    end

    idx = find(maxSF == max(maxSF));
    for i = 1:numel(idx)
         theSTF = theSTFtoAnalyze{iOri};
         maxSTF(i) = max(theSTF);
    end

    [~,idx2] = max(maxSTF);

    bestOri = idx(idx2);
    theHighestExtensionSTF = squeeze(theMeasuredSTFs(bestOri,:));
    theHighestExtensionOrientation = orientationsTested(bestOri);
end


