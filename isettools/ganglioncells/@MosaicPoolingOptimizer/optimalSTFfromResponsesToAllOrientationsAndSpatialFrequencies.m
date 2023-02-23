function theOptimalSTF = optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies(obj, ...
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
    theOptimalSTF = obj.highestExtensionSTF(theSTFsAcrossAllOrientations);
end