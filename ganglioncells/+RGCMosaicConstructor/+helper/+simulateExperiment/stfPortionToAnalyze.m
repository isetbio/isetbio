function [sfSupportForAnalysis, theSTFforAnalysis] = stfPortionToAnalyze(sfSupport, theSTF, ...
    deltaThresholdToDeclareLocalMinInSTF, minPeakProminance)

    if (~isempty(minPeakProminance))
        % Find the peaks that drop at least 0.01 on either side before the signal attains a higher value.
        minPeakProminance = 1e-2;
        [~,indicesOfPeaks] = findpeaks(theSTF, 'MinPeakProminence', minPeakProminance);
        if (~isempty(indicesOfPeaks))
            % Only consider the first peak
            iSF = indicesOfPeaks(1);
            maxSTF = theSTF(iSF);
        else
            % findpeaks did not find any peaks so go with the max
            [maxSTF,iSF] = max(theSTF(:));
        end
    else
        % max of the STF
        [maxSTF,iSF] = max(theSTF(:));
    end

    % Start at the peakSF and march upwards to find the first local min
    stillFalling = true;
    deltaThresholdToDeclareLocalMinInSTF = maxSTF * deltaThresholdToDeclareLocalMinInSTF;

    localMinSF = numel(theSTF);
    thresholdDelta = deltaThresholdToDeclareLocalMinInSTF*max(theSTF(:));

    while (iSF < numel(theSTF))&&(stillFalling)
        iSF = iSF + 1;
        delta = theSTF(iSF) - theSTF(iSF-1);
        if (delta > thresholdDelta)
            localMinSF = iSF-1;
            stillFalling = false;
        end
    end

    sfSupportForAnalysis = sfSupport(1:localMinSF);
    theSTFforAnalysis = theSTF(1:localMinSF);

    visualize = ~true;
    if (visualize)
        figure(22); clf
        plot(sfSupport, theSTF, 'ks-')
        set(gca, 'XScale', 'log', 'XLim', [0.1 100], 'YLim', [0 1]);
        hold on
        plot(sfSupportForAnalysis, theSTFforAnalysis, 'ro-');
    end
end