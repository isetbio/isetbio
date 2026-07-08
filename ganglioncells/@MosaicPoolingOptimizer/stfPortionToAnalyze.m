function [sfSupportForAnalysis, theSTFforAnalysis] = stfPortionToAnalyze(sfSupport, theSTF)

    % find the peak of the STF
    [maxSTF,idx] = max(theSTF(:));

    % Start at the peakSF and march upwards to find the first local min
    iSF = idx;
    stillFalling = true;
    deltaThresholdToDeclareLocalMinInSTF = maxSTF * MosaicPoolingOptimizer.deltaThresholdToDeclareLocalMinInSTF;

    localMinSF = numel(theSTF);
    while (iSF < numel(theSTF))&&(stillFalling)
        iSF = iSF + 1;
        delta = theSTF(iSF) - theSTF(iSF-1);
        if (delta > deltaThresholdToDeclareLocalMinInSTF)
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
        fprintf('Local min at %d\n', sfSupport(localMinSF));
        hold on
        plot(sfSupportForAnalysis, theSTFforAnalysis, 'ro-');
        pause
    end

end
