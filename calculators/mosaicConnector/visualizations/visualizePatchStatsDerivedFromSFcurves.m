function visualizePatchStatsDerivedFromSFcurves(patchDogModelParams, patchRGCeccentricityDegs)
    size(patchDogModelParams)
    size(patchRGCeccentricityDegs)
    rgcsNum = numel(patchDogModelParams);
    
    centerRadii = zeros(1, rgcsNum);
    surroundRadii = zeros(1, rgcsNum);
    centerPeakSensitivities = zeros(1, rgcsNum);
    surroundPeakSensitivities = zeros(1, rgcsNum);
    
    for iRGC = 1:rgcsNum
        p = patchDogModelParams{iRGC};
        centerRadii(iRGC) = p.rC;
        surroundRadii(iRGC) = p.rS;
        centerPeakSensitivities(iRGC) = p.kC;
        surroundPeakSensitivities(iRGC) = p.kS;
    end
    
    figure(555); 
    subplot(2,2,1);
    plot(patchRGCeccentricityDegs, centerRadii, 'ro', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 12); hold on
    plot(patchRGCeccentricityDegs, surroundRadii, 'bo', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerSize', 12);
    xlabel('eccentricity (degs)');
    ylabel('radius (degs)');
    legend({'center', 'surround'});
    set(gca, 'XScale', 'log', 'XLim', [0.01 30], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(gca, 'YScale', 'log', 'YLim', [0.01 3], 'YTick', [0.01 0.03 0.1 0.3 1 3]);
    axis 'square';
    
    subplot(2,2,2);
    plot(centerRadii,  centerPeakSensitivities, 'ro', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 12); hold on
    plot(surroundRadii, surroundPeakSensitivities, 'bo', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerSize', 12);
    ylabel('peak sensitivity');
    xlabel('radius (degs)');
    title('surround');
    set(gca, 'XScale', 'log', 'XLim', [0.01 10], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(gca, 'YScale', 'log', 'YLim', [1 100000], 'YTick', [1 3 10 30 100 300 1000 3000 10000 30000 100000]);
    axis 'square';

    subplot(2,2,3);
    plot(patchRGCeccentricityDegs, surroundPeakSensitivities./centerPeakSensitivities, 'ko', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 12); hold on
    xlabel('eccentricity (degs)');
    ylabel('peak sensitivity (surround/center)');
    legend({'center', 'surround'});
    set(gca, 'XScale', 'log', 'XLim', [0.01 30], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(gca, 'YScale', 'log', 'YLim', [1e-3 1], 'YTick', [0.0001 0.001 0.01  0.1  1]);
    axis 'square';
    
    subplot(2,2,4);
    plot(patchRGCeccentricityDegs, (surroundRadii/centerRadii).^2 .* (surroundPeakSensitivities./centerPeakSensitivities), 'ko', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 12); hold on
    xlabel('eccentricity (degs)');
    ylabel('integrated sensitivity (surround/center)');
    legend({'center', 'surround'});
    set(gca, 'XScale', 'log', 'XLim', [0.01 30], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(gca, 'YScale', 'linear', 'YLim', [0 1], 'YTick', [0:0.2:1]);
    axis 'square';
    pause
end


